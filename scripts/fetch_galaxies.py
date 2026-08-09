#!/usr/bin/env python
"""Fetch DR9 brick coadds for SGA-2020 galaxies and write training samples.

Run scripts/fetch_catalogues.py once first. Selection (redshift cut,
footprint test, brick resolution, single-galaxy isolation) is
deterministic and recomputed per run; the two heavy stages resume from
what is already on disk, so a re-run never re-fetches:

    fetch  mirrors every needed brick's coadd files (image, invvar,
           maskbits, psfsize per band, plus the tractor catalogue) from
           the static file server, sha256-verified
    cut    writes one SCI/IVAR/MASK/LAYERS sample per galaxy under
           <root>/samples/galaxies/ across --jobs worker processes
           (fz decompression holds the GIL, so threads cannot help),
           appends a manifest row per target, and reports measured
           volume plus the extrapolation to the full selection

The manifest is the resume ledger: a target with a manifest row is
skipped, anything else is (re)cut — samples are written atomically, so
a row implies a whole file. Frames are never rejected for content: the
only structural criterion is all-band coverage below the threshold
shared with the background driver.

Example:
    python scripts/fetch_galaxies.py --sample 200
"""

import argparse
import csv
import json
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

from astropy.io import fits
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import fetch
from src.bricks import Bricks
from src.catalog import Catalog
from src.cutout import Coadd, cut, write_sample
from src.mask import DEFAULT_COVERAGE_MIN, compute_layers
from src.select import galaxy_targets, subsample

MANIFEST_FIELDS = [
    "file", "sga_id", "galaxy", "ra", "dec", "brick", "hemisphere", "size_px",
    "psf_g", "psf_r", "psf_z", "valid_frac", "bright_frac", "source_frac",
    "galaxy_frac", "galaxy_bit_frac", "status", "reason",
]

# One catalogue per worker process, loaded on its first brick.
_catalog = None


def cut_brick(task):
    """Cut, mask and write every target of one brick; returns manifest rows."""
    global _catalog
    if _catalog is None:
        _catalog = Catalog()
    root, out_dir, brickname, brick_targets, coverage_min = task
    brick = Coadd(fetch.brick_dir(root, brickname))
    with fits.open(fetch.brick_dir(root, brickname)
                   / f"tractor-{brickname}.fits") as hdul:
        tractor = hdul[1].data
    records = []
    for target in brick_targets:
        path = Path(out_dir) / f"sga_{target['sga_id']}.fits"
        sample = cut(brick, target["ra"], target["dec"], target["size_px"])
        record = {
            "file": path.name, "sga_id": target["sga_id"],
            "galaxy": target["galaxy"],
            "ra": f"{target['ra']:.6f}", "dec": f"{target['dec']:.6f}",
            "brick": brickname, "hemisphere": target["hemisphere"],
            "size_px": target["size_px"],
            "psf_g": f"{sample['psf_fwhm'][0]:.4f}",
            "psf_r": f"{sample['psf_fwhm'][1]:.4f}",
            "psf_z": f"{sample['psf_fwhm'][2]:.4f}",
            "valid_frac": f"{sample['valid_frac']:.4f}",
        }
        if sample["valid_frac"] < coverage_min:
            record.update(status="rejected", reason="coverage")
            records.append(record)
            continue
        layers, verdict = compute_layers(sample, _catalog, tractor)
        record.update({key: f"{value:.4f}" for key, value in verdict.items()})
        record.update(status="written", reason="")
        write_sample(path, sample, brickname,
                     {"SGA_ID": (target["sga_id"], "SGA-2020 identifier")},
                     layers)
        records.append(record)
    return records


def report(out_dir, all_targets, bricks, transferred, written):
    """Measured volume this run, extrapolated to the full selection."""
    samples = list(out_dir.glob("sga_*.fits"))
    sample_bytes = sum(p.stat().st_size for p in samples)
    full_bricks = len({t["brick"] for t in all_targets})
    print(f"\nsamples written this run: {written} "
          f"({len(samples)} on disk, {sample_bytes / 1e9:.2f} GB)")
    print(f"bricks fetched this run: {len(bricks)} ({transferred / 1e9:.2f} GB new)")
    if samples and len(bricks):
        per_sample = sample_bytes / len(samples)
        per_brick = transferred / len(bricks) if transferred else 0
        print(f"full selection: {len(all_targets)} galaxies over {full_bricks} bricks")
        print(f"extrapolated full volume: samples "
              f"{per_sample * len(all_targets) / 1e9:.0f} GB"
              + (f", brick mirror {per_brick * full_bricks / 1e12:.2f} TB"
                 if per_brick else ""))


def main(args):
    root = Path(args.root)
    all_targets = galaxy_targets(Catalog(), Bricks(), zcut=args.zcut,
                                 min_nexp=args.min_nexp,
                                 size_mult=args.size_mult,
                                 min_size=args.min_size,
                                 sga_margin=args.sga_margin)
    print(f"{len(all_targets)} isolated galaxies selected")
    targets = subsample(all_targets, args.sample, args.limit)
    print(f"{len(targets)} targets this run")

    bricks = sorted({(t["brick"], t["hemisphere"]) for t in targets})
    transferred = 0
    if args.stage in ("all", "fetch"):
        transferred, failed = fetch.mirror_bricks(bricks, root, args.workers,
                                                  desc=f"{len(bricks)} bricks")
        print(f"fetched {transferred / 1e9:.2f} GB across {len(bricks)} bricks")
        if failed:
            raise SystemExit(f"{len(failed)} bricks failed; re-run to resume")
    if args.stage == "fetch":
        return

    out_dir = root / "samples" / "galaxies"
    out_dir.mkdir(exist_ok=True, parents=True)
    manifest_path = out_dir / "manifest.csv"
    # the manifest is only a valid ledger under the parameters that built it
    spec = {key: getattr(args, key) for key in
            ("zcut", "min_nexp", "size_mult", "min_size", "sga_margin",
             "coverage_min")}
    spec_path = out_dir / "run.json"
    if args.overwrite:
        manifest_path.unlink(missing_ok=True)
        spec_path.unlink(missing_ok=True)
    if spec_path.exists() and json.loads(spec_path.read_text()) != spec:
        raise SystemExit(f"{spec_path} records different science parameters; "
                         "wipe the samples directory or pass --overwrite")
    spec_path.write_text(json.dumps(spec, indent=1))
    done = set()
    if manifest_path.exists():
        with open(manifest_path) as f:
            done = {row["file"] for row in csv.DictReader(f)}
    todo = [t for t in targets if f"sga_{t['sga_id']}.fits" not in done]
    by_brick = {}
    for target in todo:
        by_brick.setdefault(target["brick"], []).append(target)
    tasks = [(str(root), str(out_dir), brickname, brick_targets,
              args.coverage_min)
             for brickname, brick_targets in sorted(by_brick.items())]

    written = 0
    new_manifest = not manifest_path.exists()
    with open(manifest_path, "a", newline="", buffering=1) as f:
        manifest = csv.DictWriter(f, fieldnames=MANIFEST_FIELDS)
        if new_manifest:
            manifest.writeheader()
        failed = 0
        with ProcessPoolExecutor(max_workers=args.jobs) as pool:
            futures = [pool.submit(cut_brick, task) for task in tasks]
            for future in tqdm(as_completed(futures), total=len(futures),
                               desc="cutting", unit="brick"):
                try:
                    records = future.result()
                except Exception as error:
                    failed += 1  # no manifest rows -> retried on the next run
                    tqdm.write(f"brick failed: {error}")
                    continue
                for record in records:
                    manifest.writerow(record)
                    written += record["status"] == "written"
    if failed:
        print(f"{failed} bricks failed; re-run to retry them")

    report(out_dir, all_targets, bricks, transferred, written)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__.split("\n")[0],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--root", default="data/sga",
                        help="data root: brick mirror and samples")
    parser.add_argument("--zcut", type=float, default=0.05,
                        help="upper redshift cut on Z_LEDA")
    parser.add_argument("--min-nexp", type=int, default=1,
                        help="minimum exposures per band for footprint membership")
    parser.add_argument("--size-mult", type=float, default=1.5,
                        help="stamp side as a multiple of D26")
    parser.add_argument("--min-size", type=int, default=64,
                        help="minimum stamp side (pixels)")
    parser.add_argument("--sga-margin", type=float, default=1.5,
                        help="isolation margin on other galaxies' r26")
    parser.add_argument("--coverage-min", type=float, default=DEFAULT_COVERAGE_MIN,
                        help="minimum all-band coverage fraction (structural)")
    parser.add_argument("--sample", type=int, default=None,
                        help="evenly RA-spaced subset of this many targets")
    parser.add_argument("--limit", type=int, default=None,
                        help="first N targets only")
    parser.add_argument("--stage", choices=("all", "fetch", "cut"),
                        default="all", help="run one stage only")
    parser.add_argument("--workers", type=int, default=12,
                        help="parallel download connections")
    parser.add_argument("--jobs", type=int, default=96,
                        help="cut-stage worker processes (measured peak on this "
                             "host; throughput degrades beyond ~96)")
    parser.add_argument("--overwrite", action="store_true",
                        help="reset the sample ledger and re-cut this run's targets")
    main(parser.parse_args())
