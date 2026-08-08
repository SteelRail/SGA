#!/usr/bin/env python
"""Fetch DR9 brick coadds for SGA-2020 galaxies and write training samples.

Run scripts/fetch_catalogues.py once first. Selection (redshift cut,
footprint test, brick resolution) is deterministic and recomputed per
run; the two heavy stages are resumable by file existence, so a re-run
with changed cut parameters never re-fetches. Only single-galaxy frames are selected: a target whose
stamp is reached by any other catalogued galaxy's ellipse is dropped.

    fetch  mirrors every needed brick's coadd files (image, invvar,
           maskbits, psfsize per band, plus the tractor catalogue) from
           the static file server, sha256-verified
    cut    writes one SCI/IVAR/MASK/LAYERS sample per galaxy under
           $SGADATA/samples/galaxies/, appends a manifest row per
           target, and reports measured volume plus the extrapolation
           to the full selection

The stamp side is --size-mult times D26, floored at --min-size pixels.
Frames are never rejected for content: the only structural criterion is
all-band coverage below the threshold shared with the background driver.

Example:
    python scripts/fetch_galaxies.py --sample 200
"""

import argparse
import csv
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from astropy.io import fits

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import fetch, urls
from src.bricks import Bricks
from src.catalog import Catalog
from src.cutout import PIXSCALE, Coadd, brick_dir, cut, write_sample
from src.mask import DEFAULT_COVERAGE_MIN, compute_layers
from src.select import galaxy_targets, subsample

MANIFEST_FIELDS = [
    "file", "sga_id", "galaxy", "ra", "dec", "brick", "hemisphere", "size_px",
    "psf_g", "psf_r", "psf_z", "valid_frac", "bright_frac", "source_frac",
    "galaxy_frac", "galaxy_bit_frac", "status", "reason",
]


def size_px(target, args):
    return max(int(round(args.size_mult * target["d26"] * 60.0 / PIXSCALE)),
               args.min_size)


def fetch_bricks(bricks, root, workers):
    """Mirror every file of the given (brick, hemisphere) pairs; returns bytes."""
    with ThreadPoolExecutor(max_workers=workers) as pool:
        manifests = list(pool.map(
            lambda b: fetch.fetch_checksums(urls.brick_checksums(*b)), bricks
        ))
    tasks = []
    for (brickname, hemisphere), checksums in zip(bricks, manifests, strict=True):
        directory = brick_dir(root, brickname)
        for name, url in urls.brick_files(brickname, hemisphere).items():
            tasks.append((url, directory / name, checksums.get(name)))
    transferred, failures = fetch.fetch_many(tasks, workers=workers,
                                             desc=f"{len(bricks)} bricks")
    for url, error in failures:
        print(f"FAILED {url}: {error}")
    if failures:
        raise SystemExit(f"{len(failures)} files failed; re-run to resume")
    return transferred


def cut_brick(root, brickname, brick_targets, catalog, out_dir, manifest, args):
    """Cut, mask and write every target of one brick."""
    brick = Coadd(brick_dir(root, brickname), brickname)
    with fits.open(brick_dir(root, brickname) / f"tractor-{brickname}.fits") as hdul:
        tractor = hdul[1].data
    written = 0
    for target in brick_targets:
        path = out_dir / f"sga_{target['sga_id']}.fits"
        if path.exists() and not args.overwrite:
            continue
        side = size_px(target, args)
        sample = cut(brick, target["ra"], target["dec"], side)
        record = {
            "file": path.name, "sga_id": target["sga_id"],
            "galaxy": target["galaxy"],
            "ra": f"{target['ra']:.6f}", "dec": f"{target['dec']:.6f}",
            "brick": brickname, "hemisphere": target["hemisphere"],
            "size_px": side,
            "psf_g": f"{sample['psf_fwhm'][0]:.4f}",
            "psf_r": f"{sample['psf_fwhm'][1]:.4f}",
            "psf_z": f"{sample['psf_fwhm'][2]:.4f}",
            "valid_frac": f"{sample['valid_frac']:.4f}",
        }
        if sample["valid_frac"] < args.coverage_min:
            record.update(status="rejected", reason="coverage")
            manifest.writerow(record)
            continue
        layers, verdict = compute_layers(sample, catalog, tractor)
        record.update({key: f"{value:.4f}" for key, value in verdict.items()})
        record.update(status="written", reason="")
        write_sample(path, sample, brickname,
                     {"SGA_ID": (target["sga_id"], "SGA-2020 identifier")},
                     layers=layers)
        manifest.writerow(record)
        written += 1
    return written


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
    catalog = Catalog()
    all_targets = galaxy_targets(catalog, Bricks(), zcut=args.zcut,
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
        transferred = fetch_bricks(bricks, root, args.workers)
        print(f"fetched {transferred / 1e9:.2f} GB across {len(bricks)} bricks")
    if args.stage == "fetch":
        return

    out_dir = root / "samples" / "galaxies"
    out_dir.mkdir(exist_ok=True, parents=True)
    manifest_path = out_dir / "manifest.csv"
    new_manifest = not manifest_path.exists()
    written = 0
    with open(manifest_path, "a", newline="") as f:
        manifest = csv.DictWriter(f, fieldnames=MANIFEST_FIELDS)
        if new_manifest:
            manifest.writeheader()
        by_brick = {}
        for target in targets:
            by_brick.setdefault(target["brick"], []).append(target)
        for n, (brickname, brick_targets) in enumerate(sorted(by_brick.items()), 1):
            written += cut_brick(root, brickname, brick_targets, catalog,
                                 out_dir, manifest, args)
            print(f"[{n}/{len(by_brick)}] {brickname}: "
                  f"{len(brick_targets)} targets, {written} samples so far")

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
    parser.add_argument("--overwrite", action="store_true",
                        help="rewrite existing samples")
    main(parser.parse_args())
