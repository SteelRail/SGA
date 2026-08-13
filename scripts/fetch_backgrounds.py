#!/usr/bin/env python
"""Write background training samples from positions stepped off SGA galaxies.

Run scripts/fetch_catalogues.py once first. Backgrounds are drawn from
the neighbourhood of the galaxies they will be paired with, so both
classes share observing conditions, depth and seeing, and both are cut
from the same L1 brick coadds. Per galaxy:

  1. The minimum safe separation comes from the galaxy's own
     surface-brightness profile: the SGA isophotal radii extrapolated to
     --sb-limit mag/arcsec² (where the profile falls below sky noise),
     scaled per direction by the ellipse geometry, plus the stamp
     half-diagonal. No blanket multiple of D26.
  2. --n-candidates positions are generated at randomised position
     angles (seeded per galaxy, so re-runs are reproducible; changing
     --seed, --n-candidates or --sb-limit invalidates existing samples —
     clear the samples directory first) at that separation plus a small
     uniform slack — as close as safety allows, which keeps most
     candidates inside the parent brick.
  3. A candidate is rejected only structurally: outside the footprint,
     within the extent of *any* SGA-2020 galaxy (full catalogue,
     KD-tree, plus DR9's own GALAXY maskbit), or coverage below the
     threshold shared with the galaxy driver. Faint sources never reject
     a candidate — galaxies sit in overdense regions, and pruning
     backgrounds on detections would widen a gap the consumer's
     discriminator would find. Whatever survives is kept; the count per
     galaxy is not fixed.

Bricks are opened through a funnel: candidates in bricks not yet
mirrored are pre-screened on the 0.27 MB maskbits plane alone, and full
imaging is fetched only for bricks that still host a surviving
candidate. Stamp sizes follow the parent galaxy's rule (--size-mult x
D26), so the two classes share a size distribution. Cutting runs across
--jobs worker processes; the manifest is the resume ledger (a row
implies a whole, atomically-written file).

The run reports the survival fraction overall and per brick, the
rejection-reason histogram, and the fraction of accepted backgrounds
sharing the parent's brick. All candidates land in candidates.csv,
accepted samples in manifest.csv next to the samples.

Example:
    python scripts/fetch_backgrounds.py --sample 200
"""

import argparse
import csv
import json
import sys
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import fetch, urls
from src.bricks import Bricks
from src.catalog import Catalog
from src.cutout import PIXSCALE, Coadd, cut, write_sample
from src.geometry import elliptical_radius, offset_position
from src.mask import DEFAULT_COVERAGE_MIN, GALAXY_BIT, compute_layers
from src.select import galaxy_targets, subsample

CANDIDATE_FIELDS = [
    "parent_sga_id", "k", "ra", "dec", "pa_deg", "sep_arcsec", "size_px",
    "brick", "hemisphere", "in_parent_brick", "status", "reason",
]
MANIFEST_FIELDS = [
    "file", "parent_sga_id", "k", "ra", "dec", "brick", "hemisphere",
    "size_px", "sep_arcsec", "pa_deg", "psf_g", "psf_r", "psf_z",
    "valid_frac", "bright_frac", "source_frac", "galaxy_frac",
    "galaxy_bit_frac", "status", "reason",
]

# One catalogue per worker process, loaded on its first brick.
_catalog = None


def propose(catalog, target, args):
    """Seeded candidate positions for one galaxy — pure geometry, no tests.

    Structural tests run vectorised over all candidates afterwards;
    private underscore keys carry the float values they need.
    """
    side = target["size_px"]
    half_diag = 0.5 * np.sqrt(2.0) * side * PIXSCALE

    r_iso = catalog.isophotal_radius(target["index"], args.sb_limit)
    _, ba, pa_major = catalog.ellipse(target["index"])

    rng = np.random.default_rng([args.seed, target["sga_id"]])
    pas = rng.uniform(0.0, 360.0, args.n_candidates)
    slack = rng.uniform(0.0, args.radius_slack, args.n_candidates)

    # A unit step along `pa` has elliptical radius f(pa), so the
    # iso-contour at semi-major r_iso lies at r_iso / f(pa) in that
    # direction — closer along the minor axis.
    factors = elliptical_radius(np.sin(np.radians(pas)),
                                np.cos(np.radians(pas)), pa_major, ba)
    seps = (r_iso / factors + half_diag) * (1.0 + slack)
    ras, decs = offset_position(target["ra"], target["dec"], seps, pas)

    return [
        {
            "parent_sga_id": target["sga_id"], "k": k,
            "ra": f"{ras[k]:.6f}", "dec": f"{decs[k]:.6f}",
            "pa_deg": f"{pas[k]:.2f}", "sep_arcsec": f"{seps[k]:.1f}",
            "size_px": side, "status": "candidate", "reason": "",
            "_ra": float(ras[k]), "_dec": float(decs[k]),
            "_half_diag": half_diag, "_index": target["index"],
            "_parent_brick": target["brick"],
        }
        for k in range(args.n_candidates)
    ]


def structural_tests(catalog, bricks, candidates, args):
    """Footprint and full-atlas clearance for every candidate, vectorised."""
    ras = np.array([c["_ra"] for c in candidates])
    decs = np.array([c["_dec"] for c in candidates])
    names, hemis, in_footprint = bricks.resolve(ras, decs,
                                                min_nexp=args.min_nexp)
    clear = catalog.clearance_ok(
        ras, decs,
        extra_arcsec=np.array([c["_half_diag"] for c in candidates]),
        margin=args.sga_margin,
        exclude=np.array([c["_index"] for c in candidates]),
    )
    for i, candidate in enumerate(candidates):
        candidate["brick"] = str(names[i])
        candidate["hemisphere"] = str(hemis[i])
        candidate["in_parent_brick"] = \
            int(candidate["brick"] == candidate["_parent_brick"])
        if not in_footprint[i]:
            candidate.update(status="rejected", reason="footprint")
        elif not clear[i]:
            candidate.update(status="rejected", reason="sga_overlap")


def maskbits_prescreen(root, candidates, workers, coverage_min):
    """Screen candidates on the maskbits plane alone, fetching it where needed.

    Rejects candidates whose stamp box cannot reach the coverage
    threshold on the brick pixel grid, or contains DR9 GALAXY-bit pixels
    (an SGA galaxy footprint the catalogue test missed). Costs 0.27 MB
    per new brick instead of ~90 MB.
    """
    alive = [c for c in candidates if c["status"] == "candidate"]
    bricks = sorted({(c["brick"], c["hemisphere"]) for c in alive})
    tasks, owner = [], {}
    for brickname, hemisphere in bricks:
        name = f"legacysurvey-{brickname}-maskbits.fits.fz"
        url = urls.brick_files(brickname, hemisphere)[name]
        tasks.append((url, fetch.brick_dir(root, brickname) / name))
        owner[url] = brickname
    transferred, failures = fetch.fetch_many(tasks, workers=workers,
                                             desc="maskbits funnel")
    failed = {owner[url] for url, _ in failures}

    by_brick = defaultdict(list)
    for candidate in alive:
        by_brick[candidate["brick"]].append(candidate)
    for brickname, brick_candidates in by_brick.items():
        if brickname in failed:
            for candidate in brick_candidates:
                candidate.update(status="rejected", reason="fetch")
            continue
        path = fetch.brick_dir(root, brickname) \
            / f"legacysurvey-{brickname}-maskbits.fits.fz"
        with fits.open(path) as hdul:
            hdu = hdul[1] if hdul[0].data is None else hdul[0]
            maskbits, wcs = hdu.data, WCS(hdu.header)
        for candidate in brick_candidates:
            x, y = wcs.world_to_pixel_values(candidate["_ra"],
                                             candidate["_dec"])
            half = candidate["size_px"] / 2.0
            x_lo, x_hi = int(x - half), int(x + half)
            y_lo, y_hi = int(y - half), int(y + half)
            box = maskbits[max(y_lo, 0):max(y_hi, 0),
                           max(x_lo, 0):max(x_hi, 0)]
            on_grid = box.size / float(candidate["size_px"]) ** 2
            if on_grid < coverage_min:
                candidate.update(status="rejected", reason="coverage")
            elif np.any(box & GALAXY_BIT):
                candidate.update(status="rejected", reason="sga_maskbit")
    return transferred


def cut_brick(task):
    """Cut, mask and write one brick's candidates; returns manifest rows."""
    global _catalog
    if _catalog is None:
        _catalog = Catalog()
    root, out_dir, brickname, brick_candidates, coverage_min = task
    brick = Coadd(fetch.brick_dir(root, brickname))
    with fits.open(fetch.brick_dir(root, brickname)
                   / f"tractor-{brickname}.fits") as hdul:
        tractor = hdul[1].data
    records = []
    for candidate in brick_candidates:
        path = Path(out_dir) \
            / f"bg_{candidate['parent_sga_id']}_{candidate['k']}.fits"
        record = {key: candidate[key] for key in MANIFEST_FIELDS
                  if key in candidate}
        record["file"] = path.name
        sample = cut(brick, candidate["_ra"], candidate["_dec"],
                     candidate["size_px"])
        for band, fwhm in zip("grz", sample["psf_fwhm"], strict=True):
            record[f"psf_{band}"] = f"{fwhm:.4f}"
        record["valid_frac"] = f"{sample['valid_frac']:.4f}"
        if sample["valid_frac"] < coverage_min:
            record.update(status="rejected", reason="coverage")
            records.append(record)
            continue
        layers, verdict = compute_layers(sample, _catalog, tractor)
        record.update({key: f"{value:.4f}" for key, value in verdict.items()})
        record.update(status="written", reason="")
        provenance = {
            "PARENT": (candidate["parent_sga_id"],
                       "SGA_ID of the parent galaxy"),
            "OFFSET": (float(candidate["sep_arcsec"]),
                       "separation from parent (arcsec)"),
            "OFFPA": (float(candidate["pa_deg"]),
                      "position angle of offset (deg, N->E)"),
        }
        write_sample(path, sample, brickname, provenance, layers)
        records.append(record)
    return records


def report(candidates, targets):
    """Survival fraction — the number this driver exists to measure."""
    total = len(candidates)
    accepted = [c for c in candidates if c["status"] == "accepted"]
    reasons = Counter(c["reason"] for c in candidates
                      if c["status"] == "rejected")
    print(f"\ncandidates: {total} over {len(targets)} galaxies")
    print(f"accepted: {len(accepted)} ({len(accepted) / total:.1%} survival)")
    print("rejections:", dict(reasons) or "none")
    if accepted:
        in_parent = np.mean([c["in_parent_brick"] for c in accepted])
        seps = np.array([float(c["sep_arcsec"]) for c in accepted])
        print(f"accepted in parent brick: {in_parent:.1%}")
        print(f"separations (arcsec): median {np.median(seps):.0f}, "
              f"p90 {np.percentile(seps, 90):.0f}")
        per_brick = defaultdict(lambda: [0, 0])
        for c in candidates:
            per_brick[c["brick"]][1] += 1
            if c["status"] == "accepted":
                per_brick[c["brick"]][0] += 1
        fractions = np.array([a / t for a, t in per_brick.values()])
        print(f"per-brick survival: median {np.median(fractions):.1%}, "
              f"quartiles {np.percentile(fractions, 25):.1%}-"
              f"{np.percentile(fractions, 75):.1%}")
        per_galaxy = Counter(c["parent_sga_id"] for c in accepted)
        counts = [per_galaxy.get(t["sga_id"], 0) for t in targets]
        print(f"accepted per galaxy: mean {np.mean(counts):.1f}, "
              f"zero for {np.mean(np.array(counts) == 0):.1%} of galaxies")


def main(args):
    root = Path(args.root)
    catalog = Catalog()
    bricks = Bricks()
    targets = subsample(
        galaxy_targets(catalog, bricks, zcut=args.zcut, min_nexp=args.min_nexp,
                       size_mult=args.size_mult, min_size=args.min_size,
                       sga_margin=args.sga_margin),
        args.sample, args.limit)
    print(f"{len(targets)} parent galaxies")

    candidates = []
    for target in targets:
        candidates.extend(propose(catalog, target, args))
    structural_tests(catalog, bricks, candidates, args)
    alive = sum(c["status"] == "candidate" for c in candidates)
    print(f"{len(candidates)} candidates, {alive} past catalogue tests")

    transferred = maskbits_prescreen(root, candidates, args.workers,
                                     args.coverage_min)
    alive = sum(c["status"] == "candidate" for c in candidates)
    print(f"{alive} past the maskbits funnel "
          f"({transferred / 1e6:.0f} MB of maskbits)")

    needed = sorted({(c["brick"], c["hemisphere"]) for c in candidates
                     if c["status"] == "candidate"})
    fetched, failed = fetch.mirror_bricks(needed, root, args.workers,
                                          desc=f"{len(needed)} bricks")
    transferred += fetched
    for candidate in candidates:
        if candidate["status"] == "candidate" and candidate["brick"] in failed:
            candidate.update(status="rejected", reason="fetch")
    print(f"{transferred / 1e9:.2f} GB fetched")

    out_dir = root / "samples" / "backgrounds"
    out_dir.mkdir(exist_ok=True, parents=True)
    manifest_path = out_dir / "manifest.csv"
    # the manifest is only a valid ledger under the parameters that built it
    spec = {key: getattr(args, key) for key in
            ("zcut", "min_nexp", "size_mult", "min_size", "sga_margin",
             "coverage_min", "sb_limit", "n_candidates", "radius_slack",
             "seed")}
    spec_path = out_dir / "run.json"
    if args.overwrite:
        manifest_path.unlink(missing_ok=True)
        spec_path.unlink(missing_ok=True)
    if spec_path.exists() and json.loads(spec_path.read_text()) != spec:
        raise SystemExit(f"{spec_path} records different science parameters; "
                         "wipe the samples directory or pass --overwrite")
    spec_path.write_text(json.dumps(spec, indent=1))
    done = {}
    if manifest_path.exists():
        with open(manifest_path) as f:
            done = {row["file"]: row for row in csv.DictReader(f)}

    by_brick = defaultdict(list)
    for candidate in candidates:
        if candidate["status"] != "candidate":
            continue
        name = f"bg_{candidate['parent_sga_id']}_{candidate['k']}.fits"
        if name in done:
            row = done[name]
            candidate["status"] = ("accepted" if row["status"] == "written"
                                   else "rejected")
            candidate["reason"] = row["reason"]
        else:
            by_brick[candidate["brick"]].append(candidate)

    tasks = [(str(root), str(out_dir), brickname, brick_candidates,
              args.coverage_min)
             for brickname, brick_candidates in sorted(by_brick.items())]
    outcome = {}
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
                    outcome[record["file"]] = record
                    written += record["status"] == "written"
    if failed:
        print(f"{failed} bricks failed; re-run to retry them")
    for brick_candidates in by_brick.values():
        for candidate in brick_candidates:
            record = outcome.get(
                f"bg_{candidate['parent_sga_id']}_{candidate['k']}.fits")
            if record is None:
                candidate.update(status="rejected", reason="error")
                continue
            candidate["status"] = ("accepted" if record["status"] == "written"
                                   else "rejected")
            candidate["reason"] = record["reason"]
    print(f"{written} samples written to {out_dir}")

    with open(out_dir / "candidates.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=CANDIDATE_FIELDS)
        writer.writeheader()
        writer.writerows(
            {key: c[key] for key in CANDIDATE_FIELDS} for c in candidates
        )
    report(candidates, targets)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__.split("\n")[0],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--root", default="data/sga",
                        help="data root: brick mirror and samples")
    parser.add_argument("--zcut", type=float, default=0.05,
                        help="upper redshift cut on Z_LEDA")
    parser.add_argument("--sb-limit", type=float, default=28.0,
                        help="r-band isophote (mag/arcsec²) the separation clears")
    parser.add_argument("--n-candidates", type=int, default=8,
                        help="candidate positions per galaxy")
    parser.add_argument("--radius-slack", type=float, default=0.5,
                        help="uniform outward slack on the minimum separation")
    parser.add_argument("--sga-margin", type=float, default=1.5,
                        help="clearance margin on other galaxies' r26")
    parser.add_argument("--size-mult", type=float, default=1.5,
                        help="stamp side as a multiple of the parent's D26")
    parser.add_argument("--min-size", type=int, default=64,
                        help="minimum stamp side (pixels)")
    parser.add_argument("--min-nexp", type=int, default=1,
                        help="minimum exposures per band for footprint membership")
    parser.add_argument("--coverage-min", type=float, default=DEFAULT_COVERAGE_MIN,
                        help="minimum all-band coverage fraction (structural)")
    parser.add_argument("--seed", type=int, default=20260807,
                        help="base seed; per-galaxy streams derive from it")
    parser.add_argument("--sample", type=int, default=None,
                        help="evenly RA-spaced subset of this many parents")
    parser.add_argument("--limit", type=int, default=None,
                        help="first N parents only")
    parser.add_argument("--workers", type=int, default=12,
                        help="parallel download connections")
    parser.add_argument("--jobs", type=int, default=24,
                        help="cut-stage worker processes; cold NFS reads bound "
                             "the stage, and 24 saturates the mount")
    parser.add_argument("--overwrite", action="store_true",
                        help="reset the sample ledger and re-cut this run's candidates")
    main(parser.parse_args())
