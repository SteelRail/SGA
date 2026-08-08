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
     angles (seeded per galaxy, so re-runs are reproducible) at that
     separation plus a small uniform slack — as close as safety allows,
     which keeps most candidates inside the parent brick.
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
D26), so the two classes share a size distribution.

The run reports the survival fraction overall and per brick, the
rejection-reason histogram, and the fraction of accepted backgrounds
sharing the parent's brick. All candidates land in candidates.csv,
accepted samples in manifest.csv next to the samples.

Example:
    python scripts/fetch_backgrounds.py --sample 200
"""

import argparse
import csv
import sys
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import fetch, urls
from src.bricks import Bricks
from src.catalog import Catalog
from src.cutout import PIXSCALE, Coadd, brick_dir, cut, write_sample
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


def propose(catalog, bricks, target, args):
    """Candidate positions for one galaxy, tested against catalogue and footprint."""
    side = max(int(round(args.size_mult * target["d26"] * 60.0 / PIXSCALE)),
               args.min_size)
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

    names, hemis, in_footprint = bricks.resolve(ras, decs,
                                                min_nexp=args.min_nexp)
    clear = catalog.clearance_ok(ras, decs, extra_arcsec=half_diag,
                                 margin=args.sga_margin,
                                 exclude=target["index"])
    candidates = []
    for k in range(args.n_candidates):
        candidate = {
            "parent_sga_id": target["sga_id"], "k": k,
            "ra": f"{ras[k]:.6f}", "dec": f"{decs[k]:.6f}",
            "pa_deg": f"{pas[k]:.2f}", "sep_arcsec": f"{seps[k]:.1f}",
            "size_px": side, "brick": str(names[k]), "hemisphere": str(hemis[k]),
            "in_parent_brick": int(str(names[k]) == target["brick"]),
            "status": "candidate", "reason": "",
        }
        if not in_footprint[k]:
            candidate.update(status="rejected", reason="footprint")
        elif not clear[k]:
            candidate.update(status="rejected", reason="sga_overlap")
        candidates.append(candidate)
    return candidates


def maskbits_prescreen(root, candidates, workers):
    """Screen candidates on the maskbits plane alone, fetching it where needed.

    Rejects candidates whose stamp box cannot reach the coverage
    threshold on the brick pixel grid, or contains DR9 GALAXY-bit pixels
    (an SGA galaxy footprint the catalogue test missed). Costs 0.27 MB
    per new brick instead of ~90 MB.
    """
    alive = [c for c in candidates if c["status"] == "candidate"]
    bricks = sorted({(c["brick"], c["hemisphere"]) for c in alive})
    tasks = []
    for brickname, hemisphere in bricks:
        name = f"legacysurvey-{brickname}-maskbits.fits.fz"
        url = urls.brick_files(brickname, hemisphere)[name]
        tasks.append((url, brick_dir(root, brickname) / name))
    transferred, failures = fetch.fetch_many(tasks, workers=workers,
                                             desc="maskbits funnel")
    failed = {url.rsplit("/", 2)[-2] for url, _ in failures}

    by_brick = defaultdict(list)
    for candidate in alive:
        by_brick[candidate["brick"]].append(candidate)
    for brickname, brick_candidates in by_brick.items():
        if brickname in failed:
            for candidate in brick_candidates:
                candidate.update(status="rejected", reason="fetch")
            continue
        path = brick_dir(root, brickname) \
            / f"legacysurvey-{brickname}-maskbits.fits.fz"
        with fits.open(path) as hdul:
            hdu = hdul[1] if hdul[0].data is None else hdul[0]
            maskbits, wcs = hdu.data, WCS(hdu.header)
        for candidate in brick_candidates:
            x, y = wcs.world_to_pixel_values(float(candidate["ra"]),
                                             float(candidate["dec"]))
            half = int(candidate["size_px"]) / 2.0
            x_lo, x_hi = int(x - half), int(x + half)
            y_lo, y_hi = int(y - half), int(y + half)
            box = maskbits[max(y_lo, 0):max(y_hi, 0),
                           max(x_lo, 0):max(x_hi, 0)]
            on_grid = box.size / float(candidate["size_px"]) ** 2
            if on_grid < DEFAULT_COVERAGE_MIN:
                candidate.update(status="rejected", reason="coverage")
            elif np.any(box & GALAXY_BIT):
                candidate.update(status="rejected", reason="sga_maskbit")
    return transferred


def fetch_imaging(root, candidates, workers):
    """Full coadd files for every brick still hosting a live candidate."""
    alive = [c for c in candidates if c["status"] == "candidate"]
    bricks = sorted({(c["brick"], c["hemisphere"]) for c in alive})
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
    failed = {url.rsplit("/", 2)[-2] for url, _ in failures}
    for candidate in alive:
        if candidate["brick"] in failed:
            candidate.update(status="rejected", reason="fetch")
    return transferred


def cut_candidates(root, candidates, catalog, out_dir, manifest, args):
    """Cut, mask and write every surviving candidate, brick by brick."""
    by_brick = defaultdict(list)
    for candidate in candidates:
        if candidate["status"] == "candidate":
            by_brick[candidate["brick"]].append(candidate)
    written = 0
    for n, (brickname, brick_candidates) in enumerate(sorted(by_brick.items()), 1):
        brick = Coadd(brick_dir(root, brickname), brickname)
        with fits.open(brick_dir(root, brickname)
                       / f"tractor-{brickname}.fits") as hdul:
            tractor = hdul[1].data
        for candidate in brick_candidates:
            path = out_dir \
                / f"bg_{candidate['parent_sga_id']}_{candidate['k']}.fits"
            record = {key: candidate[key] for key in MANIFEST_FIELDS
                      if key in candidate}
            record["file"] = path.name
            if path.exists() and not args.overwrite:
                candidate["status"] = "accepted"
                continue
            sample = cut(brick, float(candidate["ra"]),
                         float(candidate["dec"]), int(candidate["size_px"]))
            for band, fwhm in zip("grz", sample["psf_fwhm"], strict=True):
                record[f"psf_{band}"] = f"{fwhm:.4f}"
            record["valid_frac"] = f"{sample['valid_frac']:.4f}"
            if sample["valid_frac"] < args.coverage_min:
                candidate.update(status="rejected", reason="coverage")
                record.update(status="rejected", reason="coverage")
                manifest.writerow(record)
                continue
            layers, verdict = compute_layers(sample, catalog, tractor)
            record.update({key: f"{value:.4f}"
                           for key, value in verdict.items()})
            record.update(status="written", reason="")
            provenance = {
                "PARENT": (candidate["parent_sga_id"],
                           "SGA_ID of the parent galaxy"),
                "OFFSET": (float(candidate["sep_arcsec"]),
                           "separation from parent (arcsec)"),
                "OFFPA": (float(candidate["pa_deg"]),
                          "position angle of offset (deg, N->E)"),
            }
            write_sample(path, sample, brickname, provenance, layers=layers)
            manifest.writerow(record)
            candidate["status"] = "accepted"
            written += 1
        print(f"[{n}/{len(by_brick)}] {brickname}: "
              f"{len(brick_candidates)} candidates, {written} written so far")
    return written


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
        candidates.extend(propose(catalog, bricks, target, args))
    alive = sum(c["status"] == "candidate" for c in candidates)
    print(f"{len(candidates)} candidates, {alive} past catalogue tests")

    transferred = maskbits_prescreen(root, candidates, args.workers)
    alive = sum(c["status"] == "candidate" for c in candidates)
    print(f"{alive} past the maskbits funnel "
          f"({transferred / 1e6:.0f} MB of maskbits)")

    transferred += fetch_imaging(root, candidates, args.workers)
    print(f"{transferred / 1e9:.2f} GB fetched")

    out_dir = root / "samples" / "backgrounds"
    out_dir.mkdir(exist_ok=True, parents=True)
    manifest_path = out_dir / "manifest.csv"
    new_manifest = not manifest_path.exists()
    with open(manifest_path, "a", newline="") as f:
        manifest = csv.DictWriter(f, fieldnames=MANIFEST_FIELDS)
        if new_manifest:
            manifest.writeheader()
        written = cut_candidates(root, candidates, catalog, out_dir,
                                 manifest, args)
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
    parser.add_argument("--overwrite", action="store_true",
                        help="rewrite existing samples")
    main(parser.parse_args())
