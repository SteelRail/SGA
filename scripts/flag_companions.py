#!/usr/bin/env python
"""Measure bright-companion flux ratios for every written galaxy sample.

Sweeps the tractor catalogue of each brick (which must still be
mirrored) and writes companions.csv next to the manifest: one row per
written sample with the largest companion-to-target flux ratio and its
separation. Samples are never modified — training filters on the ratio
at whatever threshold it chooses (see src/companions.py for the
criterion).

Example:
    python scripts/flag_companions.py
"""

import argparse
import csv
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
from astropy.io import fits
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import fetch
from src.catalog import Catalog
from src.companions import companion_ratio
from src.cutout import PIXSCALE

_catalog = None


def flag_brick(task):
    """Companion ratios for one brick's samples; returns csv rows."""
    global _catalog
    if _catalog is None:
        _catalog = Catalog()
    root, brickname, brick_rows = task
    with fits.open(fetch.brick_dir(root, brickname)
                   / f"tractor-{brickname}.fits") as hdul:
        tractor = hdul[1].data
    out = []
    for r in brick_rows:
        sga_id = int(r["sga_id"])
        index = int(np.flatnonzero(_catalog.rows["SGA_ID"] == sga_id)[0])
        ratio, sep = companion_ratio(
            tractor, sga_id, float(r["ra"]), float(r["dec"]),
            0.5 * int(r["size_px"]) * PIXSCALE, _catalog.r26[index],
        )
        out.append({"file": r["file"], "companion_ratio": f"{ratio:.4f}",
                    "companion_sep_r26": f"{sep:.2f}"})
    return out


def main(args):
    root = Path(args.root)
    out_dir = root / "samples" / "galaxies"
    with open(out_dir / "manifest.csv") as f:
        rows = [r for r in csv.DictReader(f) if r["status"] == "written"]
    by_brick = {}
    for r in rows:
        by_brick.setdefault(r["brick"], []).append(r)
    tasks = [(str(root), brickname, brick_rows)
             for brickname, brick_rows in sorted(by_brick.items())]

    results = []
    with ProcessPoolExecutor(max_workers=args.jobs) as pool:
        futures = [pool.submit(flag_brick, task) for task in tasks]
        for future in tqdm(as_completed(futures), total=len(futures),
                           desc="flagging", unit="brick"):
            results.extend(future.result())

    with open(out_dir / "companions.csv", "w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=["file", "companion_ratio", "companion_sep_r26"])
        writer.writeheader()
        writer.writerows(sorted(results, key=lambda r: r["file"]))
    ratios = np.array([float(r["companion_ratio"]) for r in results])
    ok = np.isfinite(ratios)
    print(f"{len(results)} samples flagged -> {out_dir / 'companions.csv'}")
    for threshold in (0.1, 0.2, 0.5):
        print(f"  ratio > {threshold}: {(ratios[ok] > threshold).sum()} "
              f"({(ratios[ok] > threshold).mean():.2%})")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__.split("\n")[0],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--root", default="data/sga",
                        help="data root: brick mirror and samples")
    parser.add_argument("--jobs", type=int, default=16,
                        help="worker processes (tractor reads are the cost)")
    main(parser.parse_args())
