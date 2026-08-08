#!/usr/bin/env python
"""Fetch the reference catalogues: SGA-2020 and the DR9 survey-bricks tables.

A one-off prerequisite for every other script (~0.8 GB into
data/catalogues/). Existing files are never re-fetched; an interrupted
download resumes.

Example:
    python scripts/fetch_catalogues.py
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import fetch
from src.urls import CATALOGUES


def main(args):
    dest = Path(args.dest)
    tasks = [(url, dest / name) for name, url in CATALOGUES.items()]
    transferred, failures = fetch.fetch_many(tasks, workers=len(tasks),
                                             desc="catalogues")
    for url, error in failures:
        print(f"FAILED {url}: {error}")
    if failures:
        raise SystemExit(f"{len(failures)} files failed; re-run to resume")
    print(f"{transferred / 1e9:.2f} GB fetched to {dest}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__.split("\n")[0],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--dest", default="data/catalogues",
                        help="destination directory")
    main(parser.parse_args())
