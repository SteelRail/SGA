"""Galaxy target selection: redshift cut, footprint, brick resolution, isolation.

The selection is deterministic and recomputed per run from the two
catalogue sources — nothing is cached to disk. The isolation test keeps
single-galaxy frames only: a target whose stamp (`size_mult` x D26,
floored at `min_size` pixels) is reached by any *other* catalogued
galaxy's `sga_margin` x r26 ellipse is dropped — the same clearance test
the background driver applies to its candidates.
"""

import numpy as np

from .cutout import PIXSCALE


def galaxy_targets(catalog, bricks, zcut=0.05, min_nexp=1,
                   size_mult=1.5, min_size=64, sga_margin=1.5):
    """Every selected, isolated galaxy with its resolved brick, as dicts.

    Sorted by RA, so evenly spaced subsets spread over the RA range.
    """
    selected = catalog.select(zcut=zcut)
    rows = catalog.rows[selected]
    names, hemis, ok = bricks.resolve(rows["RA"], rows["DEC"], min_nexp=min_nexp)
    side = np.maximum(
        np.round(size_mult * rows["D26"].astype(float) * 60.0 / PIXSCALE)
        .astype(int), min_size)
    ok &= catalog.clearance_ok(
        rows["RA"], rows["DEC"],
        extra_arcsec=0.5 * np.sqrt(2.0) * side * PIXSCALE,
        margin=sga_margin, exclude=selected,
    )
    targets = [
        {
            "index": int(selected[i]),
            "sga_id": int(rows["SGA_ID"][i]),
            "galaxy": str(rows["GALAXY"][i]),
            "ra": float(rows["RA"][i]),
            "dec": float(rows["DEC"][i]),
            "d26": float(rows["D26"][i]),
            "size_px": int(side[i]),
            "group_name": str(rows["GROUP_NAME"][i]),
            "group_ra": float(rows["GROUP_RA"][i]),
            "brick": str(names[i]),
            "hemisphere": str(hemis[i]),
        }
        for i in np.flatnonzero(ok)
    ]
    targets.sort(key=lambda t: t["ra"])
    return targets


def subsample(targets, sample=None, limit=None):
    """Evenly RA-spaced subset of `sample` targets, then the first `limit`."""
    if sample and sample < len(targets):
        picks = np.unique(
            np.linspace(0, len(targets) - 1, sample).round().astype(int)
        )
        targets = [targets[i] for i in picks]
    if limit:
        targets = targets[:limit]
    return targets
