"""Bright-companion measurement from the tractor catalogue.

A frame has a bright companion when the tractor catalogue holds an
extended source (REX/EXP/DEV/SER — never a PSF star, never an SGA L3
model) inside the stamp, farther than `sep_min_frac` x r26 from the
target centre (closer fits are usually shredded pieces of the target
itself), with r-band flux comparable to the target's own tractor model.
Interacting pairs and satellites below the SGA angular-size threshold —
invisible to every catalogue-based test — are exactly what this catches.

Samples are never touched and nothing is rejected: the measured flux
ratio is written to a sidecar table so training can filter at whatever
threshold it chooses.
"""

import numpy as np


def companion_ratio(tractor, sga_id, ra_deg, dec_deg, half_arcsec,
                    r26_arcsec, sep_min_frac=0.25):
    """Largest companion-to-target flux ratio inside one stamp.

    Returns (ratio, separation in units of r26); (0, 0) when no
    companion qualifies, (nan, nan) when the target has no tractor model
    in this brick to compare against.
    """
    ref_cat = np.asarray(tractor["ref_cat"])
    l3 = ref_cat == "L3"
    target = l3 & (np.asarray(tractor["ref_id"]) == sga_id)
    if not target.any():
        return float("nan"), float("nan")
    target_flux = float(np.asarray(tractor["flux_r"])[target][0])
    if not target_flux > 0:
        return float("nan"), float("nan")

    types = np.char.strip(np.asarray(tractor["type"]).astype(str))
    d_north = (np.asarray(tractor["dec"]) - dec_deg) * 3600.0
    d_east = ((np.asarray(tractor["ra"]) - ra_deg + 180.0) % 360.0 - 180.0) \
        * np.cos(np.radians(dec_deg)) * 3600.0
    sep = np.hypot(d_east, d_north)
    candidates = (
        ~l3
        & ~np.isin(types, ("PSF", "DUP"))
        & (np.abs(d_east) < half_arcsec) & (np.abs(d_north) < half_arcsec)
        & (sep > sep_min_frac * r26_arcsec)
    )
    if not candidates.any():
        return 0.0, 0.0
    ratios = np.asarray(tractor["flux_r"])[candidates] / target_flux
    best = int(np.argmax(ratios))
    return (float(ratios[best]),
            float(sep[candidates][best] / r26_arcsec))
