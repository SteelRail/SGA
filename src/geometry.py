"""Sky geometry: unit vectors, ellipse radii, tangent-plane offsets."""

import numpy as np


def unit_vectors(ra_deg, dec_deg):
    """RA/Dec (deg) to unit vectors on the sphere, shape (N, 3)."""
    ra = np.radians(np.atleast_1d(ra_deg))
    dec = np.radians(np.atleast_1d(dec_deg))
    return np.stack(
        [np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)], axis=-1
    )


def chord_length(angle_deg):
    """Chord distance between unit vectors separated by `angle_deg` degrees."""
    return 2.0 * np.sin(np.radians(angle_deg) / 2.0)


def elliptical_radius(d_east, d_north, pa_deg, ba):
    """Semi-major-equivalent radius of tangent-plane offsets from an ellipse centre.

    `d_east`, `d_north` are offsets in arcsec; `pa_deg` is the position
    angle of the major axis measured from North towards East (the SGA
    convention); `ba` is the minor-to-major axis ratio. Points on the
    ellipse of semi-major axis `a` return exactly `a`.
    """
    pa = np.radians(pa_deg)
    along = d_east * np.sin(pa) + d_north * np.cos(pa)
    across = -d_east * np.cos(pa) + d_north * np.sin(pa)
    return np.hypot(along, across / ba)


def offset_position(ra_deg, dec_deg, sep_arcsec, pa_deg):
    """Position at separation `sep_arcsec` along position angle `pa_deg` (N->E)."""
    d_east = sep_arcsec * np.sin(np.radians(pa_deg)) / 3600.0
    d_north = sep_arcsec * np.cos(np.radians(pa_deg)) / 3600.0
    dec = dec_deg + d_north
    ra = (ra_deg + d_east / np.cos(np.radians(dec))) % 360.0
    return ra, np.clip(dec, -90.0, 90.0)
