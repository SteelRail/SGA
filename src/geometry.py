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
    """Position at separation `sep_arcsec` along position angle `pa_deg` (N->E).

    Exact great-circle offset, so the realised separation equals the
    requested one at any declination.
    """
    dec0 = np.radians(dec_deg)
    sep = np.radians(sep_arcsec / 3600.0)
    pa = np.radians(pa_deg)
    sin_dec = np.clip(
        np.sin(dec0) * np.cos(sep) + np.cos(dec0) * np.sin(sep) * np.cos(pa),
        -1.0, 1.0)
    dra = np.arctan2(np.sin(pa) * np.sin(sep) * np.cos(dec0),
                     np.cos(sep) - np.sin(dec0) * sin_dec)
    return ((np.asarray(ra_deg) + np.degrees(dra)) % 360.0,
            np.degrees(np.arcsin(sin_dec)))
