"""The SGA-2020 catalogue: selection, ellipse geometry, full-atlas clearance.

Isophotal extrapolation
-----------------------
The SGA ellipse fit provides semi-major axes at half-magnitude steps of
r-band surface brightness down to 26 mag/arcsec² (SMA_SB25, SMA_SB26).
Outside that range the profile is extrapolated linearly in radius per
magnitude, the exact behaviour of an exponential disk (mu grows by
1.086 r/h):

    r(mu) = r26 + (mu - 26) * dr,   dr = r26 - r25

with dr floored at 0.1 * r26 to guard against degenerate fits where the
two isophotes nearly coincide. For the ~5% of galaxies without an ellipse
fit (SMA_SB26 <= 0, D26 inherited from LEDA), dr falls back to 0.25 * r26,
the median outer slope of the fitted population.
"""

import numpy as np
from astropy.io import fits
from scipy.spatial import cKDTree

from .geometry import chord_length, elliptical_radius, unit_vectors


class Catalog:
    """The full 383,620-row SGA-2020 table, loaded once."""

    def __init__(self, path="data/catalogues/SGA-2020.fits"):
        with fits.open(path) as hdul:
            self.rows = hdul[1].data
        # Plain-array copies of the geometry columns: fancy-indexing a
        # FITS record array inside the clearance loop is orders of
        # magnitude slower.
        r26 = np.array(self.rows["SMA_SB26"], dtype=float)
        self.r26 = np.where(r26 > 0, r26, self.rows["D26"] * 30.0)  # arcsec
        self._ra = np.array(self.rows["RA"], dtype=float)
        self._dec = np.array(self.rows["DEC"], dtype=float)
        ba = np.array(self.rows["BA"], dtype=float)
        self._ba = np.where((ba > 0) & (ba <= 1), ba, 1.0)
        self._pa = np.nan_to_num(np.array(self.rows["PA"], dtype=float))
        self._tree = None

    def select(self, zcut=0.05):
        """Indices of galaxies with 0 < Z_LEDA < `zcut`."""
        z = self.rows["Z_LEDA"]
        return np.flatnonzero((z > 0) & (z < zcut))

    def ellipse(self, index):
        """(r26_arcsec, ba, pa_deg) for one row; ba and pa fall back to round."""
        return (float(self.r26[index]), float(self._ba[index]),
                float(self._pa[index]))

    def isophotal_radius(self, index, sb_limit):
        """Extrapolated semi-major axis (arcsec) at `sb_limit` mag/arcsec² in r.

        See the module docstring for the derivation and the fallbacks.
        """
        row = self.rows[index]
        r25, r26 = float(row["SMA_SB25"]), float(row["SMA_SB26"])
        if r26 > 0:
            dr = r26 - r25 if r25 > 0 else 0.25 * r26
            dr = max(dr, 0.1 * r26)
        else:
            r26 = float(row["D26"]) * 30.0
            dr = 0.25 * r26
        return r26 + max(sb_limit - 26.0, 0.0) * dr

    # -- full-atlas clearance (KD-tree over unit vectors) ------------------

    def _neighbours(self, points, reach_arcsec):
        """Catalogue indices within `reach_arcsec` of each query point."""
        if self._tree is None:
            self._tree = cKDTree(unit_vectors(self._ra, self._dec))
        return self._tree.query_ball_point(
            points, chord_length(reach_arcsec / 3600.0)
        )

    def clearance_ok(self, ra_deg, dec_deg, extra_arcsec=0.0, margin=1.5,
                     exclude=None):
        """True where a position is clear of every catalogued galaxy's extent.

        A position fails when it lies within `margin` times any galaxy's
        r26 along the ellipse direction, plus `extra_arcsec` (use the
        stamp half-diagonal so the whole stamp stays clear, not just its
        centre). `exclude` drops one catalogue index per position — the
        query galaxy itself, whose own extent must not fail its own test.
        """
        ra_deg = np.atleast_1d(np.asarray(ra_deg, dtype=float))
        dec_deg = np.atleast_1d(np.asarray(dec_deg, dtype=float))
        extra_arcsec = np.broadcast_to(
            np.asarray(extra_arcsec, dtype=float), ra_deg.shape)
        if exclude is not None:
            exclude = np.broadcast_to(np.asarray(exclude, dtype=int),
                                      ra_deg.shape)
        reach = margin * self.r26.max() + extra_arcsec.max()
        neighbour_lists = self._neighbours(unit_vectors(ra_deg, dec_deg), reach)
        ok = np.ones(len(ra_deg), dtype=bool)
        for i, neighbours in enumerate(neighbour_lists):
            indices = np.asarray(neighbours, dtype=int)
            if exclude is not None:
                indices = indices[indices != exclude[i]]
            if not len(indices):
                continue
            d_north = (dec_deg[i] - self._dec[indices]) * 3600.0
            d_east = ((ra_deg[i] - self._ra[indices] + 180.0) % 360.0 - 180.0) \
                * np.cos(np.radians(dec_deg[i])) * 3600.0
            r_ell = elliptical_radius(d_east, d_north,
                                      self._pa[indices], self._ba[indices])
            ok[i] = bool(np.all(
                r_ell > margin * self.r26[indices] + extra_arcsec[i]
            ))
        return ok

    def overlapping(self, ra_deg, dec_deg, half_diag_arcsec, margin=1.0):
        """Indices whose `margin`*r26 ellipse may reach a stamp at the position.

        Conservative (circular) pre-selection for painting the galaxy mask
        layer; the painter applies the exact ellipse geometry per galaxy.
        """
        reach = margin * self.r26.max() + half_diag_arcsec
        point = unit_vectors(ra_deg, dec_deg)
        indices = np.asarray(self._neighbours(point, reach)[0], dtype=int)
        if not len(indices):
            return indices
        cos_sep = np.clip(
            unit_vectors(self._ra[indices], self._dec[indices]) @ point[0],
            -1, 1)
        sep = np.degrees(np.arccos(cos_sep)) * 3600.0
        return indices[sep <= margin * self.r26[indices] + half_diag_arcsec]
