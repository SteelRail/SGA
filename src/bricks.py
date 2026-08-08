"""The DR9 survey-bricks tables: which brick covers a sky position."""

import numpy as np
from astropy.io import fits
from scipy.spatial import cKDTree

from .geometry import unit_vectors

# DR9 combines the two hemisphere reductions at this declination: south
# (DECam) below, north (BASS+MzLS) above. The rule is applied to the
# brick centre, so a brick name always maps to one hemisphere.
NORTH_SOUTH_DEC = 32.375


class Bricks:
    """Both hemisphere tables, loaded once; one operation: position -> brick."""

    def __init__(self,
                 north="data/catalogues/survey-bricks-dr9-north.fits.gz",
                 south="data/catalogues/survey-bricks-dr9-south.fits.gz"):
        self.tables = {}
        self._trees = {}
        for hemisphere, path in (("north", north), ("south", south)):
            with fits.open(path) as hdul:
                table = hdul[1].data
            self.tables[hemisphere] = table
            self._trees[hemisphere] = cKDTree(
                unit_vectors(table["ra"], table["dec"])
            )

    def _covered(self, points, ra_deg, dec_deg, hemisphere, min_nexp):
        """Covering-brick index and coverage flag for each position.

        Bricks tile the sky in 0.25-deg declination rows whose RA widths
        differ row to row, so the nearest brick *centre* is often in the
        adjacent row for positions near a row boundary. The containing
        brick is always among the nearest few centres, so the bounds test
        (ra1/ra2/dec1/dec2 with RA wrap-around, plus `min_nexp` exposures
        in each of g, r and z) runs over the 8 nearest and keeps the one
        that passes — brick bounds are disjoint, so at most one does.
        """
        table = self.tables[hemisphere]
        _, near = self._trees[hemisphere].query(points, k=8)
        ok = np.ones(near.shape, dtype=bool)
        for band in ("g", "r", "z"):
            ok &= table[f"nexp_{band}"][near] >= min_nexp
        ra1, ra2 = table["ra1"][near], table["ra2"][near]
        width = (ra2 - ra1 + 360.0) % 360.0
        dra = (ra_deg[:, None] - ra1 + 360.0) % 360.0
        ok &= dra <= width
        ok &= (dec_deg[:, None] >= table["dec1"][near]) \
            & (dec_deg[:, None] <= table["dec2"][near])
        first = ok.argmax(axis=1)
        idx = near[np.arange(len(near)), first]
        return idx, ok.any(axis=1)

    def resolve(self, ra_deg, dec_deg, min_nexp=1):
        """Brick name and hemisphere for each position.

        Returns (bricknames, hemispheres, ok): `ok` is False where no
        covering brick exists — the footprint test.
        """
        ra_deg = np.atleast_1d(np.asarray(ra_deg, dtype=float))
        dec_deg = np.atleast_1d(np.asarray(dec_deg, dtype=float))
        points = unit_vectors(ra_deg, dec_deg)
        n = len(ra_deg)
        names = np.full(n, "", dtype=object)
        hemis = np.full(n, "", dtype=object)
        ok = np.zeros(n, dtype=bool)
        # South first, then north overrides where the brick centre lies
        # above the combination line, so each position ends with its
        # DR9-primary reduction.
        for hemisphere in ("south", "north"):
            idx, inside = self._covered(points, ra_deg, dec_deg, hemisphere,
                                        min_nexp)
            if hemisphere == "north":
                north_dec = self.tables["north"]["dec"][idx]
                inside &= (north_dec >= NORTH_SOUTH_DEC) | ~ok
            names[inside] = self.tables[hemisphere]["brickname"][idx[inside]]
            hemis[inside] = hemisphere
            ok |= inside
        return names, hemis, ok
