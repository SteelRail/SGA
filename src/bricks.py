"""The DR9 survey-bricks tables: which brick covers a sky position."""

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits

# DR9 combines the two hemisphere reductions at this declination: south
# (DECam) below, north (BASS+MzLS) above.
NORTH_SOUTH_DEC = 32.375


class Bricks:
    """Both hemisphere tables, loaded once; one operation: position -> brick."""

    def __init__(self,
                 north="data/catalogues/survey-bricks-dr9-north.fits.gz",
                 south="data/catalogues/survey-bricks-dr9-south.fits.gz"):
        self.tables = {}
        self._coords = {}
        for hemisphere, path in (("north", north), ("south", south)):
            with fits.open(path) as hdul:
                table = hdul[1].data
            self.tables[hemisphere] = table
            self._coords[hemisphere] = SkyCoord(
                table["ra"] * u.deg, table["dec"] * u.deg
            )

    def _covered(self, coords, ra_deg, dec_deg, hemisphere, min_nexp):
        """Nearest-brick index and whether it bounds and covers each position.

        Bounds use ra1/ra2/dec1/dec2 with RA wrap-around; coverage requires
        `min_nexp` exposures in each of g, r and z.
        """
        table = self.tables[hemisphere]
        idx, _, _ = coords.match_to_catalog_sky(self._coords[hemisphere])
        ok = np.ones(len(idx), dtype=bool)
        for band in ("g", "r", "z"):
            ok &= table[f"nexp_{band}"][idx] >= min_nexp
        ra1, ra2 = table["ra1"][idx], table["ra2"][idx]
        width = (ra2 - ra1 + 360.0) % 360.0
        dra = (ra_deg - ra1 + 360.0) % 360.0
        ok &= dra <= width
        ok &= (dec_deg >= table["dec1"][idx]) & (dec_deg <= table["dec2"][idx])
        return idx, ok

    def resolve(self, ra_deg, dec_deg, min_nexp=1):
        """Brick name and hemisphere for each position.

        Returns (bricknames, hemispheres, ok): `ok` is False where no
        covering brick exists — the footprint test. Positions covered by
        both reductions follow the DR9 combination rule: south below
        Dec 32.375°, north above.
        """
        ra_deg = np.atleast_1d(np.asarray(ra_deg, dtype=float))
        dec_deg = np.atleast_1d(np.asarray(dec_deg, dtype=float))
        coords = SkyCoord(ra_deg * u.deg, dec_deg * u.deg)
        n = len(ra_deg)
        names = np.full(n, "", dtype=object)
        hemis = np.full(n, "", dtype=object)
        ok = np.zeros(n, dtype=bool)
        # South first, then north overrides above the combination line, so
        # each position ends with its DR9-primary reduction.
        for hemisphere in ("south", "north"):
            idx, inside = self._covered(coords, ra_deg, dec_deg, hemisphere,
                                        min_nexp)
            if hemisphere == "north":
                inside &= (dec_deg >= NORTH_SOUTH_DEC) | ~ok
            names[inside] = self.tables[hemisphere]["brickname"][idx[inside]]
            hemis[inside] = hemisphere
            ok |= inside
        return names, hemis, ok
