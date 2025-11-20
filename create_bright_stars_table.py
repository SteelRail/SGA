#!/usr/bin/env python3
"""
Query bright stars (g mag < 15) from Gaia DR3 with |b| galactic latitude filter.
"""

import argparse
import numpy as np
import pandas as pd
from astropy.table import Table

try:
    import sqlutilpy as sqlutil
except ImportError:
    print("Error: sqlutilpy not installed. Install with: pip install sqlutilpy")
    exit(1)


def create_bright_stars_query(mag_limit=15.0, gal_lat_min=5.0, ra_min=None, ra_max=None,
                             dec_min=None, dec_max=None):
    """
    Create SQL query for bright stars from Gaia DR3 with galactic latitude filter.

    Parameters:
    -----------
    mag_limit : float
        G-band magnitude limit (default: 15.0)
    gal_lat_min : float
        Minimum galactic latitude |b| in degrees (default: 5.0)
    ra_min, ra_max : float, optional
        RA range in degrees (for regional testing)
    dec_min, dec_max : float, optional
        Dec range in degrees (for regional testing)

    Returns:
    --------
    str : SQL query string
    """
    where_clauses = [f"phot_g_mean_mag < {mag_limit}", f"ABS(b) >= {gal_lat_min}"]

    if ra_min is not None and ra_max is not None:
        where_clauses.append(f"ra BETWEEN {ra_min} AND {ra_max}")

    if dec_min is not None and dec_max is not None:
        where_clauses.append(f"dec BETWEEN {dec_min} AND {dec_max}")

    where_clause = " AND ".join(where_clauses)

    query = f"""
    SELECT
        source_id,
        ra,
        dec,
        phot_g_mean_mag,
        b
    FROM gaia_dr3.gaia_source
    WHERE {where_clause}
    ORDER BY phot_g_mean_mag
    """
    return query


def execute_query(query):
    """
    Execute WSDB query and return results as DataFrame.

    Parameters:
    -----------
    query : str
        SQL query to execute

    Returns:
    --------
    pd.DataFrame : Query results or None on failure
    """
    try:
        data = sqlutil.get(query, asDict=True)
        df = pd.DataFrame(data)
        print(f"✓ Retrieved {len(df)} stars from Gaia DR3")
        return df
    except Exception as e:
        print(f"✗ Query failed: {e}")
        return None




def save_to_csv(df, filename):
    """Save DataFrame to CSV file."""
    df.to_csv(filename, index=False)
    print(f"✓ Saved {len(df)} stars to {filename}")


def save_to_fits(df, filename):
    """Save DataFrame to FITS file."""
    table = Table.from_pandas(df)
    table.write(filename, overwrite=True, format='fits')
    print(f"✓ Saved {len(df)} stars to {filename}")


def main():
    parser = argparse.ArgumentParser(
        description="Query bright stars (g < 15) from Gaia DR3 with |b| galactic latitude filter"
    )
    parser.add_argument('--mag-limit', type=float, default=15.0,
                        help='G-band magnitude limit (default: 15.0)')
    parser.add_argument('--gal-lat-min', type=float, default=5.0,
                        help='Minimum galactic latitude |b| in degrees (default: 5.0)')
    parser.add_argument('--ra-min', type=float, default=None,
                        help='RA minimum in degrees (for regional testing)')
    parser.add_argument('--ra-max', type=float, default=None,
                        help='RA maximum in degrees (for regional testing)')
    parser.add_argument('--dec-min', type=float, default=None,
                        help='Dec minimum in degrees (for regional testing)')
    parser.add_argument('--dec-max', type=float, default=None,
                        help='Dec maximum in degrees (for regional testing)')
    parser.add_argument('--output-csv', type=str, default=None,
                        help='Save results to CSV file')
    parser.add_argument('--output-fits', type=str, default='./data/catalogues/gaia_bright_stars_g15_b5.fits',
                        help='Save results to FITS file (default: ./data/catalogues/gaia_bright_stars_g15_b5.fits)\n(g15=g<15.0 mag, b5=|b|>=5.0 deg)')

    args = parser.parse_args()

    print("=" * 70)
    print(f"Bright Stars (g < {args.mag_limit}) with |b| >= {args.gal_lat_min}°")
    if args.ra_min is not None or args.dec_min is not None:
        print(f"Region: RA [{args.ra_min}, {args.ra_max}], Dec [{args.dec_min}, {args.dec_max}]")
    print("=" * 70)

    # Query Gaia DR3 (with |b| filter inside query)
    query = create_bright_stars_query(
        mag_limit=args.mag_limit,
        gal_lat_min=args.gal_lat_min,
        ra_min=args.ra_min,
        ra_max=args.ra_max,
        dec_min=args.dec_min,
        dec_max=args.dec_max
    )
    df = execute_query(query)

    if df is None or len(df) == 0:
        print("No results returned.")
        return

    print(f"Retrieved {len(df)} stars\n")

    # Print summary
    print(f"Brightest stars (by g magnitude):")
    print(df[['source_id', 'ra', 'dec', 'phot_g_mean_mag', 'b']].head(10).to_string(index=False))
    print(f"\nG magnitude range: {df['phot_g_mean_mag'].min():.2f} - {df['phot_g_mean_mag'].max():.2f}")
    print(f"Galactic latitude |b| range: {np.abs(df['b']).min():.2f}° - {np.abs(df['b']).max():.2f}°")

    # Save if requested
    if args.output_csv:
        save_to_csv(df, args.output_csv)

    if args.output_fits:
        save_to_fits(df, args.output_fits)

    print("=" * 70)


if __name__ == '__main__':
    main()
