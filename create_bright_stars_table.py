#!/usr/bin/env python3
"""
Create a table of bright stars (g mag < 15) from Gaia DR3 using WSDB.

Uses sqlutilpy to query and create a table in the workspace database.
Includes astrometry, parallax, proper motion, and photometry for bright stars.
Supports filtering by Galactic latitude and saving to FITS.
"""

import argparse
import pandas as pd
import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.table import Table

try:
    import sqlutilpy as sqlutil
except ImportError:
    print("Error: sqlutilpy not installed. Install with: pip install sqlutilpy")
    exit(1)


def create_bright_stars_query(mag_limit=15.0, ra_min=None, ra_max=None,
                             dec_min=None, dec_max=None):
    """
    Create SQL query for bright stars from Gaia DR3.

    Parameters:
    -----------
    mag_limit : float
        G-band magnitude limit (default: 15.0 for bright stars)
    ra_min, ra_max : float, optional
        RA range in degrees. If None, queries all-sky.
    dec_min, dec_max : float, optional
        Dec range in degrees. If None, queries all-sky.

    Returns:
    --------
    str : SQL query string
    """

    # Build WHERE clause
    where_clauses = [f"phot_g_mean_mag < {mag_limit}"]

    # Handle RA wrap-around at 0/360 degrees
    if ra_min is not None and ra_max is not None:
        if ra_min > ra_max:
            # Wrap-around case (e.g., 350-10 degrees)
            where_clauses.append(f"(ra >= {ra_min} OR ra <= {ra_max})")
        else:
            where_clauses.append(f"ra BETWEEN {ra_min} AND {ra_max}")

    if dec_min is not None and dec_max is not None:
        where_clauses.append(f"dec BETWEEN {dec_min} AND {dec_max}")

    where_clause = " AND ".join(where_clauses)

    query = f"""
    SELECT
        source_id,
        ra,
        dec,
        phot_g_mean_mag
    FROM gaia_dr3.gaia_source
    WHERE {where_clause}
    ORDER BY phot_g_mean_mag
    """

    return query


def execute_query(query, test_mode=False, max_rows=None):
    """
    Execute WSDB query and return results as DataFrame.

    Parameters:
    -----------
    query : str
        SQL query to execute
    test_mode : bool
        If True, print query before executing
    max_rows : int, optional
        Limit results to N rows (for testing)

    Returns:
    --------
    pd.DataFrame : Query results
    """
    if test_mode:
        print("Executing query...")
        print(f"Query length: {len(query)} characters\n")

    try:
        # Execute query and get results as dictionary
        data = sqlutil.get(query, asDict=True)

        # Convert to DataFrame
        df = pd.DataFrame(data)

        if max_rows:
            df = df.head(max_rows)

        print(f"✓ Retrieved {len(df)} stars")
        print(f"  G magnitude range: {df['phot_g_mean_mag'].min():.2f} - {df['phot_g_mean_mag'].max():.2f}")

        return df

    except Exception as e:
        print(f"✗ Query failed: {e}")
        return None


def filter_by_galactic_latitude(df, b_min=-90, b_max=90):
    """Filter stars by Galactic latitude.

    Parameters:
    -----------
    df : pd.DataFrame
        Must contain 'ra' and 'dec' columns
    b_min, b_max : float
        Galactic latitude range in degrees

    Returns:
    --------
    pd.DataFrame : Filtered subset
    """
    coords = SkyCoord(ra=df['ra'].values*u.deg, dec=df['dec'].values*u.deg, frame='icrs')
    galactic = coords.galactic

    mask = (np.abs(galactic.b.degree) >= b_min) & (np.abs(galactic.b.degree) <= b_max)
    df_filtered = df[mask].copy()

    # Add Galactic coordinates
    df_filtered['gal_l'] = galactic.l.degree[mask]
    df_filtered['gal_b'] = galactic.b.degree[mask]

    return df_filtered


def save_to_csv(df, filename):
    """Save DataFrame to CSV file."""
    df.to_csv(filename, index=False)
    print(f"✓ Saved {len(df)} stars to {filename}")


def save_to_fits(df, filename):
    """Save DataFrame to FITS file."""
    table = Table.from_pandas(df)
    table.write(filename, overwrite=True, format='fits')
    print(f"✓ Saved {len(df)} stars to {filename}")


def upload_to_wsdb(df, table_name):
    """
    Upload DataFrame as a table to WSDB.

    Parameters:
    -----------
    df : pd.DataFrame
        Data to upload
    table_name : str
        Name for the new table
    """
    try:
        # Convert DataFrame to dictionary for upload
        data_dict = {col: df[col].values for col in df.columns}

        # Upload to WSDB
        sqlutil.upload(table_name, data_dict)

        print(f"✓ Uploaded {len(df)} rows to WSDB table: {table_name}")

        # Verify upload
        result = sqlutil.get(f"SELECT COUNT(*) as count FROM {table_name}")
        count = result['count'][0] if isinstance(result['count'], np.ndarray) else result['count']
        print(f"✓ Verified: {count} rows in table")

    except Exception as e:
        print(f"✗ Upload failed: {e}")


def main():
    parser = argparse.ArgumentParser(
        description="Create a table of bright stars (g < 15) from Gaia DR3 using WSDB"
    )
    parser.add_argument('--mag-limit', type=float, default=15.0,
                        help='G-band magnitude limit (default: 15.0)')
    parser.add_argument('--ra-min', type=float, default=None,
                        help='RA minimum in degrees')
    parser.add_argument('--ra-max', type=float, default=None,
                        help='RA maximum in degrees')
    parser.add_argument('--dec-min', type=float, default=None,
                        help='Dec minimum in degrees')
    parser.add_argument('--dec-max', type=float, default=None,
                        help='Dec maximum in degrees')
    parser.add_argument('--gal-lat-min', type=float, default=5.0,
                        help='Minimum Galactic latitude |b| in degrees (default: 5.0)')
    parser.add_argument('--gal-lat-max', type=float, default=90.0,
                        help='Maximum Galactic latitude |b| in degrees (default: 90.0)')
    parser.add_argument('--test', action='store_true',
                        help='Test mode: print query before executing')
    parser.add_argument('--limit', type=int, default=None,
                        help='Limit to N results (for testing)')
    parser.add_argument('--output-csv', type=str, default=None,
                        help='Save results to CSV file')
    parser.add_argument('--output-fits', type=str, default=None,
                        help='Save results to FITS file')
    parser.add_argument('--upload-table', type=str, default=None,
                        help='Upload results to WSDB table with this name')

    args = parser.parse_args()

    print("=" * 70)
    print("Bright Stars (g < 15) Query from Gaia DR3")
    print("=" * 70)

    # Create query
    query = create_bright_stars_query(
        mag_limit=args.mag_limit,
        ra_min=args.ra_min,
        ra_max=args.ra_max,
        dec_min=args.dec_min,
        dec_max=args.dec_max
    )

    # Execute query
    df = execute_query(query, test_mode=args.test, max_rows=args.limit)

    if df is None or len(df) == 0:
        print("No results returned.")
        return

    n_before_gal_filter = len(df)

    # Apply Galactic latitude filter
    print(f"\nFiltering by Galactic latitude |b| > {args.gal_lat_min}°...")
    df = filter_by_galactic_latitude(df, b_min=args.gal_lat_min, b_max=args.gal_lat_max)
    n_after_gal_filter = len(df)

    print(f"  Before filter: {n_before_gal_filter}")
    print(f"  After filter:  {n_after_gal_filter}")
    print(f"  Removed: {n_before_gal_filter - n_after_gal_filter}")

    if len(df) == 0:
        print("No stars remain after Galactic latitude filter.")
        return

    # Print summary
    print("\n" + "=" * 70)
    print("Data Summary (after filtering)")
    print("=" * 70)
    print(f"Total stars: {len(df)}")
    print(f"\nBrightest stars (by g magnitude):")
    print(df[['source_id', 'ra', 'dec', 'phot_g_mean_mag']].head(10).to_string(index=False))

    print(f"\nPhotometric range:")
    print(f"  G mag: {df['phot_g_mean_mag'].min():.2f} - {df['phot_g_mean_mag'].max():.2f} mag")

    # Save to CSV if requested
    if args.output_csv:
        save_to_csv(df, args.output_csv)

    # Save to FITS if requested
    if args.output_fits:
        save_to_fits(df, args.output_fits)

    # Upload to WSDB if requested
    if args.upload_table:
        upload_to_wsdb(df, args.upload_table)

    print("\n" + "=" * 70)


if __name__ == '__main__':
    main()
