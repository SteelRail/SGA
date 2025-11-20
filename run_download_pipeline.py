#!/usr/bin/env python3
"""
Manager script for end-to-end galaxy and field cutout download and Gaia sidecar generation.

Pipeline:
  1. Download SGA galaxy central cutouts
  2. Check galaxy cutout integrity and re-download broken files
  3. Download background field cutouts
  4. Check field cutout integrity and re-download broken files
  5. Generate Gaia bright-star sidecars for all cutouts (galaxies + fields)
"""

import os
import sys
import subprocess
import argparse
from pathlib import Path


def run_command(cmd, description, capture_output=False):
    """
    Run a shell command and handle errors gracefully.

    Parameters:
    -----------
    cmd : list
        Command and arguments as list (for subprocess)
    description : str
        Human-readable description of what's running
    capture_output : bool
        If True, capture and return stdout instead of printing

    Returns:
    --------
    int or str
        Exit code (if capture_output=False) or stdout (if capture_output=True)
    """
    print(f"\n{'='*70}")
    print(f"Step: {description}")
    print(f"{'='*70}")
    print(f"Command: {' '.join(cmd)}\n")

    try:
        if capture_output:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            return result.stdout
        else:
            result = subprocess.run(cmd, check=True)
            return result.returncode
    except subprocess.CalledProcessError as e:
        print(f"✗ FAILED: {description}")
        print(f"  Exit code: {e.returncode}")
        sys.exit(1)
    except FileNotFoundError:
        print(f"✗ ERROR: Command not found. Check your PATH and venv activation.")
        sys.exit(1)


def count_cutouts(folder):
    """Count JPEG files in folder."""
    if not os.path.isdir(folder):
        return 0
    return len([f for f in os.listdir(folder) if f.lower().endswith(('.jpg', '.jpeg'))])


def main():
    parser = argparse.ArgumentParser(
        description="End-to-end pipeline: download galaxy + field cutouts → check integrity → retry → generate Gaia sidecars"
    )
    parser.add_argument('--n_galaxies', type=int, default=100,
                        help='Number of SGA galaxies to download (default: 100)')
    parser.add_argument('--n_fields', type=int, default=7,
                        help='Number of background field positions per galaxy (default: 7)')
    parser.add_argument('--cutouts_folder', type=str, default='./data/cutouts',
                        help='Parent folder for cutouts (will create galaxies/ and fields/ subfolders) (default: ./data/cutouts)')
    parser.add_argument('--catalog', type=str, default='./data/catalogues/SGA-2020.fits',
                        help='Path to SGA catalog (default: ./data/catalogues/SGA-2020.fits)')
    parser.add_argument('--gaia_catalog', type=str, default='./data/catalogues/gaia_bright_stars_g15_b5.fits',
                        help='Path to Gaia bright-star catalog (default: ./data/catalogues/gaia_bright_stars_g15_b5.fits)\n(g15=g<15.0 mag, b5=|b|>=5.0 deg)')
    parser.add_argument('--pixscale', type=float, default=0.262,
                        help='Pixel scale in arcsec/pixel (default: 0.262)')
    parser.add_argument('--size', type=int, default=-1,
                        help='Force image size in pixels, -1 for adaptive (default: -1)')
    parser.add_argument('--num_workers', type=int, default=4,
                        help='Number of parallel workers for download (default: 4)')
    parser.add_argument('--bricks_file', type=str, default='./data/catalogues/survey-bricks-dr9-south.fits',
                        help='Path to DECALS bricks file (south)')
    parser.add_argument('--bricks_file_north', type=str, default='./data/catalogues/survey-bricks-dr9-north.fits',
                        help='Path to DECALS bricks file (north)')
    parser.add_argument('--broken_size_threshold', type=int, default=500,
                        help='File size threshold (bytes) for identifying broken downloads (default: 500)')
    parser.add_argument('--skip_galaxies', action='store_true',
                        help='Skip galaxy download and retry')
    parser.add_argument('--skip_fields', action='store_true',
                        help='Skip field download and retry')
    parser.add_argument('--skip_sidecars', action='store_true',
                        help='Skip Gaia sidecar generation')

    args = parser.parse_args()

    # Validate catalog files exist
    for path, name in [
        (args.catalog, 'SGA catalog'),
        (args.gaia_catalog, 'Gaia catalog'),
        (args.bricks_file, 'DECALS bricks (south)'),
        (args.bricks_file_north, 'DECALS bricks (north)')
    ]:
        if not os.path.exists(path):
            print(f"✗ ERROR: {name} not found at {path}")
            sys.exit(1)

    # Create folder structure
    galaxies_folder = os.path.join(args.cutouts_folder, 'galaxies')
    fields_folder = os.path.join(args.cutouts_folder, 'fields')

    print("=" * 70)
    print("GALAXY + FIELD CUTOUT DOWNLOAD & GAIA SIDECAR PIPELINE")
    print("=" * 70)
    print(f"Configuration:")
    print(f"  N galaxies:        {args.n_galaxies}")
    print(f"  N fields/galaxy:   {args.n_fields}")
    print(f"  Galaxies folder:   {galaxies_folder}")
    print(f"  Fields folder:     {fields_folder}")
    print(f"  SGA catalog:       {args.catalog}")
    print(f"  Gaia catalog:      {args.gaia_catalog}")
    print(f"  Pixel scale:       {args.pixscale} arcsec/pixel")
    print(f"  Image size:        {args.size} (pixels)")
    print(f"  Workers:           {args.num_workers}")
    print(f"  Broken threshold:  {args.broken_size_threshold} bytes")
    print()

    # Create folder structure
    os.makedirs(galaxies_folder, exist_ok=True)
    os.makedirs(fields_folder, exist_ok=True)

    # =========================================================================
    # Step 1: Download galaxy cutouts
    # =========================================================================
    if not args.skip_galaxies:
        cmd_download_gal = [
            'python', 'load_sga.py',
            '--catalog', args.catalog,
            '--download_folder', galaxies_folder,
            '--n_limit', str(args.n_galaxies),
            '--size', str(args.size),
            '--pixscale', str(args.pixscale),
            '--num_workers', str(args.num_workers),
            '--mode', 'galaxy',
            '--bricks_file', args.bricks_file,
            '--bricks_file_north', args.bricks_file_north,
            '--broken_size_threshold', str(args.broken_size_threshold),
        ]
        run_command(cmd_download_gal, f"Download {args.n_galaxies} SGA galaxy cutouts")

        n_gal = count_cutouts(galaxies_folder)
        print(f"✓ Downloaded {n_gal} galaxy cutouts to {galaxies_folder}")

        # =========================================================================
        # Step 1b: Check galaxy cutout integrity and retry
        # =========================================================================
        cmd_retry_gal = [
            'python', 'load_sga.py',
            '--catalog', args.catalog,
            '--download_folder', galaxies_folder,
            '--n_limit', str(args.n_galaxies),
            '--size', str(args.size),
            '--pixscale', str(args.pixscale),
            '--num_workers', str(args.num_workers),
            '--mode', 'retry',
            '--broken_size_threshold', str(args.broken_size_threshold),
        ]
        run_command(cmd_retry_gal, "Re-download broken or missing galaxy cutouts")

        n_gal_intact = count_cutouts(galaxies_folder)
        print(f"✓ Galaxy integrity check complete: {n_gal_intact} cutouts available")

    # =========================================================================
    # Step 2: Download field cutouts
    # =========================================================================
    if not args.skip_fields:
        cmd_download_fields = [
            'python', 'load_sga.py',
            '--catalog', args.catalog,
            '--download_folder', fields_folder,
            '--n_limit', str(args.n_galaxies),
            '--n_away', str(args.n_fields),
            '--size', str(args.size),
            '--pixscale', str(args.pixscale),
            '--num_workers', str(args.num_workers),
            '--mode', 'fields',
            '--bricks_file', args.bricks_file,
            '--bricks_file_north', args.bricks_file_north,
            '--broken_size_threshold', str(args.broken_size_threshold),
        ]
        run_command(cmd_download_fields, f"Download {args.n_galaxies * args.n_fields} field cutouts ({args.n_fields} per galaxy)")

        n_fields = count_cutouts(fields_folder)
        print(f"✓ Downloaded {n_fields} field cutouts to {fields_folder}")

        # =========================================================================
        # Step 2b: Check field cutout integrity and retry
        # =========================================================================
        cmd_retry_fields = [
            'python', 'load_sga.py',
            '--catalog', args.catalog,
            '--download_folder', fields_folder,
            '--n_limit', str(args.n_galaxies),
            '--n_away', str(args.n_fields),
            '--size', str(args.size),
            '--pixscale', str(args.pixscale),
            '--num_workers', str(args.num_workers),
            '--mode', 'retry',
            '--broken_size_threshold', str(args.broken_size_threshold),
        ]
        run_command(cmd_retry_fields, "Re-download broken or missing field cutouts")

        n_fields_intact = count_cutouts(fields_folder)
        print(f"✓ Field integrity check complete: {n_fields_intact} cutouts available")

    # =========================================================================
    # Step 3: Generate Gaia bright-star sidecars for ALL cutouts (galaxies + fields)
    # =========================================================================
    if not args.skip_sidecars:
        # Generate sidecars for both galaxies and fields folders
        for folder_name, folder_path in [('galaxies', galaxies_folder), ('fields', fields_folder)]:
            if os.path.isdir(folder_path) and count_cutouts(folder_path) > 0:
                cmd_sidecars = [
                    'python', 'make_gaia_sidecars.py',
                    '--cutout_folder', folder_path,
                    '--gaia_catalog', args.gaia_catalog,
                    '--pixscale', str(args.pixscale),
                    '--overwrite',
                ]
                run_command(cmd_sidecars, f"Generate Gaia sidecars for {folder_name} cutouts")
                print(f"✓ Sidecar generation complete for {folder_name}")

    # =========================================================================
    # Summary
    # =========================================================================
    print("\n" + "=" * 70)
    print("PIPELINE COMPLETE")
    print("=" * 70)

    n_gal_final = count_cutouts(galaxies_folder)
    n_gal_sidecars = len([f for f in os.listdir(galaxies_folder) if f.endswith('_gaia.csv')]) if os.path.isdir(galaxies_folder) else 0

    n_fields_final = count_cutouts(fields_folder)
    n_fields_sidecars = len([f for f in os.listdir(fields_folder) if f.endswith('_gaia.csv')]) if os.path.isdir(fields_folder) else 0

    print(f"Galaxy cutouts:     {n_gal_final} JPEG files, {n_gal_sidecars} sidecars")
    print(f"Field cutouts:      {n_fields_final} JPEG files, {n_fields_sidecars} sidecars")
    print(f"Total cutouts:      {n_gal_final + n_fields_final} JPEG files")
    print(f"Total sidecars:     {n_gal_sidecars + n_fields_sidecars} files")
    print(f"Output location:    {args.cutouts_folder}/")
    print("=" * 70)


if __name__ == '__main__':
    main()
