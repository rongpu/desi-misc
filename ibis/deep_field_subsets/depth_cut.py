# Optimize CCD selection for uniform depth coverage in IBIS deep fields
#
# This script iteratively removes CCDs to achieve uniform effective exposure time
# across the field, using a batch removal algorithm for efficiency.
#
# Process:
# 1. Load survey-ccds table and filter for field/band
# 2. Create spatial grid over field (0.025° spacing)
# 3. Calculate effective exposure time for each CCD
# 4. Iteratively optimize CCD selection:
#    - Compute current depth at each grid point
#    - For each CCD, count grid points with excess (>110% nominal) vs deficit (<90% nominal)
#    - Score CCDs: excess_count - 2*deficit_count - 5*severe_deficit_count
#    - Identify candidates where score > 0
#    - Greedily select non-overlapping candidates and remove them as a batch
#    - Repeat until no candidates remain
# 5. Output filtered CCD table and depth grid
#
# Usage:
#   python depth_cut.py <field_index> <band>
#   field_index: 0=COSMOS, 1=XMM-LSS
#   band: M411, M438, M464, M490, M517
#
# Example:
#   python depth_cut.py 0 M411
# or run in parallel:
#   salloc -N 1 -C cpu -t 04:00:00 -q interactive
#   parallel --jobs 10 < depth_cut_tasks.txt ; exit


from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio

import argparse
from multiprocessing import Pool

sys.path.append(os.path.expanduser('~/git/Python/user_modules/'))
from match_coord import match_coord


half_width = 4094/2 * 0.262 / 3600
half_height = 2046/2 * 0.262 / 3600

field_names = ['COSMOS', 'XMM-LSS']
field_centers = {'COSMOS':[150.1, 2.182], 'XMM-LSS':[35.75, -4.75]}
field_size = 5.0  # degrees
# grid_size = 0.025  # degrees
n_points = int(2e5)

n_processes = 12

# nominal efftime for a depth of 25.0
nominal_efftime_dict = {'M411': 1278.1, 'M438': 899.1, 'M464': 786.3, 'M490': 716.1, 'M517': 696.4}

parser = argparse.ArgumentParser()
parser.add_argument('field_index')
parser.add_argument('band')
args = parser.parse_args()
field_index = int(args.field_index)
band = args.band

field = field_names[field_index]
nominal_efftime = nominal_efftime_dict[band]

# field = 'COSMOS'
# band = 'M490'
# nominal_efftime = nominal_efftime_dict[band]

# Thresholds for excess/deficit
excess_threshold = 1.1 * nominal_efftime
deficit_threshold = 0.9 * nominal_efftime
deficit_threshold_1 = 0.6 * nominal_efftime  # more penalty for big deficits

ra_center, dec_center = field_centers[field]
ramin, ramax = ra_center - field_size/2/np.cos(np.radians(dec_center)), ra_center + field_size/2/np.cos(np.radians(dec_center))
decmin, decmax = dec_center - field_size/2, dec_center + field_size/2
print(field, ramin, ramax, decmin, decmax)
# ra_list = np.arange(ramin, ramax, grid_size)
# dec_list = np.arange(decmin, decmax, grid_size)
# ra_grid, dec_grid = np.meshgrid(ra_list, dec_list)

def sample_ra_dec(ramin, ramax, decmin, decmax, n, seed=None):
    rng = np.random.default_rng(seed)
    ra = rng.uniform(ramin, ramax, n)
    # Uniform in solid angle requires uniform sampling in sin(dec)
    sin_min = np.sin(np.radians(decmin))
    sin_max = np.sin(np.radians(decmax))
    dec = np.degrees(np.arcsin(rng.uniform(sin_min, sin_max, n)))
    return ra, dec

ra_grid, dec_grid = sample_ra_dec(ramin, ramax, decmin, decmax, n=n_points)

grid = Table()
grid['ra'] = ra_grid.flatten()
grid['dec'] = dec_grid.flatten()

surveyccd_path = '/global/cfs/cdirs/desicollab/users/rongpu/data/ibis/deep_fields/survey-ccds-ibis-dr1.fits'
surveyccd_psfex_path = '/global/cfs/cdirs/desicollab/users/rongpu/data/ibis/deep_fields/survey-ccds-ibis-dr1-psfex.fits'

ccd = Table(fitsio.read(surveyccd_path))
ccd1 = Table(fitsio.read(surveyccd_psfex_path))
assert len(ccd)==len(ccd1) and np.all(ccd['expnum']==ccd1['expnum'])
for col in ccd1.colnames:
    if col in ccd.colnames:
        ccd1.remove_column(col)
ccd = hstack([ccd, ccd1])
print(len(ccd))

mask = (ccd['ccd_cuts']==0) & (~ccd['failure'])
ccd = ccd[mask]
print(len(ccd))

mask = ccd['filter']==band
ccd = ccd[mask]
print('ccd', len(ccd))

mask = (ccd['ra']>ramin) & (ccd['ra']<ramax) & (ccd['dec']>decmin) & (ccd['dec']<decmax)
ccd = ccd[mask]
print('ccd', len(ccd))

ccd['median_ccdskycounts'], ccd['median_psf_fwhm'] = 0., 0.
tmp = Table()
tmp['expnum'], idx, tmp['count'] = np.unique(ccd['expnum'], return_counts=True, return_index=True)
for expnum in tmp['expnum']:
    mask = ccd['expnum']==expnum
    ccd['median_ccdskycounts'][mask] = np.median(ccd['ccdskycounts'][mask])
    ccd['median_psf_fwhm'][mask] = np.median(ccd['psf_fwhm'][mask])

ccd['efftime'] = 10**(0.4*ccd['zpt']-9) * ccd['exptime'] / (ccd['median_ccdskycounts'] * ccd['psf_fwhm']**2)

#######################################################################################################

def get_nccd(ccd_index):
    """Get mask of grid points covered by this CCD"""
    ccd_mask = (grid['ra']>ccd['ra'][ccd_index]-half_width/np.cos(np.radians(grid['dec']))) & (grid['ra']<ccd['ra'][ccd_index]+half_width/np.cos(np.radians(grid['dec']))) \
        & (grid['dec']>ccd['dec'][ccd_index]-half_height) & (grid['dec']<ccd['dec'][ccd_index]+half_height)
    return ccd_mask

def ccds_overlap(ccd_idx1, ccd_idx2):
    """Check if two CCDs have overlapping footprints"""
    ra1, dec1 = ccd['ra'][ccd_idx1], ccd['dec'][ccd_idx1]
    ra2, dec2 = ccd['ra'][ccd_idx2], ccd['dec'][ccd_idx2]

    # Compute half-widths accounting for declination
    hw1 = half_width / np.cos(np.radians(dec1))
    hw2 = half_width / np.cos(np.radians(dec2))

    # Handle RA wraparound at 0/360 degrees
    ra_diff = abs(ra1 - ra2)
    if ra_diff > 180:
        ra_diff = 360 - ra_diff

    # Check for overlap: centers must be within combined half-widths
    ra_overlap = ra_diff <= (hw1 + hw2)
    dec_overlap = abs(dec1 - dec2) <= (half_height + half_height)

    return ra_overlap and dec_overlap

ccds_to_keep = np.arange(len(ccd))
iteration = 0

while True:
    iteration += 1
    print(f"\n=== Iteration {iteration} ===")

    # Compute current coverage
    with Pool(processes=n_processes) as pool:
        ccd_mask_cube = np.array(pool.map(get_nccd, ccds_to_keep))
    grid_efftime = np.dot(ccd_mask_cube.T, ccd['efftime'][ccds_to_keep])

    print(f"CCDs remaining: {len(ccds_to_keep)}")

    # Find candidate CCDs to remove (excess > deficit)
    candidates = []
    for i, ccd_idx in enumerate(ccds_to_keep):
        # Get grid points covered by this CCD
        covered_mask = ccd_mask_cube[i]
        grid_efftime_at_ccd = grid_efftime[covered_mask]

        # Count excess and deficit points
        excess_count = np.sum(grid_efftime_at_ccd > excess_threshold)
        deficit_count = np.sum((grid_efftime_at_ccd < deficit_threshold) & (grid_efftime_at_ccd >= deficit_threshold_1))
        deficit_count_1 = np.sum(grid_efftime_at_ccd < deficit_threshold_1)

        score = excess_count - 1 * deficit_count - 5 * deficit_count_1

        if score > 0:
            candidates.append((i, score))

    print(f"Found {len(candidates)} candidate CCDs with excess > deficit")

    if not candidates:
        print("No candidates to remove, stopping")
        break

    # Sort candidates by score in descending order
    # This prioritizes removing CCDs that contribute most to excess
    candidates.sort(key=lambda x: x[1], reverse=True)

    # Greedy selection of non-overlapping CCDs to remove
    # Only remove CCDs whose score remains non-negative after removal
    to_remove_indices = []  # indices in ccds_to_keep
    to_remove_ccd_idx = []  # actual CCD indices for overlap checking

    for i, score in candidates:
        ccd_idx = ccds_to_keep[i]

        # Check if this CCD overlaps with any already selected for removal
        overlaps = False
        for other_ccd_idx in to_remove_ccd_idx:
            if ccds_overlap(ccd_idx, other_ccd_idx):
                overlaps = True
                break

        if not overlaps:
            # Compute grid after removing this CCD
            test_grid_efftime = grid_efftime - ccd['efftime'][ccd_idx] * ccd_mask_cube[i]

            # Recompute score for this CCD's footprint in the new grid
            covered_mask = ccd_mask_cube[i]
            grid_efftime_at_ccd = test_grid_efftime[covered_mask]

            excess_count = np.sum(grid_efftime_at_ccd > excess_threshold)
            deficit_count = np.sum((grid_efftime_at_ccd < deficit_threshold) & (grid_efftime_at_ccd >= deficit_threshold_1))
            deficit_count_1 = np.sum(grid_efftime_at_ccd < deficit_threshold_1)

            new_score = excess_count - 1 * deficit_count - 5 * deficit_count_1

            # Only remove if score remains non-negative
            if new_score >= 0:
                to_remove_indices.append(i)
                to_remove_ccd_idx.append(ccd_idx)
                # print(f"  Selected CCD {i} (new_score={new_score})")
            # else:
            #     print(f"  Skipped CCD {i} (new_score={new_score} < 0)")

    if not to_remove_indices:
        print("No non-overlapping candidates found, stopping")
        break

    print(f"Removing {len(to_remove_indices)} non-overlapping CCDs")
    ccds_to_keep = np.delete(ccds_to_keep, to_remove_indices)

# Final computation for output
with Pool(processes=n_processes) as pool:
    ccd_mask_cube = np.array(pool.map(get_nccd, ccds_to_keep))
grid_efftime = np.dot(ccd_mask_cube.T, ccd['efftime'][ccds_to_keep])

grid['nccd'] = np.sum(ccd_mask_cube, axis=0)
grid['efftime'] = grid_efftime
grid.write(f'/global/cfs/cdirs/desicollab/users/rongpu/data/ibis/deep_field_subsets/misc/nominal_depth_{field}_{band}.fits', overwrite=True)

ccd = ccd[ccds_to_keep]
ccd.write(f'/global/cfs/cdirs/desicollab/users/rongpu/data/ibis/deep_field_subsets/misc/survey-ccds-ibis-dr1-nominal-wide-depth_{field}_{band}.fits', overwrite=True)

print(f"\nFinal: {len(ccd)} CCDs remaining")
