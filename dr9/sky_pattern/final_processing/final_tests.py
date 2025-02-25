# Test to run for skyscale_ccds

# 1. Output should be zero
mask = ccd_new['run']==-1
np.max(np.abs(ccd_new['skyscale'][mask]))

# 2. Output should be nonzero
mask = ccd_new['run']>=0
np.min(np.abs(ccd_new['skyscale'][mask]))

# 3
raw = Table.read('/global/cscratch1/sd/rongpu/dr9dev/sky_pattern/sky_scales/skyscales_ccds_raw.fits')
