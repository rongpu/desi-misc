ln -s /global/cfs/cdirs/desi/spectro/redux/daily/tiles/cumulative/83579/20260116/coadd-*.fits .
ln -s /global/cfs/cdirs/desi/spectro/redux/daily/tiles/cumulative/83580/20260115/coadd-*.fits .
ln -s /global/cfs/cdirs/desi/spectro/redux/daily/tiles/cumulative/83581/20260116/coadd-*.fits .
ln -s /global/cfs/cdirs/desi/spectro/redux/daily/tiles/cumulative/83582/20260116/coadd-*.fits .

################################################################################################################

source /global/common/software/desi/desi_environment.sh main
export RR_TEMPLATE_DIR=/global/cfs/cdirs/desicollab/users/rongpu/data/redrock/laelbg_templates-galaxy_only

cd /global/cfs/cdirs/desicollab/users/rongpu/data/desi2/tertiary49
salloc -N 4 -C cpu -t 04:00:00 -q interactive
parallel --jobs 4 --delay 1 < desi_redrock_lae.txt
parallel --jobs 4 --delay 1 < desi_emline.txt

################################################################################################################

source /global/common/software/desi/desi_environment.sh main
export RR_TEMPLATE_DIR=/global/cfs/cdirs/desicollab/users/rongpu/data/redrock/laelbg_templates-galaxy_only
export RR_ARCHETYPE=/global/cfs/cdirs/desicollab/users/rongpu/data/redrock/lae_archetypes_20260117/rrarchetype-lae.fits

cd /global/cfs/cdirs/desicollab/users/rongpu/data/desi2/tertiary49
salloc -N 4 -C cpu -t 04:00:00 -q interactive
parallel --jobs 4 --delay 1 < desi_redrock_archetype.txt
parallel --jobs 4 --delay 1 < desi_emline_archetype.txt

################################################################################################################

source /global/common/software/desi/desi_environment.sh main
export RR_TEMPLATE_DIR=/global/cfs/cdirs/desicollab/users/rongpu/data/redrock/lae_templates_nmf_20260120

cd /global/cfs/cdirs/desicollab/users/rongpu/data/desi2/tertiary49/ignore_z_camera
salloc -N 4 -C cpu -t 04:00:00 -q interactive
parallel --jobs 4 --delay 1 < desi_redrock_lae_nmf.txt


################################################################################################################

source /global/common/software/desi/desi_environment.sh main
export RR_TEMPLATE_DIR=/global/cfs/cdirs/desicollab/users/rongpu/data/redrock/lae_templates_20260120

PATH=/global/u2/r/rongpu/moregit/redrock/bin:$PATH
PYTHONPATH=/global/u2/r/rongpu/moregit/redrock/py:$PYTHONPATH

cd /global/cfs/cdirs/desicollab/users/rongpu/data/desi2/tertiary49/ignore_z_camera
salloc -N 4 -C cpu -t 04:00:00 -q interactive
parallel --jobs 4 --delay 1 < desi_redrock_lae_pca.txt
