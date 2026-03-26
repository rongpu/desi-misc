cd /global/cfs/cdirs/desicollab/users/rongpu/data/desi2/tertiary54/daily
ln -s /global/cfs/cdirs/desi/spectro/redux/daily/tiles/cumulative/83616/20260320/coadd-*.fits .
ln -s /global/cfs/cdirs/desi/spectro/redux/daily/tiles/cumulative/83617/20260321/coadd-*.fits .

source /global/common/software/desi/desi_environment.sh main
export RR_TEMPLATE_DIR=/global/cfs/cdirs/desicollab/users/rongpu/data/redrock/lae_templates_nmf_20260120
cd /global/cfs/cdirs/desicollab/users/rongpu/data/desi2/tertiary54/daily/rr_lae_2.25_3.4_nmf
salloc -N 4 -C cpu -t 04:00:00 -q interactive
parallel --jobs 4 --delay 1 < desi_redrock_lae_nmf.txt

