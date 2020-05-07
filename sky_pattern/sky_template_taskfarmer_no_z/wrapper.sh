source /project/projectdirs/desi/software/desi_environment.sh 19.12
export PYTHONPATH=/global/homes/r/rongpu/modules/lib/python3.6/site-packages/fitsio-1.1.0-py3.6-linux-x86_64.egg
cd $HOME/jobs/sky_template_taskfarmer_no_z
python compute_sky_templates_v2_no_z.py $1 $2
