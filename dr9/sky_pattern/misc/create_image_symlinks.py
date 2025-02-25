import glob, os

# band = 'g'
for band in ['g', 'r', 'z']:
    
    fn_list = glob.glob('/Users/rongpu/Desktop/sky_pattern/sky_templates_v2/fit_scale/{}/*'.format(band))

    for fn in fn_list:
        fn1 = os.path.join('/Users/rongpu/Desktop/sky_pattern/sky_templates_v2/original_CP_images/{}/'.format(band), os.path.basename(fn.replace('_fitscale', '')))
        try:
            os.symlink(fn1, '/Users/rongpu/Desktop/sky_pattern/sky_templates_v2/combined/fit_scale_vs_original/{}/'.format(band)+os.path.basename(fn1))
        except:
            pass
