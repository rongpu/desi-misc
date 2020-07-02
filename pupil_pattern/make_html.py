from __future__ import division, print_function
import sys, os, glob, time, warnings, gc

##########################################################################################

image_dir = '/global/cfs/cdirs/cosmo/staging/decam/CP-LG9'
image_path_list = glob.glob(os.path.join(image_dir, '*ooi*.fits.fz'))
image_path_list = sorted(image_path_list)

for band in ['g', 'r', 'z']:

    f = open("/global/cfs/cdirs/desi/users/rongpu/plots/dr9dev/pupil_pattern/{}-band-images.html".format(band), "w")
    f.write('<html>\n')
    f.write('<table>\n')
    f.write('<th>LG</th>\n')
    f.write('<th>LG pupil corrected</th> \n')
    f.write('<th>Regular CP</th>\n')
    f.write('<th>Regular CP sky corrected</th>\n')

    for image_path in image_path_list:

        if '_ooi_{}_'.format(band) in image_path:

            image_fn = os.path.basename(image_path).replace('.fits.fz', '.png')
            image_fn_cp = image_fn.replace('lg9', 'ls9')
            
            f.write('<tr>\n')

            f.write('<td><a href=\'lg_original/{}\'><img src=\'lg_original/{}\' width=\'400\'></a></td>\n'.format(image_fn, image_fn))
            f.write('<td><a href=\'lg_pupil_corrected/{}\'><img src=\'lg_pupil_corrected/{}\' width=\'400\'></a></td>\n'.format(image_fn, image_fn))
            f.write('<td><a href=\'lg_regular_cp/{}\'><img src=\'lg_regular_cp/{}\' width=\'400\'></a></td>\n'.format(image_fn_cp, image_fn_cp))
            f.write('<td><a href=\'lg_regular_cp-sky_corrected/{}\'><img src=\'lg_regular_cp-sky_corrected/{}\' width=\'400\'></a></td>\n'.format(image_fn_cp, image_fn_cp))
            
            # f.write('<td><a href=\'lg_regular_cp-sky_corrected/c4d_180221_072207_ooi_g_ls9.png\'><img src=\'lg_regular_cp-sky_corrected/c4d_180221_072207_ooi_g_ls9.png\' width=\'400\'></a></td>\n')
            
            f.write('</tr>\n')
            
    f.write('</table>\n')
    f.close()

##########################################################################################

f = open("/global/cfs/cdirs/desi/users/rongpu/plots/dr9dev/pupil_pattern/templates.html", "w")
f.write('<html>\n')
f.write('<table>\n')
f.write('<th>g-band template</th>\n')
f.write('<th>r-band template</th> \n')
f.write('<th>z-band template [cassette_3]</th>\n')
f.write('<th>z-band template [cassette_1]</th>\n')

f.write('<tr>\n')
image_fn = 'g_1_pupil_template_smooth.png'
f.write('<td><a href=\'pupil_templates/{}\'><img src=\'pupil_templates/{}\' width=\'400\'></a></td>\n'.format(image_fn, image_fn))
image_fn = 'r_2_pupil_template_smooth.png'
f.write('<td><a href=\'pupil_templates/{}\'><img src=\'pupil_templates/{}\' width=\'400\'></a></td>\n'.format(image_fn, image_fn))
image_fn = 'z_3_pupil_template_smooth.png'
f.write('<td><a href=\'pupil_templates/{}\'><img src=\'pupil_templates/{}\' width=\'400\'></a></td>\n'.format(image_fn, image_fn))
image_fn = 'z_4_pupil_template_smooth.png'
f.write('<td><a href=\'pupil_templates/{}\'><img src=\'pupil_templates/{}\' width=\'400\'></a></td>\n'.format(image_fn, image_fn))
f.write('</tr>\n')
f.write('</table>\n')
f.close()
