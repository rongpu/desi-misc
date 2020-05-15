**Steps**

1. Obtain DR8 blob mask for all CCD images ([code](https://github.com/rongpu/desi-misc/blob/master/sky_pattern/blobmask_taskfarmer/save_ccd_blob_mask.py))
2. Create a "sky run" file that groups exposures from one or more nights into one "run"; a sky template will be generated for each run (code?)
3. Compute the sky templates via median stacking of the blob-masked images ([code](https://github.com/rongpu/desi-misc/blob/master/sky_pattern/sky_template_taskfarmer/compute_sky_templates.py))
4. Visual inspection of all the templates
5. Obtain the sky coefficients for individual CCD images via template fitting ([code](https://github.com/rongpu/desi-misc/blob/master/sky_pattern/sky_fitting_taskfarmer/sky_template_fitting.py))
6. Compile a CCDs file with the coefficients ([code1](https://github.com/rongpu/desi-misc/blob/master/sky_pattern/final_processing/assemble_skyscales.py), [code2](https://github.com/rongpu/desi-misc/blob/master/sky_pattern/final_processing/create_final_skyscale_ccds.py))

Note: the actual sky subtraction uses the median of the coefficients of all CCDs in an exposure; we require that to use the median coefficient, at least 10 CCDs have their coefficients measured, otherwise use ccdskycounts as the coefficient.

**Additional fixes**

 - For templates that are particularly bad, identify the problematic exposures and blacklist them, and regenerate templates
 - For bad templates that are not salvageable, replace them with neighbors
 - For templates with bad S30, replace S30 with neighbors 
 - For exposures with only half of N10 usable, regenerate templates with the bad half masked

**QA images**

[Templates](https://portal.nersc.gov/project/desi/users/rongpu/plots/dr9dev/sky_pattern/sky_templates_v2/templates_final/) ([code](https://github.com/rongpu/desi-misc/blob/master/sky_pattern/create_images/create_images_templates.py))  
[Original CP images](https://portal.nersc.gov/project/desi/users/rongpu/plots/dr9dev/sky_pattern/sky_templates_v2/original_CP_images/) ([code](https://github.com/rongpu/desi-misc/blob/master/sky_pattern/create_images/create_images_cp_original.py))  
[Original CP images (DR8 version)](https://portal.nersc.gov/project/desi/users/rongpu/plots/dr9dev/sky_pattern/sky_templates_v2/original_CP_images_dr8/) ([code](https://github.com/rongpu/desi-misc/blob/master/sky_pattern/create_images/create_images_cp_original_dr8.py))  
[Sky-corrected images](https://portal.nersc.gov/project/desi/users/rongpu/plots/dr9dev/sky_pattern/sky_templates_v2/median_fit_scale/) ([code](https://github.com/rongpu/desi-misc/blob/master/sky_pattern/create_images/create_images_median_fitscale.py))  

