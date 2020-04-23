**Producing the fringe templates**

1. [save_ccd_blob_mask_sbatch.py](https://github.com/rongpu/desi-misc/blob/master/ccd_fringe/save_ccd_blob_mask_sbatch.py) - Map the blob mask images (in LS bricks) to individual CCD images; the per-CCD blob masks are saved on disk; they are used to mask known sources for obtaining the fringe images and the per-CCD fringe scales.
2. [compute_nightly_smooth_sky.py](https://github.com/rongpu/desi-misc/blob/master/ccd_fringe/compute_nightly_smooth_sky.py) - Compute the average per-night smooth sky image (at much larger scale than the fringe scale); this sky component will be subtracted from the image before computing the fringe image.
3. [compute_fringe_images_with_sky.py](https://github.com/rongpu/desi-misc/blob/master/ccd_fringe/compute_fringe_images.py) - Compute the raw (unnormalized and unsmoothed) fringe images via median stacking.
4. [smoothing_fringe_images.py](https://github.com/rongpu/desi-misc/blob/master/ccd_fringe/smoothing_fringe_images.py) and [convert_to_fits.py](https://github.com/rongpu/desi-misc/blob/master/ccd_fringe/convert_to_fits.py) - Convolve the raw fringe image with a Gaussian filter to produce a smoothed fringe image, and save as FITS files.
5. [normalize_fringe_templates.py](https://github.com/rongpu/desi-misc/blob/master/ccd_fringe/normalize_fringe_templates.py) - Normalize the fringe images so that a single multiplicative fringe scale can be used for all CCDs in an exposure. The fringe fitting results are used.

**Perform fringe fitting**

1. [fringe_template_fitting.py](https://github.com/rongpu/desi-misc/blob/master/ccd_fringe/fringe_template_fitting.py) Fit the (blank sky part of) images with fringe templates to obtain the per-CCD fringe scales.
2. [save_fringe_corrected_images.py](https://github.com/rongpu/desi-misc/blob/master/ccd_fringe/save_fringe_corrected_images.py) Save the fringe corrected images to disk.

**Plotting scripts**

[plot_fringe_on_focal_plane.py](https://github.com/rongpu/desi-misc/blob/master/ccd_fringe/plot_fringe_on_focal_plane.py)  
[plot_images_on_focal_plane.py](https://github.com/rongpu/desi-misc/blob/master/ccd_fringe/plot_images_on_focal_plane.py)  
(see also the Jupyter notebooks)
