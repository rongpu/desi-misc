from __future__ import division, print_function
import sys, os, glob, gc, warnings
import numpy as np
from astropy.table import Table, vstack, hstack
import fitsio
import matplotlib.pyplot as plt

import numpy.ma
from matplotlib.collections import PolyCollection
from matplotlib.colors import Normalize, colorConverter
from matplotlib.patches import Ellipse
from types import MethodType
from astropy.coordinates import SkyCoord, HeliocentricTrueEcliptic, ICRS, Longitude
import astropy.units as u


class MaskedArrayWithLimits(numpy.ma.MaskedArray):
    """Masked array with additional `vmin`, `vmax` attributes.

    This class accepts the same arguments as
    :class:`~numpy.ma.MaskedArray`.

    This is not a general-purpose subclass and is only intended to simplify
    passing `vmin`, `vmax` limits from :func:`~desiutil.plots.prepare_data` to
    the plotting utility methods defined in this module.

    Parameters
    ----------
    vmin : :class:`float`, optional
        Minimum value when used for clipping or masking.
    vmax : :class:`float`, optional
        Maximum value when used for clipping or masking.

    Attributes
    ----------
    vmin : :class:`float`
        Minimum value when used for clipping or masking.
    vmax : :class:`float`
        Maximum value when used for clipping or masking.
    """
    def __new__(cls, *args, **kwargs):
        obj = super(MaskedArrayWithLimits, cls).__new__(cls, *args, **kwargs)
        if 'vmin' in kwargs:
            obj._optinfo['vmin'] = kwargs['vmin']
        #     obj.vmin = kwargs['vmin']
        # else:
        #     obj.vmin = None
        if 'vmax' in kwargs:
            obj._optinfo['vmax'] = kwargs['vmax']
        #     obj.vmax = kwargs['vmax']
        # else:
        #     obj.vmax = None
        return obj

    @property
    def vmin(self):
        return self._optinfo.get('vmin', None)

    @property
    def vmax(self):
        return self._optinfo.get('vmax', None)


def prepare_data(data, mask=None, clip_lo=None, clip_hi=None,
                 save_limits=False):
    """Prepare array data for color mapping.

    Data is clipped and masked to be suitable for passing to matplotlib
    routines that automatically assign colors based on input values.

    Parameters
    ----------
    data : array or masked array
        Array of data values to assign colors for.
    mask : array of bool or None
        Array of bools with same shape as data, where True values indicate
        values that should be ignored when assigning colors.  When None, the
        mask of a masked array will be used or all values of an unmasked
        array will be used.
    clip_lo : float or str
        Data values below clip_lo will be clipped to the minimum color. If
        clip_lo is a string, it should end with "%" and specify a percentile
        of un-masked data to clip below.
    clip_hi : float or str
        Data values above clip_hi will be clipped to the maximum color. If
        clip_hi is a string, it should end with "%" and specify a percentile
        of un-masked data to clip above.
    save_limits : bool
        Save the calculated lo/hi clip values as attributes vmin, vmax of
        the returned masked array.  Use this flag to indicate that plotting
        functions should use these vmin, vmax values when mapping the
        returned data to colors.

    Returns
    -------
    masked array
        Masked numpy array with the same shape as the input data, with any
        input mask applied (or copied from an input masked array) and values
        clipped to [clip_lo, clip_hi].

    Examples
    --------
    If no optional parameters are specified, the input data is returned
    with only non-finite values masked:

    >>> data = np.arange(5.)
    >>> prepare_data(data)
    masked_array(data = [0.0 1.0 2.0 3.0 4.0],
                 mask = [False False False False False],
           fill_value = 1e+20)
    <BLANKLINE>

    Any mask selection is propagated to the output:

    >>> prepare_data(data, data == 2)
    masked_array(data = [0.0 1.0 -- 3.0 4.0],
                 mask = [False False  True False False],
           fill_value = 1e+20)
    <BLANKLINE>

    Values can be clipped by specifying any combination of percentiles
    (specified as strings ending with "%") and numeric values:

    >>> prepare_data(data, clip_lo='25%', clip_hi=3.5)
    masked_array(data = [1.0 1.0 2.0 3.0 3.5],
                 mask = [False False False False False],
           fill_value = 1e+20)
    <BLANKLINE>

    Clipped values are also masked when the clip value or percentile
    is prefixed with "!":

    >>> prepare_data(data, clip_lo='!25%', clip_hi=3.5)
    masked_array(data = [-- 1.0 2.0 3.0 3.5],
                 mask = [ True False False False False],
           fill_value = 1e+20)
    <BLANKLINE>

    An input masked array is passed through without any copying unless
    clipping is requested:

    >>> masked = numpy.ma.arange(5)
    >>> masked is prepare_data(masked)
    True

    Use the save_limits option to store the clipping limits as vmin, vmax
    attributes of the returned object:

    >>> d = prepare_data(data, clip_lo=1, clip_hi=10, save_limits=True)
    >>> d.vmin, d.vmax
    (1.0, 10.0)

    These attributes can then be used by plotting routines to fix the input
    range used for colormapping, independently of the actual range of data.
    """
    data = np.asanyarray(data)
    if mask is None:
        try:
            # Use the mask associated with a MaskedArray.
            cmask = data.mask
            # If no clipping is requested, pass the input through.
            if clip_lo is None and clip_hi is None:
                return data
        except AttributeError:
            # Nothing is masked by default.
            cmask = np.zeros_like(data, dtype=bool)
    else:
        #
        # Make every effort to ensure that modifying the mask of the output
        # does not modify the input mask.
        #
        cmask = np.array(mask)
        if cmask.shape != data.shape:
            raise ValueError('Invalid mask shape.')
    # Mask any non-finite values.
    cmask |= ~np.isfinite(data)
    unmasked_data = data[~cmask]

    # Convert percentile clip values to absolute values.
    def get_clip(value):
        clip_mask = False
        if isinstance(value, str):
            if value.startswith('!'):
                clip_mask = True
                value = value[1:]
            if value.endswith('%'):
                value = np.percentile(unmasked_data, float(value[:-1]))
        return float(value), clip_mask

    if clip_lo is None:
        clip_lo, mask_lo = np.min(unmasked_data), False
    else:
        clip_lo, mask_lo = get_clip(clip_lo)
    if clip_hi is None:
        clip_hi, mask_hi = np.max(unmasked_data), False
    else:
        clip_hi, mask_hi = get_clip(clip_hi)

    if save_limits:
        clipped = MaskedArrayWithLimits(
            np.clip(data, clip_lo, clip_hi), cmask, vmin=clip_lo, vmax=clip_hi)
    else:
        clipped = numpy.ma.MaskedArray(
            np.clip(data, clip_lo, clip_hi), cmask)

    # Mask values outside the clip range, if requested.  The comparisons
    # below might trigger warnings for non-finite data.
    settings = np.seterr(all='ignore')
    if mask_lo:
        clipped.mask[data < clip_lo] = True
    if mask_hi:
        clipped.mask[data > clip_hi] = True
    np.seterr(**settings)

    return clipped


def init_sky(figsize=None, projection='mollweide', ra_center=120,
             galactic_plane_color='red', ecliptic_plane_color='red',
             ax=None):
    """Initialize matplotlib axes with a projection of the full sky.

    Parameters
    ----------
    projection : :class:`str`, optional
        Projection to use. Defaults to 'mollweide'.  To show the available projections,
        call :func:`matplotlib.projections.get_projection_names`.
    ra_center : :class:`float`, optional
        Projection is centered at this RA in degrees. Default is +120°, which avoids splitting
        the DESI northern and southern regions.
    galactic_plane_color : color name, optional
        Draw a solid curve representing the galactic plane using the specified color, or do
        nothing when ``None``.
    ecliptic_plane_color : color name, optional
        Draw a dotted curve representing the ecliptic plane using the specified color, or do
        nothing when ``None``.
    ax : :class:`~matplotlib.axes.Axes`, optional
        Axes to use for drawing this map, or create new axes if ``None``.

    Returns
    -------
    :class:`~matplotlib.axes.Axes`
        A matplotlib Axes object.  Helper methods ``projection_ra()`` and ``projection_dec()``
        are added to the object to facilitate conversion to projection coordinates.

    Notes
    -----
    If requested, the ecliptic and galactic planes are plotted with ``zorder`` set to 20.
    This keeps them above most other plotted objects, but legends should be set to
    a ``zorder`` higher than this value, for example::

        leg = ax.legend(ncol=2, loc=1)
        leg.set_zorder(25)
    """
    #
    # Internal functions.
    #
    def projection_ra(self, ra):
        r"""Shift `ra` to the origin of the Axes object and convert to radians.

        Parameters
        ----------
        ra : array-like
            Right Ascension in degrees.

        Returns
        -------
        array-like
            `ra` converted to plot coordinates.

        Notes
        -----
        In matplotlib, map projections expect longitude (RA), latitude (Dec)
        in radians with limits :math:`[-\pi, \pi]`, :math:`[-\pi/2, \pi/2]`,
        respectively.
        """
        #
        # Shift RA values.
        #
        r = np.remainder(ra + 360 - ra_center, 360)
        #
        # Scale conversion to [-180, 180].
        #
        r[r > 180] -= 360
        #
        # Reverse the scale: East to the left.
        #
        r = -r
        return np.radians(r)

    def projection_dec(self, dec):
        """Shift `dec` to the origin of the Axes object and convert to radians.

        Parameters
        ----------
        dec : array-like
            Declination in degrees.

        Returns
        -------
        array-like
            `dec` converted to plot coordinates.
        """
        return np.radians(dec)
    #
    # Create ax.
    #
    if ax is None:
        if figsize is None:
            figsize = (10.0, 5.0)
        fig = plt.figure(figsize=figsize, dpi=100)
        ax = plt.subplot(111, projection=projection)
    #
    # Prepare labels.
    #
    base_tick_labels = np.array([150, 120, 90, 60, 30, 0, 330, 300, 270, 240, 210])
    base_tick_labels = np.remainder(base_tick_labels+360+ra_center, 360)
    tick_labels = np.array(['{0}°'.format(l) for l in base_tick_labels])
    #
    # Galactic plane.
    #
    if galactic_plane_color is not None:
        galactic_l = np.linspace(0, 2 * np.pi, 1000)
        galactic = SkyCoord(l=galactic_l*u.radian, b=np.zeros_like(galactic_l)*u.radian,
                            frame='galactic').transform_to(ICRS)
        #
        # Project to map coordinates and display.  Use a scatter plot to
        # avoid wrap-around complications.
        #
        paths = ax.scatter(projection_ra(0, galactic.ra.degree),
                           projection_dec(0, galactic.dec.degree),
                           marker='.', s=20, lw=0, alpha=0.75,
                           c=galactic_plane_color, zorder=20)
        # Make sure the galactic plane stays above other displayed objects.
        # paths.set_zorder(20)
    #
    # Ecliptic plane.
    #
    if ecliptic_plane_color is not None:
        ecliptic_l = np.linspace(0, 2 * np.pi, 50)
        ecliptic = SkyCoord(lon=ecliptic_l*u.radian, lat=np.zeros_like(ecliptic_l)*u.radian, distance=1 * u.Mpc,
                            frame='heliocentrictrueecliptic').transform_to(ICRS)
        #
        # Project to map coordinates and display.  Use a scatter plot to
        # avoid wrap-around complications.
        #
        paths = ax.scatter(projection_ra(0, ecliptic.ra.degree),
                           projection_dec(0, ecliptic.dec.degree),
                           marker='.', s=20, lw=0, alpha=0.75,
                           c=ecliptic_plane_color, zorder=20)
        # paths.set_zorder(20)
    #
    # Set RA labels.
    #
    labels = ax.get_xticklabels()
    for l, item in enumerate(labels):
        item.set_text(tick_labels[l])
    ax.set_xticklabels(labels)
    #
    # Set axis labels.
    #
    ax.set_xlabel('R.A. [deg]')
    # ax.xaxis.label.set_fontsize(12)
    ax.set_ylabel('Dec. [deg]')
    # ax.yaxis.label.set_fontsize(12)
    ax.grid(True)
    #
    # Attach helper methods.
    #
    if hasattr(ax, '_ra_center'):
        warnings.warn("Attribute '_ra_center' detected.  Will be overwritten!")
    ax._ra_center = ra_center
    if hasattr(ax, 'projection_ra'):
        warnings.warn("Attribute 'projection_ra' detected.  Will be overwritten!")
    ax.projection_ra = MethodType(projection_ra, ax)
    if hasattr(ax, 'projection_dec'):
        warnings.warn("Attribute 'projection_dec' detected.  Will be overwritten!")
    ax.projection_dec = MethodType(projection_dec, ax)
    return ax


def plot_healpix_map(data, nest=False, cmap='viridis', colorbar=True,
                     label=None, ax=None, figsize=None, **kwargs):
    """Plot a healpix map using an all-sky projection.

    Pass the data array through :func:`prepare_data` to select a subset to plot
    and clip the color map to specified values or percentiles.

    This function is similar to :func:`plot_grid_map` but is generally slower
    at high resolution and has less elegant handling of pixels that wrap around
    in RA, which are not drawn.

    Requires that matplotlib and healpy are installed.

    Additional keyword parameters will be passed to :func:`init_sky`.

    Parameters
    ----------
    data : array or masked array
        1D array of data associated with each healpix.  Must have a size that
        exactly matches the number of pixels for some NSIDE value. Use the
        output of :func:`prepare_data` as a convenient way to specify
        data cuts and color map clipping.
    nest : :class:`bool`, optional
        If ``True``, assume NESTED pixel ordering.  Otheriwse, assume RING pixel
        ordering.
    cmap : colormap name or object, optional
        Matplotlib colormap to use for mapping data values to colors.
    colorbar : :class:`bool`, optional
        Draw a colorbar below the map when ``True``.
    label : :class:`str`, optional
        Label to display under the colorbar.  Ignored unless colorbar is ``True``.
    ax : :class:`~matplotlib.axes.Axes`, optional
        Axes to use for drawing this map, or create default axes using
        :func:`init_sky` when ``None``.

    Returns
    -------
    :class:`~matplotlib.axes.Axes`
        The axis object used for the plot.
    """
    import healpy as hp

    data = prepare_data(data)
    if len(data.shape) != 1:
        raise ValueError('Invalid data array, should be 1D.')
    nside = hp.npix2nside(len(data))
    #
    # Create axes.
    #
    if ax is None:
        ax = init_sky(figsize=figsize, **kwargs)
    proj_edge = ax._ra_center - 180
    #
    # Find the projection edge.
    #
    while proj_edge < 0:
        proj_edge += 360
    #
    # Get pixel boundaries as quadrilaterals.
    #
    corners = hp.boundaries(nside, np.arange(len(data)), step=1, nest=nest)
    corner_theta, corner_phi = hp.vec2ang(corners.transpose(0, 2, 1))
    corner_ra, corner_dec = (np.degrees(corner_phi),
                             np.degrees(np.pi/2-corner_theta))
    #
    # Convert sky coords to map coords.
    #
    x, y = ax.projection_ra(corner_ra), ax.projection_dec(corner_dec)
    #
    # Regroup into pixel corners.
    #
    verts = np.array([x.reshape(-1, 4), y.reshape(-1, 4)]).transpose(1, 2, 0)
    #
    # Find and mask any pixels that wrap around in RA.
    #
    uv_verts = np.array([corner_phi.reshape(-1, 4),
                         corner_theta.reshape(-1, 4)]).transpose(1, 2, 0)
    theta_edge = np.unique(uv_verts[:, :, 1])
    phi_edge = np.radians(proj_edge)
    eps = 0.1 * np.sqrt(hp.nside2pixarea(nside))
    wrapped1 = hp.ang2pix(nside, theta_edge, phi_edge - eps, nest=nest)
    wrapped2 = hp.ang2pix(nside, theta_edge, phi_edge + eps, nest=nest)
    wrapped = np.unique(np.hstack((wrapped1, wrapped2)))
    data.mask[wrapped] = True
    #
    # Normalize the data using its vmin, vmax attributes, if present.
    #
    try:
        norm = Normalize(vmin=data.vmin, vmax=data.vmax)
    except AttributeError:
        norm = None
    #
    # Make the collection and add it to the plot.
    #
    collection = PolyCollection(verts, array=data, cmap=cmap, norm=norm,
                                edgecolors='none')
    ax.add_collection(collection)
    ax.autoscale_view()

    if colorbar:
        bar = plt.colorbar(collection, ax=ax,
                           orientation='horizontal', spacing='proportional',
                           pad=0.11, fraction=0.05, aspect=50)
        if label:
            bar.set_label(label)

    return ax


def plot_sky_binned(nside, ra, dec, weights=None, data=None, figsize=(10.0, 5.0),
                    clip_lo=None, clip_hi=None, verbose=False,
                    cmap='viridis', colorbar=True, label=None, ax=None,
                    return_grid_data=False, **kwargs):
    """Show objects on the sky using a binned plot.

    Bin values either show object counts per unit sky area or, if an array
    of associated data values is provided, mean data values within each bin.
    Objects can have associated weights.

    Requires that matplotlib is installed. When plot_type is
    "healpix", healpy must also be installed.

    Additional keyword parameters will be passed to :func:`init_sky`.

    Parameters
    ----------
    nside: healpy nside
    ra : array
        Array of object RA values in degrees. Must have the same shape as
        dec and will be flattened if necessary.
    dec : array
        Array of object Dec values in degrees. Must have the same shape as
        ra and will be flattened if necessary.
    weights : array, optional
        Optional of weights associated with each object.  All objects are
        assumed to have equal weight when this is None.
    data : array, optional
        Optional array of scalar values associated with each object. The
        resulting plot shows the mean data value per bin when data is
        specified.  Otherwise, the plot shows counts per unit sky area.
    clip_lo : :class:`float` or :class:`str`, optional
        Clipping is applied to the plot data calculated as counts / area
        or the mean data value per bin. See :func:`prepare_data` for
        details.
    clip_hi : :class:`float` or :class:`str`, optional
        Clipping is applied to the plot data calculated as counts / area
        or the mean data value per bin. See :func:`prepare_data` for
        details.
    verbose : :class:`bool`, optional
        Print information about the automatic bin size calculation.
    cmap : colormap name or object, optional
        Matplotlib colormap to use for mapping data values to colors.
    colorbar : :class:`bool`, optional
        Draw a colorbar below the map when True.
    label : :class:`str`, optional
        Label to display under the colorbar.  Ignored unless colorbar is ``True``.
    ax : :class:`~matplotlib.axes.Axes`, optional
        Axes to use for drawing this map, or create default axes using
        :func:`init_sky` when ``None``.
    return_grid_data : :class:`bool`, optional
        If ``True``, return (ax, grid_data) instead of just ax.

    Returns
    -------
    :class:`~matplotlib.axes.Axes` or (ax, grid_data)
        The axis object used for the plot, and the grid_data if
        `return_grid_data` is ``True``.
    """
    ra = np.asarray(ra).reshape(-1)
    dec = np.asarray(dec).reshape(-1)
    if len(ra) != len(dec):
        raise ValueError('Arrays ra,dec must have same size.')

    if data is not None and weights is None:
        weights = np.ones_like(data)

    import healpy as hp

    bin_area = hp.nside2pixarea(nside, degrees=True)
    npix = hp.nside2npix(nside)
    nest = False
    if verbose:
        print('Using healpix map with NSIDE={0}'.format(nside),
              'and pixel area {:.3f} sq.deg.'.format(bin_area))

    pixels = hp.ang2pix(nside, np.radians(90 - dec), np.radians(ra), nest)
    counts = np.bincount(pixels, weights=weights, minlength=npix)
    if data is None:
        grid_data = counts / bin_area
    else:
        sums = np.bincount(pixels, weights=weights * data, minlength=npix)
        grid_data = np.zeros_like(sums, dtype=float)
        nonzero = counts > 0
        grid_data[nonzero] = sums[nonzero] / counts[nonzero]

    grid_data = prepare_data(grid_data, clip_lo=clip_lo, clip_hi=clip_hi)

    ax = plot_healpix_map(grid_data, nest=nest,
                          cmap=cmap, colorbar=colorbar, label=label,
                          ax=ax, figsize=figsize, **kwargs)

    if return_grid_data:
        return (ax, grid_data)
    else:
        return ax


