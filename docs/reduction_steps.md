(reduction_steps)=

# Reduction steps

Some sections of the configuration file are dedicated to specific steps of the
reduction process, and contain important parameters that can ultimately affect
the final transmission spectrum. All the relevant parameters can be specified by
the user. These sections can be duplicated under a specific dataset or
instrument if a different treatment is required for a given dataset.

### refraction

In this section/subsection the user can decide the ``method`` for the
**differential refraction** modeling, whether with a ``polynomial`` or with a
``spline``. A series of keywords (e.g., ``fit_order``, ``knots_spacing``) help the user to
set the parameters of the modeling.
In addition, the user can also decide the ``approach`` to compute the differential refraction correction, whether over an *individual_order*, or over the *full_spectrum* at once, iteratively (the number of iterations should be set) and with a sigma-clipping removal of the outliers. 

### CLV_RM_correction

Here some parameters related to the correction of the center-to-limb variation
and the Rossiter-McLaughlin effect are listed:

- ``synthesis_files``: the path + the root name of the synthesis files (i.e., the stellar model, the list with the limb-angle values, and the synthetic stellar spectra at each limb angle);
- ``n_gridpoints``: the number of cell per side to sample the stellar surface (it must be odd);
- ``time_step``: the oversampling step for the computation of the CLV+RM
  correction (in seconds);
- ``rebinning_step``: the rebinning step used during the synthesis of the surface element;
- ``continuum_range``: the wavelength region for the computation of the continuum level;
- ``radius_factor``: interval to set up the grid of r values for the rescaling
  factor of the planetary radius; the user must indicate the start and the end
  of the interval, and the spacing between values (e.g., *[0.5, 2.5, 0.1]*).
