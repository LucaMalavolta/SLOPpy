(configuration_file)=

# Configuration File

In order to successfully perform the extraction of a transmission spectrum,
SLOPpy requires a configuration file, written in ``YAML``, for each target.

The configuration file is divided in several sections:

#### *pipeline* and *plots* 


Here the user can list which analysis modules need to be executed (e.g., ``sky_correction``, ``master_out``) and which plots to be shown, respectively.
Each data reduction step is accompanied by a plotting module with the same name, so that checking the outcome of each reduction step is straightforward. It is important to remember that analyses and plots are performed independently, that is to say, the user can modify the list of plots anytime without the need of performing the time-consuming analysis again, even for intermediate steps. 

#### *nights* 

Assuming that during a night only one transit of a given target can be observed,
in this section the characteristics of each dataset gathered during a night are
detailed, such as the lists of ``all`` spectra, ``in_transit`` spectra,
``full_transit`` spectra and out-of-transit (``out_transit``) spectra, the
``instrument`` used to gather the observations (e.g, *HARPS*, *HARPS-N*), the kind of ``mask`` used for the calibration, and the time of mid-transit in BJD (time_of_transit).
If some residuals are still present after telluric correction,, the user can
decide to execute an additional correction by setting the keyword
``spline_residuals`` equal to *True*.

```{note} It is possible to analyzed data obtained during the same night by two independent instruments, however some night-specific parameters (as the time of transit) will be treated independently. To do so, just use two different *labels* for the nights.
```

### master-out 
In this section the user should set the ``wavelength_range`` and the
``wavelength_step`` (in Ångström), and the ``binning_factor`` (e.g., 20) of the
master-out. In addition, the user can decide whether or not to use the composite
master-out (e.g., ``use_composite``: *False*), which is obtained by combining
all the master-out spectra of the various nights analyzed.

### instruments
In this section the instrument properties, such as the ``resolution``, the
number of echelle ``orders`` (i.e., *red_ccd* or *blue_ccd*) and the wavelength
range for the rescaling (``wavelength_rescaling``), are listed. The path to the
archive where the spectra to be analyzed are located (``data_archive``) is also
listed here.

### spectral_lines 
In this section the spectral lines to analyze are listed (e.g., ``Na``,
``Halpha``).

For each line the user should indicate: the spectral ``range`` over which to
calculate the transmission spectrum (it must be wide enough to contain both the
stellar lines and the continuum); the wavelength [Å] of the spectral ``lines``
to analyze; the width of the central ``passbands`` around each spectral line for
the calculation of the relative absorption depth (AD) in Å (e.g., ``c075``:
*0.750*); the left/blue (``B``) and right/red (``R``) reference bands for the
calculation of the relative AD; the polynomial degree for the ``normalization``
of spectra and models.

The ``fit_parameters`` and the ``sampler_parameters`` for the Markov chain Monte Carlo (MCMC) analysis to fit the transmission spectrum are also listed in this section. Here, the user should set the number of steps and walkers (``n_step`` and ``n_walkers`` respectively), the range where to compute the fit, the binning step (``bin_step``) and the parameters ``priors``. Besides, there is a flag to decide whether to leave free the effective planet radius scale factor (``free_Rp``) and the radial velocity shift due to the wind (``free_winds``). If the MCMC analysis is to be performed on several spectral lines at the same time (e.g., the two lines of the sodium doublet or the magnesium triplet), the user can also decide whether the lines should share the same FWHM (``shared_fwhm``: *True*) or the same radial velocity shift due to the wind (``shared_winds``: *True*).

### planets 
In this section the planetary parameters are listed: the kind of ``orbit``
(‘circular’ if the user assumes e=0); the orbital ``eccentricity``; the orbital
``period`` [days]; the ``RV_semiamplitude`` of the planet [km/s]; the
``total_transit_duration`` (between the first and the fourth contact); the
``full_transit_duration`` (between the second and the third contact); the orbital
``inclination`` [deg]; the ``reference_time_of_transit`` [BJD]; the ``radius_ratio`` (ratio
of planet radius to stellar radius); the ``semimajor_axis_ratio`` (ratio of
semi-major axis to stellar radius); the ``impact_parameter``.

### star
In this section the stellar parameters are listed: the stellar ``mass`` [solar masses]); the stellar ``radius`` [solar radii]; the differential rotation rate (``alpha``); the projected obliquity [deg] (``lambda``); the projected rotational velocity [km/s] (``vsini``); the stellar ``inclination`` [deg]; the RV_semiamplitude of the star [km/s]; the systemic velocity  [km/s] (``RV_gamma``).
