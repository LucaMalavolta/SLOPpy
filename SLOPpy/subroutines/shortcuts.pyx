from __future__ import print_function, division
from SLOPpy.subroutines.io_subroutines import *
from SLOPpy.subroutines.common import *

def retrieve_observations(output_name, night, observations,
                          use_refraction=True, use_telluric=True, use_interstellar=True,
                          use_telluric_spline= False):

    #   config_in['output'], night, lists['observations']

    """ Retrieving the observations"""
    try:
        sky_corrected = load_from_cpickle('skycorrected_fibA', output_name, night)
        input_data = load_from_cpickle('input_dataset_fibA', output_name, night)
        for obs in observations:
            input_data[obs]['e2ds'] = sky_corrected[obs]['e2ds'].copy()
        print("  Sky contamination correction: applied ")
    except:
        input_data = load_from_cpickle('input_dataset_fibA', output_name, night)

    try:
        """ Retrieving the telluric correction"""
        telluric = load_from_cpickle('telluric', output_name, night)
        correct_telluric = True
    except:
        correct_telluric = False
        if use_telluric:
            print("  No telluric correction available - is this expected? ")

    """ Retrieval of refraction correction, if present """
    try:
        try:
            refraction = load_from_cpickle('refraction_update', output_name, night)
            correct_refraction = True
        except:
            refraction = load_from_cpickle('refraction', output_name, night)
            correct_refraction = True
    except:
        correct_refraction = False
        if use_refraction:
            print("  No refraction correction available - is this expected? ")

    """ Retrieval of refraction correction, if present """
    try:
        interstellar = load_from_cpickle('interstellar_lines', output_name, night)
        correct_interstellar = True
    except:
        correct_interstellar = False
        if use_interstellar:
            print("  No interstellar lines correction available - is this expected?  ")

    for obs in observations:

        """ Spectra are collected for variations of differential refraction,
        if the correction has been computed  before"""
        if correct_refraction and use_refraction:
            try:
                input_data[obs]['e2ds'] /= refraction[obs]['fit_e2ds']
                input_data[obs]['e2ds_err'] /= refraction[obs]['fit_e2ds']
            except:
                input_data[obs]['e2ds'] /= refraction[obs]['polyfit_e2ds']
                input_data[obs]['e2ds_err'] /= refraction[obs]['polyfit_e2ds']

        if correct_telluric and use_telluric:
            if use_telluric_spline:
                input_data[obs]['e2ds'] /= telluric[obs]['spline']
                input_data[obs]['e2ds_err'] /= telluric[obs]['spline']
            else:
                input_data[obs]['e2ds'] /= telluric[obs]['spectrum']
                input_data[obs]['e2ds_err'] /= telluric[obs]['spectrum']

        if correct_interstellar and use_interstellar:
            input_data[obs]['e2ds'] /= interstellar[obs]['correction']
            input_data[obs]['e2ds_err'] /= interstellar[obs]['correction']


    if correct_refraction and use_refraction:
        print("  Differential refraction correction: applied ")
    if correct_telluric and use_telluric:
        print("  Telluric correction: applied ")
    if correct_interstellar and use_interstellar:
        print("  Interstellar lines correction: applied ")

    print()
    return input_data


def compute_rescaling(wave, flux, wavelength_range):

    wl_med = np.average(wavelength_range)
    wl_dif = wavelength_range[1]-wl_med
    sel = np.where(np.abs(wave-wl_med) < wl_dif)
    return np.median(flux[sel])


def perform_rescaling(wave, flux, ferr, wavelength_range):

    factor = compute_rescaling(wave, flux, wavelength_range)
    # flux_out = flux / factor
    # ferr_out = ferr / factor

    return factor, flux / factor, ferr / factor


def replace_values(vals, threshold, replacement=None):
    if not replacement:
        replacement = np.median(vals)

    null = (vals <= threshold)
    vals[null] = replacement
    return vals, null


def replace_values_errors(vals, errs, threshold, replacement=None, error_replacement=None):
    if not replacement:
        replacement = np.median(vals)
    if not error_replacement:
        error_replacement = np.median(errs)*10.0

    null = (vals <= threshold)
    vals[null] = replacement
    errs[null] = error_replacement
    return vals, errs, null


def replace_values_errors_with_interpolation_1d(vals,
                                                errs=None,
                                                less_than=None,
                                                greater_than=None,
                                                force_positive=False,
                                                error_replacement=None,
                                                sigma_iter=1):
    """
    Replace outliers in the input array and the associated arrors array (if provided)
    Written as a fix to out-of-boundaries rebinning
    At least one between  less_than  or  greater_than  must be provided

    :param vals: input array to be fixed
    :param errs: error array associated to input array
    :param less_than: lower threshold for considering a value as outlier.
    :param greater_than: upper threshold for considering a value as outlier.
    :param error_replacement: optional, error to be associated to the outlier.
        if not provided, ten times the error of the interpolated points is taken
    :param sigma_iter: optional, if provided it will consider less_than/greater_than as sigma multiplicators
    :return: vals, errs, null
    """

    if errs is None:
        errt = vals*0.0
    else:
        errt = errs

    if sigma_iter <= 1:
        use_sigma = False
    else:
        use_sigma = True

    if force_positive:
        null = (vals <= 0.00000000001)
    else:
        null = (~np.isfinite(vals))

    for ii in range(0, sigma_iter):

        if use_sigma:
            median = np.median(vals)
            sigma = np.std(vals)
            if less_than:
                threshold_less = median - sigma*less_than
            if greater_than:
                threshold_more = median + sigma*greater_than
        else:
            if less_than:
                threshold_less = less_than
            if greater_than:
                threshold_more = greater_than

        if less_than and greater_than:
            ind_sel = np.where((vals < threshold_less) | (vals > threshold_more) | null)[0]
        elif less_than:
            ind_sel = np.where((vals < threshold_less) | null)[0]
        elif greater_than:
            ind_sel = np.where((vals > threshold_more) | null)[0]
        elif force_positive:
            ind_sel = np.where((vals < 0.00001) | null)[0]
        else:
            raise ValueError('Provide at least on condition')


        null[ind_sel] = True
        i_e = 0
        for jj in ind_sel:
            if jj < i_e: continue
            i_s = jj
            i_e = jj

            while i_s in ind_sel:
                i_s -= 1
                if i_s < 0:
                    break

            while i_e in ind_sel:
                i_e += 1
                if i_e >= vals.size:
                    break

            if i_s < 0 and i_e >= vals.size:
                raise ValueError('Invalid input, something must be wrong in the input array')
            elif i_s < 0:
                vals[:i_e] = vals[i_e]
                errt[:i_e] = errt[i_e]*10.
            elif i_e >= vals.size:
                vals[i_s:] = vals[i_s]
                errt[i_s:] = errt[i_s]*10.
            else:
                vals[i_s+1:i_e] = (vals[i_e]+vals[i_s])/2.
                errt[i_s+1:i_e] = (vals[i_e]+vals[i_s])*10.

    if errs is None:
        return vals, null
    else:
        if error_replacement:
            errs[ind_sel] = error_replacement
        else:
            errs[ind_sel] = errt[ind_sel]
        return vals, errs, null


def replace_values_errors_with_interpolation_2d(vals,
                                                errs=None,
                                                less_than=None,
                                                greater_than=None,
                                                force_positive=False,
                                                error_replacement=None,
                                                sigma_iter=1,
                                                axis=0):

    if errs is None:
        errt = vals*0.0
    else:
        errt = errs

    null = (~np.isfinite(vals))

    if axis ==0:
        for n_order in range(0, np.size(vals[:,0])):
            vals[n_order, :], errt[n_order,:], null[n_order, :] = replace_values_errors_with_interpolation_1d(
                vals[n_order, :],
                errt[n_order, :],
                less_than,
                greater_than,
                force_positive,
                error_replacement,
                sigma_iter)
    else:
        for n_order in range(0, np.size(vals[0,:])):
            vals[:, n_order], errt[:, n_order], null[:, n_order] = replace_values_errors_with_interpolation_1d(
                vals[:, n_order],
                errt[:, n_order],
                less_than,
                greater_than,
                force_positive,
                error_replacement,
                sigma_iter)

    if errs is None:
        return vals, null
    else:
        return vals, errs, null


def compute_spline(wave_in, flux_in, knots_step, knot_order=3, use_null=None, use_selection=None):
    """ Sometimes the interpolation range is just outside the spline evaluation range. The spline procedure will
    fail even if the affected point are at the very extremes of the spectral range, i.e., points that would not be
    inin the analysis anyway. We trick the spline procedure by adding fake points
    """

    if use_null is not None:
        selection = (~use_null)
    elif use_selection is not None:
        selection = use_selection
    else:
        selection = (wave_in > 0.000)

    len_wave = np.sum(selection) + 6
    wave = np.empty(len_wave, dtype=np.double)
    spec = np.empty(len_wave, dtype=np.double)

    wave[:3] = wave_in[0] - np.linspace(1.0, 0.5, 3)
    wave[3:-3] = wave_in[selection]
    wave[-3:] = wave_in[-1] + np.linspace(0.5, 1.0, 3)

    spec[:3] = flux_in[selection][0]
    spec[3:-3] = flux_in[selection]
    spec[-3:] = flux_in[selection][-1]

    """ computing the spline approximation of the rescaled master out"""

    """ picking the number of knots """
    nknots = (np.amax(wave) - np.amin(wave)) /knots_step
    """ picking the indices of the knots"""
    idx_knots = (np.arange(1, len(wave) - 1, (len(wave) - 2.) / nknots)).astype('int')

    """ passing from indices to knots values """
    spline_knots = wave[idx_knots]

    spline_coeff = sci_int.splrep(wave, spec, task=-1, k=knot_order, t=spline_knots)
    spline_eval = sci_int.splev(wave_in, spline_coeff)

    return spline_eval, spline_coeff, spline_knots


def f_telluric_rescaling_factor(x, x_range, tel1, tel2):
    return np.sum((tel1-(1.+x[0]*x_range)-(tel2-1.)*x[1])**2)


def find_telluric_rescaling_factor(tel1, tel2):
    fit_selection = (tel2 < 1.0)

    x_range = np.arange(0.0, len(tel1), 1.0)/len(tel1)
    x_start = np.zeros(2)
    x_start[1] = np.median((tel1[fit_selection]-1.)/(tel2[fit_selection]-1.))
    x_result = sci_optimize.minimize(f_telluric_rescaling_factor,
                                     x_start,
                                     args=(x_range, tel1, tel2),
                                     method='Nelder-Mead',
                                     options=dict(maxiter=10000))
    return x_result['x'][1], x_result['x'][0], x_result['success']



def f_telluric_xslope(x, x_range, tel1, tel2, std):
    return np.sum((tel1-(1.+x*x_range))**2/std**2)


def f_telluric_factor(x, tel1, tel2, std):
    print(x, np.sum(((tel1-1.)-(tel2-1.)*x)**2/std**2))
    return np.sum(((tel1-1.)-(tel2-1.)*x)**2/std**2)


def find_telluric_rescaling_factor_2steps(tel1, tel2):
    fit_xslope = (tel2 > 0.99999)
    fit_factor = (tel2 < 0.99999)
    std_all = np.std(tel1)

    x_range = np.arange(0.0, len(tel1), 1.0)/len(tel1)
    x_start = 0.000

    x_result = sci_optimize.minimize(f_telluric_xslope,
                                     x_start,
                                     args=(x_range[fit_xslope], tel1[fit_xslope], tel2[fit_xslope], std_all),
                                     method='Nelder-Mead',
                                     options=dict(maxiter=10000))

    x_xslope = x_result['x'][0]
    x_start = np.median((tel1[fit_factor]-1.-x_range[fit_factor]*x_xslope)/(tel2[fit_factor]-1.))
    print('---> ', x_start)
    x_result = sci_optimize.minimize(f_telluric_factor,
                                     x_start,
                                     args=(tel1[fit_factor]-x_range[fit_factor]*x_xslope, tel2[fit_factor], std_all),
                                     method='Nelder-Mead',
                                     options=dict(maxiter=10000))

    x_factor = x_result['x'][0]

    return x_factor, x_xslope, x_result['success']




def write_molecfit_v1_par(filename_par, filename_data,  filename_output, filename_include, molecfit_dict, observing_dict):


    """ The observing_dict must contains these keywords:
    'MJD': Observing date in years or MJD in days
    'UTC': UTC in s
    'ELEVATION': Telescope altitude angle in deg
    'HUMIDITY': Humidity in %
    'PRESSURE': Pressure in hPa
    'TEMPERATURE_EN': Ambient temperature in deg C
    'TEMPERATURE_M1': Mirror temperature in deg C
    'GEOELEV': Elevation above sea level in m (default is Paranal: 2635m)
    'GEOLONG': Longitude
    'GEOLAT': Latitude
    """

    fileout = open(filename_par, 'w')
    fileout.write("### Driver for MOLECFIT\n")

    # user working directory only important for REFLEX workflow and GUI
    # not used by molecfit itself.
    fileout.write("user_workdir:./\n")

    ## INPUT DATA
    # Data file name (path relative to the current directory or absolute path)
    fileout.write("filename: " + filename_data + "\n")

    # ASCII list of files to be corrected for telluric absorption using the
    # transmission curve derived from the input reference file (path of list and
    # listed files relative to the current directory or absolute path; default: "none")
    fileout.write("listname: none\n")

    # Type of input spectrum -- 1 = transmission (default); 0 = emission
    fileout.write("trans: 1\n")

    # Names of the file columns (table) or extensions (image) containing:
    # Wavelength  Flux  Flux_Err  Mask
    # - Flux_Err and/or Mask can be avoided by writing 'NULL'
    # - 'NULL' is required for Wavelength if it is given by header keywords
    # - parameter list: col_lam, col_flux, col_dflux, and col_mask
    fileout.write("columns: Wavelength Flux NULL NULL\n")

    # Default error relative to mean for the case that the error column is missing
    fileout.write("default_error: 0.001\n")

    # Multiplicative factor to convert wavelength to micron
    # (e.g. nm -> wlgtomicron = 1e-3)
    fileout.write("wlgtomicron: 0.0001\n")

    # Wavelengths in vacuum (= vac) or air (= air)
    fileout.write("vac_air: air\n")

    # TODO: input from configuration file for molecfit installation path
    # ASCII or FITS table for wavelength ranges in micron to be fitted
    # (path relative to the current directory or absolute path; default: "none")
    fileout.write("wrange_include: " + filename_include + "\n")

    # ASCII or FITS table for wavelength ranges in micron to be excluded from the
    # fit (path relative to the current directory or absolute path; default: "none")
    # wrange_exclude: /Users/malavolta/Astro/ExoAtmospheres/molecfit_test//HIP63901_exclude_w.dat

    # ASCII or FITS table for pixel ranges to be excluded from the fit
    # (path relative to the current directory or absolute path; default: "none")
    # prange_exclude: /Users/malavolta/Astro/ExoAtmospheres/molecfit_test//HIP63901_exclude_p.dat

    ## RESULTS

    # Directory for output files (path relative to the current directory or absolute path)
    fileout.write("output_dir:./\n")

    # Name for output files
    # (supplemented by "_fit" or "_tac" as well as ".asc", ".atm", ".fits",
    # ".par, ".ps", and ".res")
    fileout.write("output_name: "+ filename_output + "\n")

    # Plot creation: gnuplot is used to create control plots
    # W - screen output only (incorporating wxt terminal in gnuplot)
    # X - screen output only (incorporating x11 terminal in gnuplot)
    # P - postscript file labelled '<output_name>.ps', stored in <output_dir>
    # combinations possible, i.e. WP, WX, XP, WXP (however, keep the order!)
    # all other input: no plot creation is performed
    fileout.write("plot_creation: none\n")

    # Create plots for individual fit ranges? -- 1 = yes; 0 = no
    fileout.write("plot_range: 0\n")

    ## FIT PRECISION

    # Relative chi2 convergence criterion
    fileout.write("ftol: " + molecfit_dict['ftol'] + "\n")

    # Relative parameter convergence criterion
    fileout.write("xtol: " + molecfit_dict['xtol'] + "\n")

    ## MOLECULAR COLUMNS

    # List of molecules to be included in the model
    # (default: 'H2O', N_val: nmolec)
    molecules_list = "list_molec:"
    for mol in molecfit_dict['molecules']:
        molecules_list += " " + mol
    fileout.write(molecules_list + "\n")

    # Fit flags for molecules -- 1 = yes; 0 = no (N_val: nmolec)
    fileout.write("fit_molec: 1 1\n")

    # Values of molecular columns, expressed relatively to the input ATM profile
    # columns (N_val: nmolec) [1 = 100%]
    fileout.write("relcol: 1.0 1.0\n")

    ## BACKGROUND AND CONTINUUM

    # Conversion of fluxes from phot/(s*m2*mum*as2) (emission spectrum only) to
    # flux unit of observed spectrum:
    # 0: phot/(s*m^2*mum*as^2) [no conversion]
    # 1: W/(m^2*mum*as^2)
    # 2: erg/(s*cm^2*A*as^2)
    # 3: mJy/as^2
    # For other units the conversion factor has to be considered as constant term
    # of the continuum fit.
    fileout.write("flux_unit: 0\n")

    # Fit of telescope background -- 1 = yes; 0 = no (emission spectrum only)
    fileout.write("fit_back: 0\n")

    # Initial value for telescope background fit (range: [0,1])
    fileout.write("telback: 0.1\n")

    # Polynomial fit of continuum --> degree: cont_n
    fileout.write("fit_cont: 1\n")

    # Degree of coefficients for continuum fit
    fileout.write("cont_n: {0:1.0f}".format(molecfit_dict['cont_n']) + "\n")

    # Initial constant term for continuum fit (valid for all fit ranges)
    # (emission spectrum: about 1 for correct flux_unit)
    fileout.write("cont_const: {0:1.0f}".format(molecfit_dict['cont_const']) + "\n")

    ## WAVELENGTH SOLUTION

    # Refinement of wavelength solution using a polynomial of degree wlc_n
    fileout.write("fit_wlc: 1\n")

    # Polynomial degree of the refined wavelength solution
    fileout.write("wlc_n: {0:1.0f}".format(molecfit_dict['wlc_n']) + "\n")

    # Initial constant term for wavelength correction (shift relative to half
    # wavelength range)
    fileout.write("wlc_const: {0:1.0f}".format(molecfit_dict['wlc_const']) + "\n")

    ## RESOLUTION

    # Fit resolution by boxcar -- 1 = yes; 0 = no
    fileout.write("fit_res_box: 0\n")

    # Initial value for FWHM of boxcar relative to slit width (>= 0. and <= 2.)
    fileout.write("relres_box: 0.0\n")

    # Voigt profile approximation instead of independent Gaussian and Lorentzian
    # kernels? -- 1 = yes; 0 = no
    fileout.write("kernmode: 0\n")

    # Fit resolution by Gaussian -- 1 = yes; 0 = no
    fileout.write("fit_res_gauss: 1\n")

    # Initial value for FWHM of Gaussian in pixels
    fileout.write("res_gauss: {0:3.1f}".format(molecfit_dict['res_gauss']) + "\n")

    # Fit resolution by Lorentzian -- 1 = yes; 0 = no
    fileout.write("fit_res_lorentz: 0\n")

    # Initial value for FWHM of Lorentzian in pixels
    fileout.write("res_lorentz: 0.0\n")

    # Size of Gaussian/Lorentzian/Voigtian kernel in FWHM
    fileout.write("kernfac: {0:3.0f}".format(molecfit_dict['kernfac']) + "\n")

    # Variable kernel (linear increase with wavelength)? -- 1 = yes; 0 = no
    fileout.write("varkern: 0\n")

    # ASCII file for kernel elements (one per line; normalisation not required)
    # instead of synthetic kernel consisting of boxcar, Gaussian, and Lorentzian
    # components (path relative to the current directory or absolute path; default: "none\n")
    fileout.write("kernel_file: none\n")

    ## AMBIENT PARAMETERS

    # If the input data file contains a suitable FITS header, the keyword names of
    # the following parameters will be read, but the corresponding values will not
    # be used. The reading of parameter values from this file can be forced by
    # setting keywords to NONE.

    # Observing date in years or MJD in days
    # changed the format to inter in order to not have molecfit crashed
    #fileout.write("obsdate: {0:13.5f}".format(observing_dict['MJD']) + "\n")
    fileout.write("obsdate: {0:13.50}".format(observing_dict['MJD']) + "\n")
    fileout.write("obsdate_key: NONE\n")

    # UTC in s
    fileout.write("utc: {0:8.0f}".format(observing_dict['UTC']) + "\n")
    fileout.write("utc_key: NONE\n")

    # Telescope altitude angle in deg
    fileout.write("telalt: {0:13.5f}".format(observing_dict['ELEVATION']) + "\n")
    fileout.write("telalt_key: NONE\n")

    # Humidity in %
    fileout.write("rhum: {0:13.5f}".format(observing_dict['HUMIDITY']) + "\n")
    fileout.write("rhum_key: NONE\n")

    # Pressure in hPa
    fileout.write("pres: {0:5.1f}".format(observing_dict['PRESSURE']) + "\n")
    fileout.write("pres_key: NONE\n")

    # Ambient temperature in deg C
    fileout.write("temp: {0:4.1f}".format(observing_dict['TEMPERATURE_EN']) + "\n")
    fileout.write("temp_key: NONE\n")

    # Mirror temperature in deg C
    fileout.write("m1temp: {0:4.1f}".format(observing_dict['TEMPERATURE_M1']) + "\n")
    fileout.write("m1temp_key: NONE\n")

    # Elevation above sea level in m (default is Paranal: 2635m)
    fileout.write("geoelev: {0:4.0f}".format(observing_dict['GEOELEV']) + "\n")
    fileout.write("geoelev_key: NONE\n")

    # Longitude (default is Paranal: -70.4051)
    fileout.write("longitude: {0:9.4f}".format(observing_dict['GEOLONG']) + "\n")
    fileout.write("longitude_key: NONE\n")

    # Latitude (default is Paranal: -24.6276)
    fileout.write("latitude: {0:9.4f}".format(observing_dict['GEOLAT']) + "\n")
    fileout.write("latitude_key: NONE\n")

    ## INSTRUMENTAL PARAMETERS

    # Slit width in arcsec (taken from FITS header if present)
    fileout.write("slitw: {0:3.1f}".format(molecfit_dict['slitwidth']) + "\n")
    fileout.write("slitw_key: NONE\n")

    # Pixel scale in arcsec (taken from this file only)
    fileout.write("pixsc: {0:4.2f}".format(molecfit_dict["pixelscale"]) + "\n")
    fileout.write("pixsc_key: NONE\n")

    ## ATMOSPHERIC PROFILES

    # Reference atmospheric profile
    fileout.write("ref_atm: equ.atm\n")

    # Specific GDAS-like input profile (P[hPa] HGT[m] T[K] RELHUM[%]) (path
    # relative to the installation directory or absolute path). In the case of "none", no GDAS
    # profiles will be considered. The default "auto" performs an automatic
    # retrieval.
    fileout.write("gdas_dir: data/profiles/grib\n")
    fileout.write("gdas_prof: auto\n")

    # Grid of layer heights for merging ref_atm and GDAS profile. Fixed grid = 1
    # (default) and natural grid = 0.
    fileout.write("layers: 0\n")

    # Upper mixing height in km (default: 5) for considering data of a local meteo
    # station. If emix is below geoelev, rhum, pres, and temp are not used for
    # modifying the corresponding profiles.
    fileout.write("emix: 5.0\n")

    # PWV value in mm for the input water vapour profile. The merged profile
    # composed of ref_atm, GDAS, and local meteo data will be scaled to this value
    # if pwv > 0 (default: -1 -> no scaling).
    fileout.write("pwv: -1.\n")

    # internal GUI specific parameter
    fileout.write("clean_mflux: 1\n")

    fileout.write("end\n")
    fileout.close()



def write_molecfit_par(filename_par, wave_include, molecfit_dict, observing_dict):


    """ The observing_dict must contains these keywords:
    'MJD': Observing date in years or MJD in days
    'UTC': UTC in s
    'ELEVATION': Telescope altitude angle in deg
    'HUMIDITY': Humidity in %
    'PRESSURE': Pressure in hPa
    'TEMPERATURE_EN': Ambient temperature in deg C
    'TEMPERATURE_M1': Mirror temperature in deg C
    'GEOELEV': Elevation above sea level in m (default is Paranal: 2635m)
    'GEOLONG': Longitude
    'GEOLAT': Latitude
    """

    fileout = open(filename_par, 'w')

    # File: molecfitConfigFiles/model.rc
    #
    # Note: This configuration file has been automatically
    #       generated by the esorex (v3.13.5) program.
    #
    # Date: 25-Jun-2022 16:52:38
    #
    #

    # --USE_ONLY_INPUT_PRIMARY_DATA
    # Value=TRUE implies that only the fits primary contains the input science flux
    # data.
    # Value=FALSE implies that the fits extensions also contains input science
    # flux data.
    fileout.write('USE_ONLY_INPUT_PRIMARY_DATA=FALSE\n')

    # --USE_DATA_EXTENSION_AS_DFLUX
    # Only valid if USE_ONLY_INPUT_PRIMARY_DATA=TRUE. The fits extension index that
    # contains the
    # errors of the science flux data (DFLUX). A value of 0 implies that there is
    # no DFLUX.
    fileout.write('USE_DATA_EXTENSION_AS_DFLUX=0\n')

    # --USE_DATA_EXTENSION_AS_MASK
    # Only valid if USE_ONLY_INPUT_PRIMARY_DATA=TRUE. The fits extension index that
    # contains the
    # mask associated with the science flux data. A value of 0 implies that there
    # is no mask data.
    fileout.write('USE_DATA_EXTENSION_AS_MASK=0\n')

    # --USE_INPUT_KERNEL
    # If TRUE, use the kernel library if it is provided.
    fileout.write('USE_INPUT_KERNEL=TRUE\n')

    # --MODEL_MAPPING_KERNEL
    # Mapping 'STD_MODEL/SCIENCE' - 'MODEL_KERNEL_LIBRARY' [string with ext_number
    # comma separated (int)] :
    # If set to NULL, check if the TAG[MODEL_MAPPING_KERNEL] FITS BINTABLE values
    # is provided.
    # The FITS BINTABLE have to one column [KERNEL_LIBRARY_EXT].
    fileout.write('MODEL_MAPPING_KERNEL=NULL\n')


    ## MOLECULAR COLUMNS

    # --LIST_MOLEC
    # List of molecules to be included in the model. Represented as a comma
    # separated
    # string of molecule names, e.g. "H2O,CO2,O3".
    # If set to NULL, the input TAG[MOLECULES] FITS BINTABLE values have to be
    # provided
    # where the FITS BINTABLE specified contains the three columns:
    # LIST_MOLEC; FIT_MOLEC; and REL_COL.
    molecules_list = 'LIST_MOLEC="'

    # --FIT_MOLEC
    # List of flags that specify which of the listed molecules are to be fitted for.
    # Flag=1 implies yes. Flag=0 implies no. Represented as a string of comma
    # separated
    # integers in the same order as the listed molecules. For example: if
    # LIST_MOLEC="H2O,CO2,O3", then
    # FIT_MOLEC="1,0,1" implies that only H2O and O3 should be fitted for.
    # If set to NULL, the input TAG[MOLECULES] FITS BINTABLE values have to be
    # provided where the FITS
    # BINTABLE specified contains the three columns: LIST_MOLEC; FIT_MOLEC; and
    # REL_COL.
    molecules_flag = 'FIT_MOLEC="'

    # --REL_COL
    # List of the intial values of fitting of the molecular columns expressed
    # relatively to the input
    #  ATM profile columns. Represented as a comma separated list of doubles in
    # the same order as the
    # listed molecules. For example, if LIST_MOLEC="H2O,CO2,O3", then
    # REL_COL="1.0,1.2,0.8"
    # implies that H2O, CO2 and O3 have initial relative values of 1.0, 1.2 and
    # 0.8 respectively.
    # If set to NULL, the input TAG[MOLECULES] FITS BINTABLE values have to be
    # provided where the FITS
    # BINTABLE specified contains the three columns: LIST_MOLEC; FIT_MOLEC; and
    # REL_COL.
    molecules_rel = 'REL_COL="'
    for mol in molecfit_dict['molecules']:
        molecules_list += mol
        molecules_list += ','
        molecules_flag += '1,'
        molecules_rel += '1.0,'

    molecules_list = molecules_list[:-1] + '"'
    molecules_flag = molecules_flag[:-1] + '"'
    molecules_rel = molecules_rel[:-1] + '"'

    fileout.write(molecules_list + '\n')
    fileout.write(molecules_flag + '\n')
    fileout.write(molecules_rel + '\n')

    # --WAVE_INCLUDE
    # Wavelength ranges to be included. Represented as a string of comma separated
    # doubles in pairs
    # specifying the start and end wavelengths of a range. The wavelength units
    # are always in microns.
    # For example a KMOS sample data in the range of 1.11um to 1.67um may have
    # WAVE_INCLUDE="1.773,1.78633,1.79098,1.80434,1.187691,1.189937" to represent
    # three inclusion regions:
    # [1.773,1.78633], [1.79098,1.80434] and [1.187691,1.189937].
    # If set to NULL, molecfit will check if the TAG[WAVE_INCLUDE] FITS BINTABLE
    # values is provided where
    # the FITS BINTABLE specified has the two columns: LOWER_LIMIT; and
    # UPPER_LIMIT.
    fileout.write('WAVE_INCLUDE='+wave_include+'\n')

    # --WAVE_EXCLUDE
    # Wavelength ranges to be excluded. Represented as a string of comma separated
    # doubles in pairs
    # specifying the start and end wavelengths of a range. The wavelength units
    # are always in microns.
    # as the input science data. For example a KMOS sample data in the range of
    # 1.11um to 1.67um may have
    # WAVE_EXCLUDE="1.773,1.78633,1.79098,1.80434,1.187691,1.189937" to represent
    # three exclusion regions:
    # [1.773,1.78633], [1.79098,1.80434] and [1.187691,1.189937].
    # If set to NULL, molecfit will check if the TAG[WAVE_EXCLUDE] FITS BINTABLE
    # values is provided where
    # the FITS BINTABLE specified has the two columns: LOWER_LIMIT; and
    # UPPER_LIMIT.
    fileout.write('WAVE_EXCLUDE=NULL\n')


    # --PIXEL_EXCLUDE
    # Pixel ranges to be excluded. Represented as a string of comma separated
    # integers in pairs specifying the
    # start and end pixel of a range. For example:
    # PIXEL_EXCLUDE="54,128,512,514,1020,1024" represents three
    # exclusion regions: [54,128], [512,514] and [1020,1024].
    # If set to NULL, molecfit will check if the TAG[PIXEL_EXCLUDE] FITS BINTABLE
    # values is provided where the
    # FITS BINTABLE specified has the two columns: LOWER_LIMIT; and UPPER_LIMIT.
    fileout.write('PIXEL_EXCLUDE=NULL\n')

    # --TELLURICCORR_PATH
    # Installation directory.
    fileout.write('TELLURICCORR_PATH=TELLURICCORR_PARAMETER_DEFAULT\n')

    # --TELLURICCORR_DATA_PATH
    # Data directory.
    fileout.write('TELLURICCORR_DATA_PATH=TELLURICCORR_PARAMETER_DEFAULT\n')

    # --TMP_PATH
    # Temporary directory.
    fileout.write('TMP_PATH=TELLURICCORR_PARAMETER_DEFAULT\n')

    # --SILENT_EXTERNAL_BINS
    # Silent the output of the external binaries.
    fileout.write('SILENT_EXTERNAL_BINS=TRUE\n')

    # --TRANSMISSION
    # Type of input spectrum : 0 = Emission(radiance); 1 = Transmission.
    fileout.write('TRANSMISSION=TRUE\n')

    # --COLUMN_LAMBDA
    # Wavelength column ('NULL' can be used if the file is an image and that
    # the data are in the primary
    # (data are given by the FITS header keywords [CRVAL1=wave_ini, CD1_1=step])
    # f CD1_1 is absent, then the DEPRECATED CDELT1 keyword will be used.
    fileout.write('COLUMN_LAMBDA=lambda\n')

    # --COLUMN_FLUX
    # Flux column.
    fileout.write('COLUMN_FLUX=flux\n')

    # --COLUMN_DFLUX
    # Flux error column (Avoided by writing 'NULL') : 1-sigma error on the flux.
    fileout.write('COLUMN_DFLUX=NULL\n')

    # --COLUMN_MASK
    # Mask column (Avoided by writing 'NULL') : Indicates if a pixel is invalid.
    fileout.write('COLUMN_MASK=NULL\n')

    # --DEFAULT_ERROR
    # Default error relative to mean for the case that the error column
    # is not provided.
    fileout.write('DEFAULT_ERROR=0.01\n')

    # --WLG_TO_MICRON
    # Multiplicative factor applied to the wavelength to express is in micron.
    # E.g.: if wavelength is given in nm, the value should be 0.001.
    fileout.write('WLG_TO_MICRON=1.\n')

    # --WAVELENGTH_FRAME
    # Wavelength in vacuum                                      = 'VAC'.
    # Wavelength in air    with the observatory reference frame = 'AIR'.
    # Wavelength in vacuum with another         reference frame = 'VAC_RV'.
    #   (typically the  sun or the barycenter of the solar system).
    # In the latter case, the radial velocity of the observatory relative
    #   to the external reference frame must be provided in the parameter obs_RV.
    fileout.write('WAVELENGTH_FRAME=AIR\n')
    # TODO should I convert everything in vacuum or air??


    # --OBS_ERF_RV_KEY
    # The radial velocity of the observatory in km/s
    # relative to the external reference frame;
    # It is positive if the distance between the science target and the Earth
    # increases along the line-of-sight to the science target.
    # It must be provided if MF_PARAMETERS_WAVELENGTH_FRAME = 'VAC_RV'.
    fileout.write('OBS_ERF_RV_KEY=NONE\n')

    # --OBS_ERF_RV_VALUE
    # If OBS_ERF_RV_KEYWORD=='NONE' take this value.
    fileout.write('OBS_ERF_RV_VALUE=0.0\n')

    # --CLEAN_MODEL_FLUX
    # Set model flux to 0 for non-fitted pixels.
    fileout.write('CLEAN_MODEL_FLUX=FALSE\n')

    # --FTOL
    # Relative chi-square convergence criterion.
    # FTOL=1e-10
    try:
        fileout.write('FTOL={0:f}\n'.format(molecfit_dict['ftol']))
    except ValueError:
        fileout.write("FTOL=" + molecfit_dict['ftol'] + "\n")

    # --XTOL
    # Relative parameter convergence criterion.
    # XTOL=1e-10
    try:
        fileout.write('XTOL={0:f}\n'.format(molecfit_dict['xtol']))
    # Relative chi2 convergence criterion
    except ValueError:
        fileout.write("XTOL=" + molecfit_dict['xtol'] + "\n")



    # --FLUX_UNIT
    # Conversion of fluxes from phot/(s*m2*mum*as2) (emission spectrum only)
    #    to flux unit of observed spectrum:
    # 0: phot / (s *  m^2 * mum * as^2) [no conversion]
    # 1:    W / (     m^2 * mum * as^2)
    # 2:  erg / (s * cm^2 *   A * as^2)
    # 3:  mJy / (                 as^2)
    # For other units, the conversion factor has to be considered
    #  as constant term of the continuum fit.
    fileout.write('FLUX_UNIT=0\n')

    # --FIT_TELESCOPE_BACKGROUND
    # Fit of telescope background -- 1 = yes; 0 = no (emission spectrum only).
    fileout.write('FIT_TELESCOPE_BACKGROUND=TRUE\n')
    #todo check this, it was zero in previous file format

    # --TELESCOPE_BACKGROUND_CONST
    # Initial value for telescope background fit.
    fileout.write('TELESCOPE_BACKGROUND_CONST=0.1\n')

    # --FIT_CONTINUUM
    # Comma deliminated string of flags (1=true, 0=false) for fitting continuum in
    # specific regions.
    # If set to NULL, check if the TAG[WAVE_INCLUDE] points to a FITS BINTABLE
    # with column CONT_FIT_FLAG provided.
    fileout.write('FIT_CONTINUUM=1\n')

    # --CONTINUUM_N
    # Polynomial order for continuum model for each region. Presented as a comma
    # deliminated string.
    # If set to NULL, check if the TAG[WAVE_INCLUDE] points to a FITS BINTABLE
    # with column CONT_POLY_ORDER provided.
    fileout.write('CONTINUUM_N={0:1.0f}\n'.format(molecfit_dict['cont_n']))

    # --CONTINUUM_CONST
    # Initial constant term for continuum fit (valid for all fit ranges)
    # [emission spectrum: about 1 for correct flux_unit].
    fileout.write('CONTINUUM_CONST={0:10f}\n'.format(molecfit_dict['cont_const']))

    # --FIT_WLC
    # Flags for including regions in wavelength corrections.
    fileout.write('FIT_WLC=1\n')

    # --WLC_N
    # Polynomial degree of the refined wavelength solution.
    fileout.write('WLC_N={0:1.0f}\n'.format(molecfit_dict['wlc_n']))

    # --WLC_CONST
    # Initial constant term for wavelength adjustment
    # (shift relative to half wavelength range).
    fileout.write('WLC_CONST={0:10f}\n'.format(molecfit_dict['wlc_const']))

    # --FIT_RES_BOX
    # Fit resolution by Boxcar LSF.
    fileout.write('FIT_RES_BOX=FALSE\n')

    # --RES_BOX
    # Initial value for FWHM of Boxcar rel. to slit width
    # at the centre of the spectrum.
    fileout.write('RES_BOX=0.0\n')

    # --FIT_RES_GAUSS
    # Fit resolution by Gaussian.
    fileout.write('FIT_RES_GAUSS=TRUE\n')


    # --RES_GAUSS
    # Initial value for FWHM of the Gaussian in pixels
    # at the centre of the spectrum.
    fileout.write('RES_GAUSS={0:10f}\n'.format(molecfit_dict['res_gauss']))

    # --FIT_RES_LORENTZ
    # Fit resolution by Lorentzian.
    fileout.write('FIT_RES_LORENTZ=FALSE\n')

    # --RES_LORENTZ
    # Initial value for FWHM of the Lorentz in pixels
    # at the centre of the spectrum.
    fileout.write('RES_LORENTZ=0.0\n')

    # --KERNMODE
    # Voigtian profile approximation instead of independent Gaussian and
    # Lorentzian?.
    fileout.write('KERNMODE=FALSE\n')

    # --KERNFAC
    # Size of Voigtian/Gaussian/Lorentzian kernel in FWHM.
    fileout.write('KERNFAC={0:10f}\n'.format(molecfit_dict['kernfac']))

    # --VARKERN
    # Does the kernel size increase linearly with wavelength?.
    fileout.write('VARKERN=TRUE\n')
    #TODO: check why it was FALSE in the previous iteration

    # --OBSERVING_DATE_KEYWORD
    # Observing date in years or MJD in days (not string).
    fileout.write('OBSERVING_DATE_KEYWORD=NONE\n')

    # --OBSERVING_DATE_VALUE
    # If OBSERVING_DATE_KEYWORD=='NONE' take this value.
    fileout.write('OBSERVING_DATE_VALUE={0:30f}\n'.format(observing_dict['MJD']))

    # --UTC_KEYWORD
    # UTC in s.
    fileout.write('UTC_KEYWORD=NONE\n')

    # --UTC_VALUE
    # If UTC_KEYWORD=='NONE' take this value.
    fileout.write('UTC_VALUE={0:30f}\n'.format(observing_dict['UTC']))

    # --TELESCOPE_ANGLE_KEYWORD
    # Telescope altitude angle in deg.
    fileout.write('TELESCOPE_ANGLE_KEYWORD=NONE\n')

    # --TELESCOPE_ANGLE_VALUE
    # If TELESCOPE_ANGLE_KEYWORD=='NONE' take this value.
    fileout.write('TELESCOPE_ANGLE_VALUE={0:30f}\n'.format(observing_dict['ELEVATION']))

    # --RELATIVE_HUMIDITY_KEYWORD
    # Relative humidity in %.
    fileout.write('RELATIVE_HUMIDITY_KEYWORD=NONE\n')

    # --RELATIVE_HUMIDITY_VALUE
    # If RELATIVE_HUMIDITY_KEYWORD=='NONE' take this value.
    fileout.write('RELATIVE_HUMIDITY_VALUE={0:30f}\n'.format(observing_dict['HUMIDITY']))

    # --PRESSURE_KEYWORD
    # Pressure in hPa.
    fileout.write('PRESSURE_KEYWORD=NONE\n')

    # --PRESSURE_VALUE
    # If PRESSURE_KEYWORD=='NONE' take this value.
    fileout.write('PRESSURE_VALUE={0:5.1f}\n'.format(observing_dict['PRESSURE']))

    # --TEMPERATURE_KEYWORD
    # Ambient temperature in deg C.
    fileout.write('TEMPERATURE_KEYWORD=NONE\n')

    # --TEMPERATURE_VALUE
    # If TEMPERATURE_KEYWORD=='NONE' take this value.
    fileout.write('TEMPERATURE_VALUE={0:4.1f}\n'.format(observing_dict['TEMPERATURE_EN']))

    # --MIRROR_TEMPERATURE_KEYWORD
    # Mirror temperature in deg C.
    fileout.write('MIRROR_TEMPERATURE_KEYWORD=NONE\n')

    # --MIRROR_TEMPERATURE_VALUE
    # If MIRROR_TEMPERATURE_KEYWORD=='NONE' take this value.
    fileout.write('MIRROR_TEMPERATURE_VALUE={0:4.1f}\n'.format(observing_dict['TEMPERATURE_M1']))

    # --ELEVATION_KEYWORD
    # Elevation above sea level in m (default is Paranal: 2635. m).
    fileout.write('ELEVATION_KEYWORD=NONE\n')

    # --ELEVATION_VALUE
    # If ELEVATION_KEYWORD=='NONE' take this value.
    fileout.write('ELEVATION_VALUE={0:30f}\n'.format(observing_dict['GEOELEV']))

    # --LONGITUDE_KEYWORD
    # Longitude (default is Paranal: -70.4051 deg).
    fileout.write('LONGITUDE_KEYWORD=NONE\n')

    # --LONGITUDE_VALUE
    # If LONGITUDE_KEYWORD=='NONE' take this value.
    fileout.write('LONGITUDE_VALUE={0:9.4f}\n'.format(observing_dict['GEOLONG']))

    # --LATITUDE_KEYWORD
    # Latitude (default is Paranal: -24.6276 deg).
    fileout.write('LATITUDE_KEYWORD=NONE\n')

    # --LATITUDE_VALUE
    # If LATITUDE_KEYWORD=='NONE' take this value.
    fileout.write('LATITUDE_VALUE={0:9.4f}\n'.format(observing_dict['GEOLAT']))

    # --SLIT_WIDTH_KEYWORD
    # Slit width in arcsec (taken from FITS header if present).
    fileout.write('SLIT_WIDTH_KEYWORD=NONE\n')

    # --SLIT_WIDTH_VALUE
    # If SLIT_WIDTH_KEYWORD=='NONE' take this value.
    fileout.write('SLIT_WIDTH_VALUE={0:3.1f}\n'.format(molecfit_dict['slitwidth']))

    # --PIX_SCALE_KEYWORD
    # Pixel scale in arcsec (taken from this file only).
    fileout.write('PIX_SCALE_KEYWORD=NONE\n')

    # --PIX_SCALE_VALUE
    # If PIX_SCALE_KEYWORD=='NONE' take this value.
    fileout.write('PIX_SCALE_VALUE={0:4.12}\n'.format(molecfit_dict['pixelscale']))

    # --REFERENCE_ATMOSPHERIC
    # Reference atmospheric profile. Possible values:
    # - equ_atm (default; equatorial atmosphere, valid for Paranal);
    # - tro_atm (tropical atmosphere);
    # - std_atm (standard atmosphere);
    # - Other file located in :
    #   ({TELLURICCORR_DATA_PATH}/profiles/mipas/).
    fileout.write('REFERENCE_ATMOSPHERIC=equ.fits\n')

    # --GDAS_PROFILE
    # Specific GDAS-like input profile (P[hPa] HGT[m] T[K] RELHUM[%])
    # (if starts with /, absolute path, otherwise path relative to basedir).
    # In the case of 'none', no GDAS profiles will be considered.
    # The value 'auto' performs an automatic retrieval.
    fileout.write('GDAS_PROFILE=auto\n')

    # --LAYERS
    # Grid of layer heights for merging ref_atm and GDAS profile.
    # Fixed grid = CPL_TRUE and natural grid = CPL_FALSE.
    fileout.write('LAYERS=TRUE\n')

    # --EMIX
    # Upper mixing height in km for considering data of a local meteo station.
    # If emix is below geoelev, rhum, pres, and temp are not used
    # for modifying the corresponding profiles.
    fileout.write('EMIX=5.0\n')

    # --PWV
    # PWV value in mm for the input water vapor profile.
    # The merged profile composed of ref_atm, GDAS, and local meteo data
    # will be scaled to this value if pwv > 0 (default: -1 -> no scaling).
    fileout.write('PWV=-1.0\n')

    # --LNFL_LINE_DB
    # File name of the line list (must be stored in the directory :
    # ({TELLURICCORR_DATA_PATH}/hitran/).
    if molecfit_dict.get('aer_version', False):
        fileout.write('LNFL_LINE_DB=aer_v_{0:2.1f}\n'.format(molecfit_dict['aer_version']))
    else:
        fileout.write('LNFL_LINE_DB=aer_v_3.8\n')

    # --LNFL_LINE_DB_FORMAT
    # Format of the line file: gives the length in terms of characters per line.
    fileout.write('LNFL_LINE_DB_FORMAT=100.0\n')

    # --LBLRTM_ICNTNM
    # Continua and Rayleigh extinction [0,1,2,3,4,5].
    fileout.write('LBLRTM_ICNTNM=5\n')

    # --LBLRTM_IAERSL
    # Aerosols [0,1].
    fileout.write('LBLRTM_IAERSL=0\n')

    # --LBLRTM_MPTS
    # Number of optical depth values.
    fileout.write('LBLRTM_MPTS=5\n')

    # --LBLRTM_NPTS
    # Number of values for each panel.
    fileout.write('LBLRTM_NPTS=5\n')

    # --LBLRTM_V1
    # Beginning wavenumber value for the calculation.
    fileout.write('LBLRTM_V1=1.9\n')

    # --LBLRTM_V2
    # Ending wavenumber value for the calculation.
    fileout.write('LBLRTM_V2=2.4\n')

    # --LBLRTM_SAMPLE
    # Number of sample points per mean halfwidth [between 1 to 4, default=4].
    fileout.write('LBLRTM_SAMPLE=4\n')

    # --LBLRTM_ALFAL0
    # Average collision broadened halfwidth [cm-1/atm].
    fileout.write('LBLRTM_ALFAL0=0.0\n')

    # --LBLRTM_AVMASS
    # Average molecular mass [amu] for Doppler halfwidth.
    fileout.write('LBLRTM_AVMASS=0.0\n')

    # --LBLRTM_DPTMIN
    # Minimum molecular optical depth below which lines will be rejected.
    fileout.write('LBLRTM_DPTMIN=0.0002\n')

    # --LBLRTM_DPTFAC
    # Factor multiplying molecular continuum optical depth.
    fileout.write('LBLRTM_DPTFAC=0.001\n')

    # --LBLRTM_TBOUND
    # Temperature of boundary [K].
    fileout.write('LBLRTM_TBOUND=0.0\n')

    # --LBLRTM_SREMIS1
    # Emissivity coefficient 1.
    fileout.write('LBLRTM_SREMIS1=0.0\n')

    # --LBLRTM_SREMIS2
    # Emissivity coefficient 2.
    fileout.write('LBLRTM_SREMIS2=0.0\n')

    # --LBLRTM_SREMIS3
    # Emissivity coefficient 3.
    fileout.write('LBLRTM_SREMIS3=0.0\n')

    # --LBLRTM_SRREFL1
    # Reflectivity coefficient 1.
    fileout.write('LBLRTM_SRREFL1=0.0\n')

    # --LBLRTM_SRREFL2
    # Reflectivity coefficient 2.
    fileout.write('LBLRTM_SRREFL2=0.0\n')

    # --LBLRTM_SRREFL3
    # Reflectivity coefficient 3.
    fileout.write('LBLRTM_SRREFL3=0.0\n')

    # --LBLRTM_MODEL
    # Atmospheric profile [0,1,2,3,4,5,6].
    fileout.write('LBLRTM_MODEL=0\n')

    # --LBLRTM_ITYPE
    # Type of path [1,2,3].
    fileout.write('LBLRTM_ITYPE=3\n')

    # --LBLRTM_NOZERO
    # Zeroing of small amounts of absorbers [0,1].
    fileout.write('LBLRTM_NOZERO=0\n')

    # --LBLRTM_NOPRNT
    # Do not print output? [0,1].
    fileout.write('LBLRTM_NOPRNT=0\n')

    # --LBLRTM_IPUNCH
    # Write out layer data to TAPE7 [0,1].
    fileout.write('LBLRTM_IPUNCH=0\n')

    # --LBLRTM_RE
    # Radius of earth [km].
    fileout.write('LBLRTM_RE=0.0\n')

    # --LBLRTM_HSPACE
    # Altitude definition for space [km].
    fileout.write('LBLRTM_HSPACE=120.0\n')

    # --LBLRTM_H2
    # Upper height limit [km].
    fileout.write('LBLRTM_H2=0.0\n')

    # --LBLRTM_RANGE
    # Length of a straight path from H1 to H2 [km].
    fileout.write('LBLRTM_RANGE=0.0\n')

    # --LBLRTM_BETA
    # Earth centered angle from H1 to H2 [degrees].
    fileout.write('LBLRTM_BETA=0.0\n')

    # --LBLRTM_LEN
    # Path length [0,1].
    fileout.write('LBLRTM_LEN=0\n')

    # --LBLRTM_HOBS
    # Height of observer.
    fileout.write('LBLRTM_HOBS=0.0\n')

    # --LBLRTM_AVTRAT
    # Maximum Voigt width ratio across a layer.
    fileout.write('LBLRTM_AVTRAT=2.0\n')

    # --LBLRTM_TDIFF1
    # Maximum layer temperature difference at ALTD1 [K].
    fileout.write('LBLRTM_TDIFF1=5.0\n')

    # --LBLRTM_TDIFF2
    # Maximum layer temperature difference at ALTD2 [K].
    fileout.write('LBLRTM_TDIFF2=8.0\n')

    # --LBLRTM_ALTD1
    # Altitude of TDIFF1 [km].
    fileout.write('LBLRTM_ALTD1=0.0\n')

    # --LBLRTM_ALTD2
    # Altitude of TDIFF2 [km].
    fileout.write('LBLRTM_ALTD2=0.0\n')

    # --LBLRTM_DELV
    # Number of wavenumbers [cm-1] per major division.
    fileout.write('LBLRTM_DELV=1.0\n')

    # --EXPERT_MODE
    # If set to true, will check if TAG[INIT_FIT_PARAMETERS] points to a fits file
    # with a bintable of parameter values to use as initial values for the
    # fitting process.
    fileout.write('EXPERT_MODE=FALSE\n')

    # --CHIP_EXTENSIONS
    # Flag that determines if image extensions are to be treated as independent
    # science data to be fitted for independently or as CHIP specific subranges
    # of a single observation to be fitted for as a single combined spectrum.
    # Value = TRUE implies to treat as CHIPS to be combined. Value = FALSE
    # implies
    # to treat as independent. [FALSE].
    fileout.write('CHIP_EXTENSIONS=FALSE\n')

    #
    # End of file


    fileout.close()


def write_molecfit_input_spectrum(wave, flux, filename):
    #writing the molecfit input fits files

    good=range(len(wave))

    for i in range(5):

        medFlux=np.median(flux[good])
        madFlux=np.median(np.abs(flux[good]-medFlux))
        good=np.where(flux>(medFlux-5*madFlux))[0]

    dflux=np.std(flux[good])
    dflux=dflux*np.sqrt(flux)#to scale the error to the signal assuming poissonian noise

    #To create a table from scratch, we need to define columns first, by constructing the Column objects and their data.
    #Suppose we have two columns, the first containing strings, and the second containing floating point numbers:

    col1 = fits.Column(name='lambda', format='E', array=wave*1e-4)
    col2 = fits.Column(name='flux', format='E', array=flux)
    col3 = fits.Column(name='dflux', format='E', array=dflux)
    col4 = fits.Column(name='mask', format='I', array=np.zeros(len(good)))

    #Now, create a new binary table HDU object by using the BinTableHDU.from_columns() function:
    #hdu = fits.BinTableHDU.from_columns(cols)
    hdu = fits.BinTableHDU.from_columns([col1, col2, col3, col4])

    #The data structure used to represent FITS tables is called a FITS_rec and is derived from the numpy.recarray interface.
    #When creating a new table HDU the individual column arrays will be assembled into a single FITS_rec array.

    #Now you may write this new table HDU directly to a FITS file like so:
    hdu.writeto(filename, overwrite=True)


def write_calctrans_par(filename_par):
    fileout = open(filename_par, 'w')


    # File: /archive/molecfitRootDataDir/calctrans.rc
    #
    # Note: This configuration file has been automatically
    #       generated by the esorex (v3.13.3) program.
    #
    # Date: 17-Jan-2021 18:30:01
    #
    #

    # --USE_ONLY_INPUT_PRIMARY_DATA
    # Value=TRUE implies that only the fits primary contains the input science flux
    # data.
    # Value=FALSE implies that the fits extensions also contains input science
    # flux data.
    fileout.write('USE_ONLY_INPUT_PRIMARY_DATA=FALSE\n')

    # --USE_DATA_EXTENSION_AS_DFLUX
    # Only valid if USE_ONLY_INPUT_PRIMARY_DATA=TRUE. The fits extension index that
    # contains the
    # errors of the science flux data (DFLUX). A value of 0 implies that there is
    # no DFLUX.
    fileout.write('USE_DATA_EXTENSION_AS_DFLUX=0\n')

    # --USE_DATA_EXTENSION_AS_MASK
    # Only valid if USE_ONLY_INPUT_PRIMARY_DATA=TRUE. The fits extension index that
    # contains the
    # mask associated with the science flux data. A value of 0 implies that there
    # is no mask data.
    fileout.write('USE_DATA_EXTENSION_AS_MASK=0\n')

    # --USE_INPUT_KERNEL
    # If TRUE, use the kernel library if it is provided.
    fileout.write('USE_INPUT_KERNEL=TRUE\n')

    # --CALCTRANS_MAPPING_KERNEL
    # Mapping 'SCIENCE' - 'CALCTRANS_KERNEL_LIBRARY' [string with ext_number comma
    # separated (int)] :
    # If set to NULL, check if the TAG[CALCTRANS_MAPPING_KERNEL] FITS BINTABLE
    # values is provided.
    # The FITS BINTABLE have to one column [KERNEL_LIBRARY_EXT].
    fileout.write('CALCTRANS_MAPPING_KERNEL=NULL\n')

    # --MAPPING_ATMOSPHERIC
    # Mapping 'SCIENCE' - 'ATM_PARAMETERS' input [string with ext_number comma
    # separated (int)] :
    # If set to NULL, check if the TAG[MAPPING_ATMOSPHERIC] FITS BINTABLE value
    # is provided.
    # The FITS BINTABLE have to one column [ATM_PARAMETERS_EXT].
    fileout.write('MAPPING_ATMOSPHERIC=0,1\n')

    # --MAPPING_CONVOLVE
    # Mapping 'LBLRTM_RESULTS' - 'TELLURIC_CORR' output [string with ext_number
    # comma separated (int)] :
    # If set to NULL, check if the TAG[MAPPING_CONVOLVE] FITS BINTABLE value is
    # provided.
    # The FITS BINTABLE have to one column [LBLRTM_RESULTS_EXT].
    fileout.write('MAPPING_CONVOLVE=0,1\n')

    # --CHIP_EXTENSIONS
    # Flag that determines if image extensions are to be treated as independant
    # science data to be fitted for independently or as CHIP specific subranges
    # of a single observation to be fitted for as a single combined spectrum.
    # Value = TRUE implies to treat as CHIPS to be combined. Value = FALSE
    # implies
    # to treat as independent. [FALSE].
    fileout.write('CHIP_EXTENSIONS=FALSE\n')

    #
    # End of file

    fileout.close()

def write_correct_par(filename_par):

    fileout = open(filename_par, 'w')

    # File: /archive/molecfitRootDataDir/correct.rc
    #
    # Note: This configuration file has been automatically
    #       generated by the esorex (v3.13.3) program.
    #
    # Date: 17-Jan-2021 18:30:12
    #
    #

    # --USE_ONLY_INPUT_PRIMARY_DATA
    # Value=TRUE implies that only the fits primary contains the input science flux
    # data.
    # Value=FALSE implies that the fits extensions also contains input science
    # flux data.
    fileout.write('USE_ONLY_INPUT_PRIMARY_DATA=FALSE\n')

    # --USE_DATA_EXTENSION_AS_DFLUX
    # Only valid if USE_ONLY_INPUT_PRIMARY_DATA=TRUE. The fits extension index that
    # contains the
    # errors of the science flux data (DFLUX). A value of 0 implies that there is
    # no DFLUX.
    fileout.write('USE_DATA_EXTENSION_AS_DFLUX=0\n')

    # --USE_DATA_EXTENSION_AS_MASK
    # Only valid if USE_ONLY_INPUT_PRIMARY_DATA=TRUE. The fits extension index that
    # contains the
    # mask associated with the science flux data. A value of 0 implies that there
    # is no mask data.
    fileout.write('USE_DATA_EXTENSION_AS_MASK=0\n')

    # --SUPPRESS_EXTENSION
    # Suppress arbitrary filename extension : TRUE (apply) or FALSE (don't apply).
    fileout.write('SUPPRESS_EXTENSION=FALSE\n')

    # --MAPPING_CORRECT
    # Mapping 'SCIENCE' - 'TELLURIC_CORR' [string with ext_number comma separated
    # (int)] :
    # If set to NULL, check if the TAG[MAPPING_CORRECT] FITS BINTABLE value is
    # provided.
    # The FITS BINTABLE have to one column [TELLURIC_CORR_EXT].
    fileout.write('MAPPING_CORRECT=0,1\n')

    # --CHIP_EXTENSIONS
    # Flag that determines if image extensions are to be treated as independant
    # science data to be fitted for independently or as CHIP specific subranges
    # of a single observation to be fitted for as a single combined spectrum.
    # Value = TRUE implies to treat as CHIPS to be combined. Value = FALSE
    # implies
    # to treat as independent. [FALSE].
    fileout.write('CHIP_EXTENSIONS=FALSE\n')

    #
    # End of file

    fileout.close()
