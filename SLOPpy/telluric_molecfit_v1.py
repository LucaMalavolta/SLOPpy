from __future__ import print_function, division
from SLOPpy.subroutines.common import *
from SLOPpy.subroutines.spectral_subroutines import *
from SLOPpy.subroutines.io_subroutines import *
from SLOPpy.subroutines.fit_subroutines import *
from SLOPpy.subroutines.plot_subroutines import *
from SLOPpy.subroutines.shortcuts import *

__all__ = ["compute_telluric_molecfit_v1",
           "plot_telluric_molecfit_v1"]

def compute_telluric_molecfit_v1(config_in):
    """
    Lazy workaround
    :param config_in:
    :param kwargs:
    :return:
    """

    print('UNTESTED PRROCEDUE - I QUIT')
    quit()
    night_dict = from_config_get_nights(config_in)
    instrument_dict = from_config_get_instrument(config_in)

    for night in night_dict:

        instrument_name = night_dict[night]['instrument']
        template_dict = instrument_dict[instrument_name]['telluric_template']
        print()
        print("compute_telluric_molecfit                  Night: ", night)
        print()

        try:
            telluric = load_from_cpickle('telluric', config_in['output'], night)
            continue
        except:
            print("  No telluric correction file found, computing now ")
            print()

        print('    instrument :', instrument_name)
        print('    template   :', template_dict['file'])
        print('    fit_range  :', template_dict['fit_range'])
        print()

        """ Retrieving the list of observations"""
        lists = load_from_cpickle('lists', config_in['output'], night)

        """ Retrieving the observations"""
        calib_data = load_from_cpickle('calibration_fibA', config_in['output'], night)
        input_data = retrieve_observations(config_in['output'], night, lists['observations'], use_telluric=False)
        observational_pams = load_from_cpickle('observational_pams', config_in['output'], night)

        processed = {
            'subroutine': 'telluric_molecfit',
            'n_orders': 0,
            'n_pixels': 0,
        }

        telluric = {
            'subroutine': 'telluric_molecfit',
            'reference_frame': 'observer'
        }


        processed['airmass_ref'] = 0.000
        processed['telluric'] = {}
        processed['rebin'] = {}

        """
        Molecfit works on pixel grid, so we must ensure that the spectra are rebinned always on the same wavelength 
        scale and same wavelength step. We use local arrays for this purpose 
        """

        rebin_step_unit = 0.01000000
        processed['rebin']['wave'] = np.arange(input_data['coadd']['wavelength_range'][0],
                                               input_data['coadd']['wavelength_range'][1],
                                               rebin_step_unit,
                                               dtype=np.double)

        processed['rebin']['size'] = np.size(processed['rebin']['wave'])
        processed['rebin']['step'] = np.ones(processed['rebin']['size'], dtype=np.double) * rebin_step_unit

        processed['rebin'] = {
            'wave': input_data['coadd']['wave'],
            'size': input_data['coadd']['size'],
            'step': input_data['coadd']['step'],
        }

        print('  Writing data and configuration files for molecfit+calctrans')
        print()
        """
        We store all the molecfit files in a subdirectory
        We save the path of the main directory to a temporary file
        """

        os.system('mkdir -p molecfit_'+night)
        os.system('mkdir -p molecfit_'+night + '/output/')

        # There must be a more elegant way to do this, but I'm, not aware of it
        for n_obs, obs in enumerate(lists['observations']):

            processed[obs] = {
                'n_orders': input_data[obs]['n_orders'],
                'n_pixels': input_data[obs]['n_pixels']
            }



            """ e2ds spectra are rescaled and then rebinned while keeping them in the Observer Reference Frame"""

            processed[obs]['e2ds_rescaling'], processed[obs]['e2ds_rescaled'], processed[obs]['e2ds_rescaled_err'] = \
                perform_rescaling(input_data[obs]['wave'],
                                  input_data[obs]['e2ds'],
                                  input_data[obs]['e2ds_err'],
                                  observational_pams['wavelength_rescaling'])

            processed[obs]['rebin_ORF'] = \
            rebin_2d_to_1d(input_data[obs]['wave'],
                           input_data[obs]['step'],
                           processed[obs]['e2ds_rescaled'],
                           calib_data['blaze'],
                           processed['rebin']['wave'],
                           processed['rebin']['step'],
                           rv_shift=0.00)

            """ Molecfit analysis is skipped if the telluric computation has been computed already"""
            if os.path.isfile('./molecfit_'+night +'/output/'+obs+'_ORF_s1d_TAC.dat'):
                print('    molecfit+calctrans results for ' + obs + ' already available')

                continue

            """ the spectra is save onto an ASCII file in a format suitable for molecfit """
            fileout = open('./molecfit_'+night +'/'+obs+'_ORF_s1d.dat', 'w')
            for w, f in zip(processed['rebin']['wave'], processed[obs]['rebin_ORF']):
                fileout.write('{0:12.6f} {1:12.6f} \n'.format(w, f))
            fileout.close()

            """
            processed[obs]['rebin_SRF'] = \
            rebin_2d_to_1d(input_data[obs]['wave'],
                           input_data[obs]['step'],
                           processed[obs]['e2ds_rescaled'],
                           calib_data['blaze'],
                           processed['rebin']['wave'],
                           processed['rebin']['step'],
                           rv_shift = observational_pams[obs]['rv_shift_ORF2SRF'])

            fileout = open('./molecfit_'+night +'/'+obs+'_SRF_s1d.dat','w')
            for w, f in zip(processed['rebin']['wave'], processed[obs]['rebin_SRF']):
                fileout.write('{0:12.6f} {1:12.6f} \n'.format(w, f))
            fileout.close()
            """

            # TODO: input from configuration file for molecfit installation path
            bash_script = open('./molecfit_'+night +'/molecfit_exec_' + obs + '.source', 'w')
            bash_script.write('#!/bin/bash \n')

            bash_script.write('echo  "   " executing molecfit+calctrans on '+obs+' \n')
            bash_script.write('/usr/local/eso/bin/molecfit '+obs+'.par > ' + obs +'_molecfit.log\n')
            bash_script.write('/usr/local/eso/bin/calctrans '+obs+'.par > ' + obs +'_calctrans.log\n')
            bash_script.close()

            fileout = open('./molecfit_'+night +'/'+obs+'.par', 'w')
            fileout.write("### Driver for MOLECFIT\n")

            # user working directory only important for REFLEX workflow and GUI
            # not used by molecfit itself.
            fileout.write("user_workdir:./\n")

            ## INPUT DATA
            # Data file name (path relative to the current directory or absolute path)
            fileout.write("filename: " + obs +"_ORF_s1d.dat\n")

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
            fileout.write("wrange_include: include_"+night+".dat\n")

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
            fileout.write("output_name: "+ obs + "\n")

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
            fileout.write("ftol: " + input_data[obs]['molecfit']['ftol'] + "\n")

            # Relative parameter convergence criterion
            fileout.write("xtol: " + input_data[obs]['molecfit']['xtol'] + "\n")

            ## MOLECULAR COLUMNS

            # List of molecules to be included in the model
            # (default: 'H2O', N_val: nmolec)
            molecules_list = "list_molec:"
            for mol in input_data[obs]['molecfit']['molecules']:
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
            fileout.write("cont_n: {0:1.0f}".format(input_data[obs]['molecfit']['cont_n']) + "\n")

            # Initial constant term for continuum fit (valid for all fit ranges)
            # (emission spectrum: about 1 for correct flux_unit)
            fileout.write("cont_const: {0:1.0f}".format(input_data[obs]['molecfit']['cont_const']) + "\n")

            ## WAVELENGTH SOLUTION

            # Refinement of wavelength solution using a polynomial of degree wlc_n
            fileout.write("fit_wlc: 1\n")

            # Polynomial degree of the refined wavelength solution
            fileout.write("wlc_n: {0:1.0f}".format(input_data[obs]['molecfit']['wlc_n']) + "\n")

            # Initial constant term for wavelength correction (shift relative to half
            # wavelength range)
            fileout.write("wlc_const: {0:1.0f}".format(input_data[obs]['molecfit']['wlc_const']) + "\n")

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
            fileout.write("res_gauss: {0:3.1f}".format(input_data[obs]['molecfit']['res_gauss']) + "\n")

            # Fit resolution by Lorentzian -- 1 = yes; 0 = no
            fileout.write("fit_res_lorentz: 0\n")

            # Initial value for FWHM of Lorentzian in pixels
            fileout.write("res_lorentz: 0.0\n")

            # Size of Gaussian/Lorentzian/Voigtian kernel in FWHM
            fileout.write("kernfac: {0:3.0f}".format(input_data[obs]['molecfit']['kernfac']) + "\n")

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
            fileout.write("obsdate: {0:13.5f}".format(input_data[obs]['MJD']) + "\n")
            fileout.write("obsdate_key: NONE\n")

            # UTC in s
            fileout.write("utc: {0:8.1f}".format(input_data[obs]['UTC']) + "\n")
            fileout.write("utc_key: NONE\n")

            # Telescope altitude angle in deg
            fileout.write("telalt: {0:13.5f}".format(input_data[obs]['ELEVATION']) + "\n")
            fileout.write("telalt_key: NONE\n")

            # Humidity in %
            fileout.write("rhum: {0:13.5f}".format(input_data[obs]['HUMIDITY']) + "\n")
            fileout.write("rhum_key: NONE\n")

            # Pressure in hPa
            fileout.write("pres: {0:5.1f}".format(input_data[obs]['PRESSURE']) + "\n")
            fileout.write("pres_key: NONE\n")

            # Ambient temperature in deg C
            # temp: 15.0
            fileout.write("temp: {0:4.1f}".format(input_data[obs]['TEMPERATURE_EN']) + "\n")
            fileout.write("temp_key: NONE\n")

            # Mirror temperature in deg C
            # m1temp: 15.0
            fileout.write("m1temp: {0:4.1f}".format(input_data[obs]['TEMPERATURE_M1']) + "\n")
            fileout.write("m1temp_key: NONE\n")

            # Elevation above sea level in m (default is Paranal: 2635m)
            # geoelev: 2387.2
            fileout.write("geoelev: {0:4.0f}".format(input_data[obs]['GEOELEV']) + "\n")
            fileout.write("geoelev_key: NONE\n")

            # Longitude (default is Paranal: -70.4051)
            # longitude: -17.889
            fileout.write("longitude: {0:9.4f}".format(input_data[obs]['GEOLONG']) + "\n")
            fileout.write("longitude_key: NONE\n")

            # Latitude (default is Paranal: -24.6276)
            # latitude: 28.754
            fileout.write("latitude: {0:9.4f}".format(input_data[obs]['GEOLAT']) + "\n")
            fileout.write("latitude_key: NONE\n")

            ## INSTRUMENTAL PARAMETERS

            # Slit width in arcsec (taken from FITS header if present)
            fileout.write("slitw: {0:3.1f}".format(input_data[obs]['molecfit']['slitwidth']) + "\n")
            fileout.write("slitw_key: NONE\n")

            # Pixel scale in arcsec (taken from this file only)
            fileout.write("pixsc: {0:4.2f}".format(input_data[obs]['molecfit']["pixelscale"]) + "\n")
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

            os.system('cd  molecfit_' + night + '/ && . ./molecfit_exec_' + obs + '.source')

        print()
        print('  molecfit+calcatrans completed')

        for n_obs, obs in enumerate(lists['observations']):

            telluric[obs] = {}

            """ Loading the telluric spectrum from the output directory of molecfit """
            telluric_molecfit = np.genfromtxt('./molecfit_'+night +'/output/'+obs+'_ORF_s1d_TAC.dat', usecols=2)
            """ rebinning onto the e2ds wave scale"""

            telluric[obs]['spectrum'] = \
                rebin_1d_to_2d(processed['rebin']['wave'],
                               processed['rebin']['step'],
                               telluric_molecfit,
                               input_data[obs]['wave'],
                               input_data[obs]['step'],
                               preserve_flux=False)

            try:
                telluric[obs]['spectrum'] = np.nan_to_num(nan=1.0, posinf=1.0, neginf=1.0)
            except:
                temp = np.isfinite(telluric[obs]['spectrum'])
                telluric[obs]['spectrum'][temp] = 1.0


            telluric[obs]['airmass'] = input_data[obs]['AIRMASS']


            " for compatibilty to some plots, even if it doesn't make any sense"
            telluric[obs]['airmass_ref'] = 0.000
            telluric[obs]['spectrum_noairmass'] = np.power(telluric[obs]['spectrum'],
                                    telluric[obs]['airmass_ref'] - input_data[obs]['AIRMASS'])
            telluric[obs]['null'] = telluric[obs]['spectrum_noairmass'] < 0.001
            telluric[obs]['spectrum_noairmass'][telluric[obs]['null']] = 1.0
            # we just copy the spectrum file, it's it's a model itself
            telluric[obs]['spline'] = telluric[obs]['spectrum'].copy()

            processed[obs]['e2ds_corrected'] = processed[obs]['e2ds_rescaled'] / telluric[obs]['spectrum']
            processed[obs]['e2ds_corrected_err'] = processed[obs]['e2ds_rescaled_err'] / telluric[obs]['spectrum']

        save_to_cpickle('telluric', telluric, config_in['output'], night)
        save_to_cpickle('telluric_processed', processed, config_in['output'], night)

        print()
        print("Night ", night, " completed")

        #
        #""" After being rescaled for the proper factor, the template telluric spectrum is rebinned onto the 2D
        #scale of the observations """
        #
        #telluric['template']['rebinned']['flux'] = \
        #    rebin_1d_to_2d(telluric['template']['input']['wave'],
        #                   telluric['template']['input']['step'],
        #                   telluric['template']['input']['flux'],
        #                   telluric['template']['rebinned']['wave'],
        #                   telluric['template']['rebinned']['step'],
        #                   preserve_flux=False)
        #
        #telluric['template']['rebinned']['ferr'] = \
        #    rebin_1d_to_2d(telluric['template']['input']['wave'],
        #                   telluric['template']['input']['step'],
        #                   telluric['template']['input']['ferr'],
        #                   telluric['template']['rebinned']['wave'],
        #                   telluric['template']['rebinned']['step'],
        #                   preserve_flux=False,
        #                   is_error=True)
        #
        #
        #sel_out_of_range = ~((telluric['template']['rebinned']['wave'] > telluric['template']['input']['range'][0]+1.) \
        #                    & (telluric['template']['rebinned']['wave'] < telluric['template']['input']['range'][1]-1.))
        #telluric['template']['rebinned']['flux'][sel_out_of_range] = 1.
        #telluric['template']['rebinned']['ferr'][sel_out_of_range] = 0.1
        #
        #processed['telluric']['spectrum_noairmass'] = \
        #    (telluric['template']['rebinned']['flux'] - 1.) * telluric_factor + 1.0
        #
        #telluric['airmass_ref'] = processed['airmass_ref']
        #
        #for obs in lists['observations']:
        #    """ Correction of telluric lines for the average airmass value, following Wyttenbach et al. 2015 """
        #    processed[obs]['e2ds_corrected'] = processed[obs]['e2ds_rescaled'] / \
        #                                            np.power(processed['telluric']['spectrum_noairmass'],
        #                                                     input_data[obs]['AIRMASS'] - processed['airmass_ref'])
        #    processed[obs]['e2ds_corrected_err'] = processed[obs]['e2ds_rescaled_err'] / \
        #                                            np.power(processed['telluric']['spectrum_noairmass'],
        #                                                     input_data[obs]['AIRMASS'] - processed['airmass_ref'])
        #
        #for obs in lists['observations']:
        #    # Correction of telluric lines
        #
        #    telluric[obs] = {}
        #
        #    telluric[obs]['spectrum_noairmass'] = processed['telluric']['spectrum_noairmass']
        #
        #    telluric[obs]['airmass'] = input_data[obs]['AIRMASS']
        #    telluric[obs]['airmass_ref'] = processed['airmass_ref']
        #
        #    """ Set anomalosly low point to one (e.g. when the template is not computed)"""
        #    telluric[obs]['null'] = telluric[obs]['spectrum_noairmass'] < 0.001
        #    telluric[obs]['spectrum_noairmass'][telluric[obs]['null']] = 1.0
        #
        #    telluric[obs]['spectrum'] = np.power(processed['telluric']['spectrum_noairmass'],
        #                               input_data[obs]['AIRMASS'] - processed['airmass_ref'])
        #
        #    telluric[obs]['spline_noairmass'] = telluric[obs]['spectrum_noairmass'].copy()
        #
        #    """ No need to compute the spline approximation since we are already dealing with a very high SNR template"""
        #    telluric[obs]['spline'] = np.power(telluric[obs]['spline_noairmass'],
        #                               input_data[obs]['AIRMASS'] - processed['airmass_ref'])
        #
        #    """ copy the keyword for future use"""
        #    telluric[obs]['airmass'] = input_data[obs]['AIRMASS']
        #
        #    telluric[obs]['telluric_corrected'] = processed[obs]['e2ds_corrected']
        #    telluric[obs]['telluric_corrected_err'] = processed[obs]['e2ds_corrected_err']
        #
        #    save_to_cpickle('telluric', telluric, config_in['output'], night)
        #    save_to_cpickle('telluric_processed', processed, config_in['output'], night)

        print()
        print("Night ", night, " completed")

    quit()

def plot_telluric_molecfit_v1(config_in, night_input=''):
    import matplotlib.pyplot as plt

    night_dict = from_config_get_nights(config_in)
    instrument_dict = from_config_get_instrument(config_in)
    system_dict = from_config_get_system(config_in)

    if night_input == '':
        night_list = night_dict
    else:
        night_list = np.atleast_1d(night_input)

    for night in night_list:

        #plt.scatter(rescaling_array, computed_std, c='C0', zorder=1)
        #plt.scatter(sel_factor, sel_stdev, c='C1', zorder=2)
        #plt.plot(rescaling_array, np.polyval(coeff, rescaling_array))
        #plt.plot(rescaling_array, 2*rescaling_array*coeff[0] + coeff[1] )
        #plt.plot()

        print("plot_telluric_template                     Night: ", night)

        """ Retrieving the list of observations"""
        lists = load_from_cpickle('lists', config_in['output'], night)
        observational_pams = load_from_cpickle('observational_pams', config_in['output'], night)

        """ Retrieving the analysis"""
        try:
            processed = load_from_cpickle('telluric_processed', config_in['output'], night)
            telluric = load_from_cpickle('telluric', config_in['output'], night)
        except:
            print()
            print("No telluric correction, no plots")
            continue

        colors, cmap, line_colors = make_color_array(lists, observational_pams)

        fig = plt.figure(figsize=(12, 6))
        gs = GridSpec(2, 2, width_ratios=[50, 1])

        ax1 = plt.subplot(gs[0, 0])
        ax2 = plt.subplot(gs[1, 0], sharex=ax1)
        cbax1 = plt.subplot(gs[:, 1])

        lift_spectrum = 0.25

        for i, obs in enumerate(lists['observations']):
            color_array = cmap(i / len(lists['observations']))

            _, e2ds_rescaled , _ = \
                perform_rescaling(processed[obs]['wave'],
                                  processed[obs]['e2ds'],
                                  processed[obs]['e2ds_err'],
                                  observational_pams['wavelength_rescaling'])

            e2ds_rescaled_corrected_spectrum = e2ds_rescaled / telluric[obs]['spectrum']
            e2ds_rescaled_corrected_spline = e2ds_rescaled / telluric[obs]['spline']

            for order in range(0, processed[obs]['n_orders']):

                if order == 0 and i==0:
                    ax1.plot(processed[obs]['wave'][order, :],
                                e2ds_rescaled[order, :],
                                c=color_array, lw=1, alpha=0.5, label='uncorrected')
                    ax1.scatter(processed[obs]['wave'][order, :],
                                e2ds_rescaled_corrected_spectrum[order, :],
                                s=1, c=np.atleast_2d(color_array), label='corrected')
                else:
                    ax1.plot(processed[obs]['wave'][order, :],
                                e2ds_rescaled[order, :],
                                c=color_array, lw=1, alpha=0.5)
                    ax1.scatter(processed[obs]['wave'][order, :],
                                e2ds_rescaled_corrected_spectrum[order, :],
                                s=1, c=np.atleast_2d(color_array))

                #ax1.plot(processed[obs]['wave'][order, :],
                #            e2ds_rescaled[order, :]+lift_spectrum,
                #            c=color_array, lw=1, alpha=0.5)
                #ax1.scatter(processed[obs]['wave'][order, :],
                #            e2ds_rescaled_corrected_spline[order, :]+lift_spectrum,
                #            s=1, c=np.atleast_2d(color_array))

                ax2.plot(processed[obs]['wave'][order, :],
                         telluric[obs]['spectrum'][order, :],
                         c=color_array)
                ax2.axhline(1.00, c='k')

                #ax2.plot(processed[obs]['wave'][order, :],
                #         telluric[obs]['spline'][order, :]+lift_spectrum,
                #         c=color_array)
                #ax2.axhline(1.00+lift_spectrum, c='k')

        #ax2.plot(input_data['coadd']['wave'],telluric['stellarRF']['spline_eval']+0.1,c='k')
        #ax2.scatter(input_data['coadd']['wave'],telluric['stellarRF']['spectrum']+0.1,c='r', s=2)

        ax1.legend(loc=3)
        ax1.set_title('Night: ' + night)

        ax2.set_xlabel('$\lambda$ [$\AA$]')

        try:
            instrument = night_dict[night]['instrument']
            comparison_file = config_in['instruments'][instrument]['telluric_comparison']
            comparison_data = np.genfromtxt(comparison_file, skip_header=1)
            if comparison_data[0,0]<1000.0:
                nm2Ang = 10.
            else:
                nm2Ang = 1.
            ax1.plot(comparison_data[:, 0]*nm2Ang, comparison_data[:, 1], c='C0', zorder=1000)
            ax2.plot(comparison_data[:, 0]*nm2Ang, comparison_data[:, 1], c='C0', zorder=1000)
        except:
            pass

        sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=colors[0], vmax=colors[-1]))
        sm.set_array([])  # You have to set a dummy-array for this to work...
        cbar = plt.colorbar(sm, cax=cbax1)
        cbar.set_label('BJD - 2450000.0')
        fig.subplots_adjust(wspace=0.05, hspace=0.4)
        plt.show()