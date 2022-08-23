from __future__ import print_function, division
from SLOPpy.subroutines.common import *
from SLOPpy.subroutines.spectral_subroutines import *
from SLOPpy.subroutines.io_subroutines import *
from SLOPpy.subroutines.fit_subroutines import *
from SLOPpy.subroutines.plot_subroutines import *
from SLOPpy.subroutines.shortcuts import *

__all__ = ["compute_telluric_template_alternative", "plot_telluric_template_alternative"]

def compute_telluric_template_alternative(config_in):

    compute_telluric_template_alternative_routine(config_in,
                                n_iterations=1,
                                use_berv=False,
                                use_reference_airmass=False,
                                use_template=True,
                                subroutine_name='telluric_template')


def compute_telluric_template_alternative_routine(config_in, **kwargs):
    """
    Lazy workaround
    :param config_in:
    :param kwargs:
    :return:
    """

    night_dict = from_config_get_nights(config_in)
    instrument_dict = from_config_get_instrument(config_in)

    shared_data = load_from_cpickle('shared', config_in['output'])

    for night in night_dict:

        instrument_name = night_dict[night]['instrument']

        print()
        print("compute_telluric_template                  Night: ", night)
        print()

        try:
            telluric = load_from_cpickle('telluric', config_in['output'], night)
            continue
        except:
            print("  No telluric correction file found, computing now ")
            print()



        """ Retrieving the list of observations"""
        lists = load_from_cpickle('lists', config_in['output'], night)

        """ Retrieving the observations"""
        calib_data = load_from_cpickle('calibration_fibA', config_in['output'], night)
        input_data = retrieve_observations(config_in['output'], night, lists['observations'], use_telluric=False)
        observational_pams = load_from_cpickle('observational_pams', config_in['output'], night)

        processed = {
            'subroutine': kwargs['subroutine_name'],
            'n_orders': 0,
            'n_pixels': 0,
            'telluric': {}
        }

        telluric = {
            'subroutine': kwargs['subroutine_name'],
            'reference_frame': 'observer',
            'template': {},
            'linear' : {}
        }

        # There must be a more elegant way to do this, but I'm, not aware of it
        """ computation of the rescaled spectra, for later use in the analysis and in the plotting subroutines"""
        for obs in lists['observations']:
            processed[obs] = {
                'n_orders': input_data[obs]['n_orders'],
                'n_pixels': input_data[obs]['n_pixels']
            }

            """ for plotting purpose only"""
            processed[obs]['wave'] = input_data[obs]['wave']
            processed[obs]['e2ds'] = input_data[obs]['e2ds']
            processed[obs]['e2ds_err'] = input_data[obs]['e2ds_err']

            processed[obs]['e2ds_rescaling'], processed[obs]['e2ds_rescaled'], processed[obs]['e2ds_rescaled_err'] = \
                perform_rescaling(input_data[obs]['wave'],
                                  input_data[obs]['e2ds'],
                                  input_data[obs]['e2ds_err'],
                                  observational_pams['wavelength_rescaling'])

        """ Reference airmass for iterative correction of airmass - disabled here"""
        processed['airmass_ref'] = 0.000

        """ 
        Definiton of the wavelength scale for the output telluric spectrum
        We assume that the wavelength solution did not change during the night 
        """
        obs_reference = lists['observations'][0]
        telluric['rebinned'] = {
                'wave': input_data[obs_reference]['wave'],
                'step': input_data[obs_reference]['step'],
                'flux': np.ones(input_data[obs_reference]['step'].shape),
                'ferr': np.ones(input_data[obs_reference]['step'].shape)*0.001
        }

        processed['telluric']['spectrum_noairmass'] = telluric['rebinned']['flux'].copy()
        processed['telluric']['spectrum_noairmass_err'] = telluric['rebinned']['ferr'].copy()
        processed['n_orders'] = input_data[obs_reference]['orders']
        processed['n_pixels'] = input_data[obs_reference]['wave_size']

        """ 
        Computation of the rescaling factor for the telluric template. This factor has no physical value, it's just 
        the ratio of the observed telluric features (in a normalized spectra) to the input template.
        """
        try:

            if 'telluric_template' not in instrument_dict[instrument_name]:
                raise MissingKeywordException()

            template_dict = instrument_dict[instrument_name]['telluric_template']

            if template_dict['fit_range'][0] > shared_data['coadd']['wavelength_range'][1] or \
                    template_dict['fit_range'][1] < shared_data['coadd']['wavelength_range'][0]:
                raise OutOfRangeException()

            print('    rescaled telluric template')
            print('    instrument :', instrument_name)
            print('    template   :', template_dict['file'])
            print('    fit_range  :', template_dict['fit_range'])
            print()

            """ Retrieving the template data"""
            telluric_template_data = np.genfromtxt(template_dict['file'])

            telluric['template']['input'] = {
                'range': [np.amin(telluric_template_data[:, 0]), np.amax(telluric_template_data[-1, 0])],
                'wave': telluric_template_data[:, 0],
                'flux': telluric_template_data[:, 1],
                'ferr': telluric_template_data[:, 2],
                'step': telluric_template_data[:, 3]
            }

            processed['template'] = {
                'array_wave': [],
                'array_data': {},
                'array_move': {}
            }

            """ Again we assume that the wavelength solution did not change during the night  """

            tel_selection = (input_data[obs_reference]['wave'] > template_dict['fit_range'][0]) & \
                            (input_data[obs_reference]['wave'] < template_dict['fit_range'][1])
            processed['template']['wave_sel'] = input_data[obs_reference]['wave'][tel_selection]
            processed['template']['step_sel'] = input_data[obs_reference]['step'][tel_selection]

            for obs in lists['telluric']:
                processed['template'][obs] = processed[obs]['e2ds_rescaled'][tel_selection]
                """ saving the wavelength array for plotting purpose"""
                processed['template']['array_wave'].extend(processed['template']['wave_sel'])

            rescaling_array = np.arange(0.05, 2.0, 0.05)
            computed_std = np.zeros(len(rescaling_array))

            for n_rescaling_factor, rescaling_factor in enumerate(rescaling_array):

                """ 
                Template spectrum is rebinned onto the observations wavelength scale, using only the wavelength 
                range selected for the computation of the rescaling factor 
                """
                template_rebinned_flux = \
                    rebin_1d_to_1d(telluric['template']['input']['wave'],
                                   telluric['template']['input']['step'],
                                   telluric['template']['input']['flux'],
                                   processed['template']['wave_sel'],
                                   processed['template']['step_sel'],
                                   preserve_flux=False)

                """ Saving the outcome to dictionary """
                template_rebinned_flux -= 1.00
                template_rebinned_flux *= rescaling_factor
                template_rebinned_flux += 1.00

                e2ds_corrected = []

                for obs in lists['telluric']:

                    e2ds_corrected.extend(processed['template'][obs] /
                                            np.power(template_rebinned_flux, input_data[obs]['AIRMASS']))

                computed_std[n_rescaling_factor] = np.nanstd(e2ds_corrected)

                #print(n_rescaling_factor, rescaling_factor, computed_std[n_rescaling_factor])
                #plt.scatter(processed['template']['array_wave'], e2ds_corrected)
                #plt.plot(processed['template']['wave_sel'], template_rebinned_flux)
                #plt.show()


                label_dict = '{0}'.format(rescaling_factor)
                processed['template']['array_data'][label_dict] = e2ds_corrected
                processed['template']['array_move'][label_dict] = rescaling_factor


            #plt.scatter(rescaling_array, computed_std)
            #plt.show()

            """ Selection of the rescaling factor with the lowest scatter """
            ind_factor = np.argmin(computed_std)
            ind_range = 3
            if ind_factor < ind_range:
                sel_factor = rescaling_array[0:ind_factor+ind_range]
                sel_stdev = computed_std[0:ind_factor+ind_range]
            elif ind_factor > len(rescaling_array) - ind_range:
                sel_factor = rescaling_array[ind_factor-ind_range:]
                sel_stdev = computed_std[ind_factor-ind_range:]
            else:
                sel_factor = rescaling_array[ind_factor-ind_range:ind_factor+ind_range]
                sel_stdev = computed_std[ind_factor-ind_range:ind_factor+ind_range]

            coeff = np.polyfit(sel_factor, sel_stdev, 2)

            telluric_factor = - coeff[1] / (2*coeff[0])

            print('   telluric factor: {0:7f}'.format(telluric_factor))
            print()

            processed['template']['telluric_factor'] = telluric_factor
            processed['template']['rescaling_factor'] = sel_factor
            processed['template']['computed_std'] = computed_std
            processed['template']['polyfit'] = {
                'package': 'numpy', # in case we forget what we used...
                'order': 2,
                'coeff': coeff,
                'sel_factor': sel_factor,
                'sel_stdev': sel_stdev
             }

            """ 
            The template telluric spectrum is rebinned onto the 2D scale of the observations.
            Then it is rescaled according to the computed factor
            We assume that the wavelength solution did not change during the night
            """

            processed['template']['rebinned'] = {}
            processed['template']['rebinned']['flux'] = \
                rebin_1d_to_2d(telluric['template']['input']['wave'],
                               telluric['template']['input']['step'],
                               telluric['template']['input']['flux'],
                               telluric['rebinned']['wave'],
                               telluric['rebinned']['step'],
                               preserve_flux=False)

            processed['template']['rebinned']['ferr'] = \
                rebin_1d_to_2d(telluric['template']['input']['wave'],
                               telluric['template']['input']['step'],
                               telluric['template']['input']['ferr'],
                               telluric['rebinned']['wave'],
                               telluric['rebinned']['step'],
                               preserve_flux=False,
                               is_error=True)

            sel_out_of_range = ~((telluric['rebinned']['wave'] > telluric['template']['input']['range'][0]+1.) \
                                & (telluric['rebinned']['wave'] < telluric['template']['input']['range'][1]-1.))
            processed['template']['rebinned']['flux'][sel_out_of_range] = 1.
            processed['template']['rebinned']['ferr'][sel_out_of_range] = 0.1

            processed['telluric']['spectrum_noairmass'] = \
                (processed['template']['rebinned']['flux'] - 1.) * telluric_factor + 1.0

            processed['telluric']['spectrum_noairmass_err'] = processed['template']['rebinned']['ferr'] * telluric_factor

        except MissingFileException:
            print(' *** Missing telluric_template keyword in configuration file ***')
            print()

        except OutOfRangeException:
            print(' *** Wavelength range for the calculation of the rescaling factor is outside the boundaries ***')
            print('     Rescaling factor wavelength range: {0:7.2f} to {1:7.2f}'.format(
                template_dict['fit_range'][0], template_dict['fit_range'][1]))
            print('     Shared data wavelength range     : {0:7.2f} to {1:7.2f}'.format(
                shared_data['coadd']['wavelength_range'][0],
                shared_data['coadd']['wavelength_range'][1]))

            print()

        """
        Introduction of a second telluric correction, where the telluric spectrum depends linearly on the precipitable 
        water vapour (PWV). As such, the values stored in the files are the coefficient of the PWV term, while the 
        baseline is given by the outcome of the previous step, i.e., the spectrum_noairmass array
        """
        try:

            """ 
            The algorith is essentially the same as in the previous step, with the exception of the calculation
            of the telluric spectrum at each iteration of the chi-square minimization
            """

            if 'telluric_linear_term' not in instrument_dict[instrument_name]:
                raise MissingKeywordException()

            linear_dict = instrument_dict[instrument_name]['telluric_linear_term']

            if linear_dict['fit_range'][0] > shared_data['coadd']['wavelength_range'][1] or \
                    linear_dict['fit_range'][1] < shared_data['coadd']['wavelength_range'][0]:
                raise OutOfRangeException()

            print('    PWV-dependent telluric spectrum')
            print('    instrument :', instrument_name)
            print('    linear     :', linear_dict['file'])
            print('    fit_range  :', linear_dict['fit_range'])
            print()

            """ Retrieving the linear data"""
            telluric_linear_data = np.genfromtxt(linear_dict['file'])

            telluric['linear']['input'] = {
                'range': [np.amin(telluric_linear_data[:, 0]), np.amax(telluric_linear_data[-1, 0])],
                'wave': telluric_linear_data[:, 0],
                'coef': telluric_linear_data[:, 1],
                'cerr': telluric_linear_data[:, 2],
                'step': telluric_linear_data[:, 3]
            }

            processed['linear'] = {
                'array_wave': [],
                'array_data': {},
                'array_move': {}
            }

            """ Again we assume that the wavelength solution did not change during the night  """
            tel_selection = (input_data[obs_reference]['wave'] > linear_dict['fit_range'][0]) & \
                            (input_data[obs_reference]['wave'] < linear_dict['fit_range'][1])
            processed['linear']['wave_sel'] = input_data[obs_reference]['wave'][tel_selection]
            processed['linear']['step_sel'] = input_data[obs_reference]['step'][tel_selection]

            for obs in lists['telluric']:
                processed['linear'][obs] = processed[obs]['e2ds_rescaled'][tel_selection]
                """ saving the wavelength array for plotting purpose"""
                processed['linear']['array_wave'].extend(processed['linear']['wave_sel'])

            rescaling_array = 10**np.arange(-1, np.log10(50), 0.1)
            computed_std = np.zeros(len(rescaling_array))

            for n_rescaling_factor, rescaling_factor in enumerate(rescaling_array):

                """ 
                Template spectrum is rebinned onto the observations wavelength scale, using only the wavelength 
                range selected for the computation of the rescaling factor 
                """
                linear_rebinned_flux = \
                    rebin_1d_to_1d(telluric['linear']['input']['wave'],
                                   telluric['linear']['input']['step'],
                                   telluric['linear']['input']['coef'],
                                   processed['linear']['wave_sel'],
                                   processed['linear']['step_sel'],
                                   preserve_flux=False)

                """ Saving the outcome to dictionary """
                linear_rebinned_flux *= rescaling_factor
                linear_rebinned_flux += processed['telluric']['spectrum_noairmass'][tel_selection]

                e2ds_corrected = []

                for obs in lists['telluric']:

                    e2ds_corrected.extend(processed['linear'][obs] /
                                            np.power(linear_rebinned_flux, input_data[obs]['AIRMASS']))

                computed_std[n_rescaling_factor] = np.nanstd(e2ds_corrected)
                label_dict = '{0}'.format(rescaling_factor)
                processed['linear']['array_data'][label_dict] = e2ds_corrected
                processed['linear']['array_move'][label_dict] = rescaling_factor


                #print(n_rescaling_factor, rescaling_factor, computed_std[n_rescaling_factor])
                #plt.scatter(processed['linear']['array_wave'], e2ds_corrected)
                #plt.plot(processed['linear']['wave_sel'], linear_rebinned_flux)
                #plt.show()


            #plt.scatter(rescaling_array, computed_std)
            #plt.show()

            """ Selection of the PWV value with the lowest scatter """
            ind_factor = np.argmin(computed_std)
            ind_range = 3
            if ind_factor < ind_range:
                sel_factor = rescaling_array[0:ind_factor+ind_range]
                sel_stdev = computed_std[0:ind_factor+ind_range]
            elif ind_factor > len(rescaling_array) - ind_range:
                sel_factor = rescaling_array[ind_factor-ind_range:]
                sel_stdev = computed_std[ind_factor-ind_range:]
            else:
                sel_factor = rescaling_array[ind_factor-ind_range:ind_factor+ind_range]
                sel_stdev = computed_std[ind_factor-ind_range:ind_factor+ind_range]

            coeff = np.polyfit(sel_factor, sel_stdev, 2)

            PWV_value = - coeff[1] / (2*coeff[0])

            print('   PWV value      : {0:7f}'.format(PWV_value))
            print()

            processed['linear']['PWV_value'] = PWV_value
            processed['linear']['PWV_closest'] = sel_factor
            processed['linear']['computed_std'] = computed_std
            processed['linear']['polyfit'] = {
                'package': 'numpy', # in case we forget what we used...
                'order': 2,
                'coeff': coeff,
                'sel_factor': sel_factor,
                'sel_stdev': sel_stdev
             }

            """ 
            The linear coefficient for the PWV-dependent part are rebinned onto the 2D  scale of the observations. 
            Then the teluric spectrum is computed, using also the baseline from the previous step 
            """

            processed['linear']['rebinned'] = {}
            processed['linear']['rebinned']['coef'] = \
                rebin_1d_to_2d(telluric['linear']['input']['wave'],
                               telluric['linear']['input']['step'],
                               telluric['linear']['input']['coef'],
                               telluric['rebinned']['wave'],
                               telluric['rebinned']['step'],
                               preserve_flux=False)

            processed['linear']['rebinned']['cerr'] = \
                rebin_1d_to_2d(telluric['linear']['input']['wave'],
                               telluric['linear']['input']['step'],
                               telluric['linear']['input']['cerr'],
                               telluric['rebinned']['wave'],
                               telluric['rebinned']['step'],
                               preserve_flux=False,
                               is_error=True)

            sel_out_of_range = ~((telluric['rebinned']['wave'] > telluric['linear']['input']['range'][0]+1.) \
                                & (telluric['rebinned']['wave'] < telluric['linear']['input']['range'][1]-1.))
            processed['linear']['rebinned']['coef'][sel_out_of_range] = 0.
            processed['linear']['rebinned']['cerr'][sel_out_of_range] = 0.1

            processed['telluric']['spectrum_noairmass'] += processed['linear']['rebinned']['coef'] * PWV_value
            processed['telluric']['spectrum_noairmass_err'] = np.sqrt(
                processed['telluric']['spectrum_noairmass_err']**2 +
                (processed['linear']['rebinned']['cerr'] * PWV_value)**2)


        except MissingFileException:
            print(' *** Missing telluric_linear_coeff keyword in configuration file ***')
            print()

        except OutOfRangeException:
            print(' *** Wavelength range for the calculation of the PWV-dependent telluric spectrum is outside the boundaries ***')
            print()

        for obs in lists['observations']:
            """ Correction of telluric lines for the average airmass value, following Wyttenbach et al. 2015 """
            processed[obs]['e2ds_corrected'] = processed[obs]['e2ds_rescaled'] / \
                                                    np.power(processed['telluric']['spectrum_noairmass'],
                                                             input_data[obs]['AIRMASS'])
            processed[obs]['e2ds_corrected_err'] = processed[obs]['e2ds_rescaled_err'] / \
                                                    np.power(processed['telluric']['spectrum_noairmass'],
                                                             input_data[obs]['AIRMASS'])

        for obs in lists['observations']:
            # Correction of telluric lines

            telluric[obs] = {}

            telluric[obs]['spectrum_noairmass'] = processed['telluric']['spectrum_noairmass']

            telluric[obs]['airmass'] = input_data[obs]['AIRMASS']
            telluric[obs]['airmass_ref'] = processed['airmass_ref']

            """ Set anomalosly low point to one (e.g. when the template is not computed)"""
            telluric[obs]['null'] = telluric[obs]['spectrum_noairmass'] < 0.001
            telluric[obs]['spectrum_noairmass'][telluric[obs]['null']] = 1.0

            telluric[obs]['spectrum'] = np.power(processed['telluric']['spectrum_noairmass'],
                                       input_data[obs]['AIRMASS'] - processed['airmass_ref'])

            telluric[obs]['spline_noairmass'] = telluric[obs]['spectrum_noairmass'].copy()

            """ No need to compute the spline approximation since we are already dealing with a very high SNR template"""
            telluric[obs]['spline'] = np.power(telluric[obs]['spline_noairmass'],
                                       input_data[obs]['AIRMASS'] - processed['airmass_ref'])

            """ copy the keyword for future use"""
            telluric[obs]['airmass'] = input_data[obs]['AIRMASS']

            telluric[obs]['telluric_corrected'] = processed[obs]['e2ds_corrected']
            telluric[obs]['telluric_corrected_err'] = processed[obs]['e2ds_corrected_err']

            save_to_cpickle('telluric', telluric, config_in['output'], night)
            save_to_cpickle('telluric_processed', processed, config_in['output'], night)

        print()
        print("Night ", night, " completed")


def plot_telluric_template_alternative(config_in, night_input=''):
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