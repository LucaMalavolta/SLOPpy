from __future__ import print_function, division
from SLOPpy.subroutines.common import *
from SLOPpy.subroutines.spectral_subroutines import *
from SLOPpy.subroutines.io_subroutines import *
from SLOPpy.subroutines.fit_subroutines import *
from SLOPpy.subroutines.plot_subroutines import *
from SLOPpy.subroutines.shortcuts import *

__all__ = ["compute_telluric_template",
           "plot_telluric_template",
           "compute_telluric_template_reference",
           "plot_telluric_template_reference"]

def compute_telluric_template(config_in):

    compute_telluric_template_routine(config_in,
                                n_iterations=1,
                                use_berv=True,
                                use_reference_airmass=False,
                                use_template=True,
                                subroutine_name='telluric_template')


def compute_telluric_template_reference(config_in):

    compute_telluric_template_routine(config_in,
                                n_iterations=1,
                                use_berv=True,
                                use_reference_airmass=True,
                                use_template=True,
                                subroutine_name='telluric_template')


def plot_telluric_template_reference(config_in, night_input=''):

    plot_telluric_template(config_in, night_input=night_input)


def compute_telluric_template_routine(config_in, **kwargs):
    """
    Lazy workaround
    :param config_in:
    :param kwargs:
    :return:
    """

    night_dict = from_config_get_nights(config_in)
    instrument_dict = from_config_get_instrument(config_in)

    for night in night_dict:

        instrument_name = night_dict[night]['instrument']
        template_dict = instrument_dict[instrument_name]['telluric_template']
        print()
        print("compute_telluric_template                  Night: ", night)
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
            'subroutine': kwargs['subroutine_name'],
            'n_orders': 0,
            'n_pixels': 0
        }

        telluric = {
            'subroutine': kwargs['subroutine_name'],
            'reference_frame': 'observer'
        }

        """ Retrieving the template data"""
        telluric_template_data = np.genfromtxt(template_dict['file'])

        obs_reference = lists['observations'][0]
        telluric['template'] = {
            'input': {
                'range': [np.amin(telluric_template_data[:, 0]), np.amax(telluric_template_data[-1, 0])],
                'wave': telluric_template_data[:, 0],
                'flux': telluric_template_data[:, 1],
                'ferr': telluric_template_data[:, 2],
                'step': telluric_template_data[:, 3]
            },
            'rebinned':{
                'wave': input_data[obs_reference]['wave'],
                'step': input_data[obs_reference]['step']
            }
        }

        """ Reference airmass for iterative correction of airmass"""
        if kwargs['use_reference_airmass']:
            airmass_temp = np.zeros(lists['n_transit_in'])
            for n_obs, obs in enumerate(lists['transit_in']):
                # This is to ensure that airmass, berv and rvc are associated to the correct spectra
                airmass_temp[n_obs] = input_data[obs]['AIRMASS']
            processed['airmass_ref'] = np.average(airmass_temp)
        else:
            processed['airmass_ref'] = 0.000
        processed['telluric'] = {}

        airmass = np.zeros(lists['n_observations'], dtype=np.double)
        berv = np.zeros(lists['n_observations'], dtype=np.double)
        rvc = np.zeros(lists['n_observations'], dtype=np.double)

        # There must be a more elegant way to do this, but I'm, not aware of it
        for n_obs, obs in enumerate(lists['observations']):

            tel_selection = (input_data[obs]['wave'] > template_dict['fit_range'][0]) & \
                            (input_data[obs]['wave'] < template_dict['fit_range'][1])

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

            processed[obs]['e2ds_sel'] = processed[obs]['e2ds_rescaled'][tel_selection]

            if processed['n_orders'] == 0:
                processed['wave_sel'] = input_data[obs]['wave'][tel_selection]
                processed['step_sel'] = input_data[obs]['step'][tel_selection]
                processed['n_orders'] = input_data[obs]['orders']
                processed['n_pixels'] = input_data[obs]['wave_size']

            # This is to ensure that airmass, berv and rvc are associated to the correct spectra
            processed['telluric'][obs] = {'n_obs': n_obs}
            airmass[n_obs] = input_data[obs]['AIRMASS']
            berv[n_obs] = input_data[obs]['BERV']
            rvc[n_obs] = input_data[obs]['RVC']

        processed['template_fit'] = {}

        rescaling_array = np.arange(0.05, 2.0, 0.05)
        computed_std = np.zeros(len(rescaling_array))

        processed['template_fit']['array_data'] = {}
        processed['template_fit']['array_move'] = {}
        processed['template_fit']['array_wave'] = []

        """ saving the wavelength array for plotting purpose"""
        for obs in lists['telluric']:
            processed['template_fit']['array_wave'].extend(processed['wave_sel'])

        for n_rescaling_factor, rescaling_factor in enumerate(rescaling_array):

            """ Template spectrum is rebinned onto the observations wavelength scale, using only the wavelength range
            selected for the computation of the rescaling factor
            """
            template_rebinned_flux = \
                rebin_1d_to_1d(telluric['template']['input']['wave'],
                               telluric['template']['input']['step'],
                               telluric['template']['input']['flux'],
                               processed['wave_sel'],
                               processed['step_sel'],
                               preserve_flux=False)

            """ Saving the outcome to dictionary """
            template_rebinned_flux -= 1.00
            template_rebinned_flux *= rescaling_factor
            template_rebinned_flux += 1.00

            e2ds_corrected = []

            ## DEBUG
            # wave_corrected = []
            # e2ds_original = []

            for obs in lists['telluric']:

                e2ds_corrected.extend(processed[obs]['e2ds_sel'] /
                                        np.power(template_rebinned_flux, input_data[obs]['AIRMASS']))

                # wave_corrected.extend(processed['wave_sel'])
                # e2ds_original.extend(processed[obs]['e2ds_sel'])

            computed_std[n_rescaling_factor] = np.std(e2ds_corrected)

            label_dict = '{0}'.format(rescaling_factor)
            # print(label_dict, computed_std[n_rescaling_factor])

            # processed['template_fit']['array_data'][label_dict] = e2ds_corrected
            # processed['template_fit']['array_move'][label_dict] = rescaling_factor
            # plt.scatter(wave_corrected, e2ds_original, s=2, alpha=0.5)
            # plt.scatter(wave_corrected, e2ds_corrected, s=2, alpha=0.5)
            # plt.plot(processed['wave_sel'],template_rebinned_flux)
            # plt.show()

        """ selection of the rescaling factor with the lowest scatter """
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

        processed['template_fit']['telluric_factor'] = telluric_factor
        processed['template_fit']['rescaling_factor'] = rescaling_factor
        processed['template_fit']['computed_std'] = computed_std
        processed['template_fit']['polyfit'] = {
            'package': 'numpy', # in case we forget what we used...
            'order': 2,
            'coeff': coeff,
            'sel_factor': sel_factor,
            'sel_stdev': sel_stdev
         }

        """ After being rescaled for the proper factor, the template telluric spectrum is rebinned onto the 2D
        scale of the observations """

        telluric['template']['rebinned']['flux'] = \
            rebin_1d_to_2d(telluric['template']['input']['wave'],
                           telluric['template']['input']['step'],
                           telluric['template']['input']['flux'],
                           telluric['template']['rebinned']['wave'],
                           telluric['template']['rebinned']['step'],
                           preserve_flux=False)

        telluric['template']['rebinned']['ferr'] = \
            rebin_1d_to_2d(telluric['template']['input']['wave'],
                           telluric['template']['input']['step'],
                           telluric['template']['input']['ferr'],
                           telluric['template']['rebinned']['wave'],
                           telluric['template']['rebinned']['step'],
                           preserve_flux=False,
                           is_error=True)

        sel_out_of_range = ~((telluric['template']['rebinned']['wave'] > telluric['template']['input']['range'][0]+1.) \
                            & (telluric['template']['rebinned']['wave'] < telluric['template']['input']['range'][1]-1.))
        telluric['template']['rebinned']['flux'][sel_out_of_range] = 1.
        telluric['template']['rebinned']['ferr'][sel_out_of_range] = 0.1


        processed['telluric']['spectrum_noairmass'] = \
            (telluric['template']['rebinned']['flux'] - 1.) * telluric_factor + 1.0

        telluric['airmass_ref'] = processed['airmass_ref']

        for obs in lists['observations']:
            """ Correction of telluric lines for the average airmass value, following Wyttenbach et al. 2015 """
            processed[obs]['e2ds_corrected'] = processed[obs]['e2ds_rescaled'] / \
                                                    np.power(processed['telluric']['spectrum_noairmass'],
                                                             input_data[obs]['AIRMASS'] - processed['airmass_ref'])
            processed[obs]['e2ds_corrected_err'] = processed[obs]['e2ds_rescaled_err'] / \
                                                    np.power(processed['telluric']['spectrum_noairmass'],
                                                             input_data[obs]['AIRMASS'] - processed['airmass_ref'])

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


def plot_telluric_template(config_in, night_input=''):
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
