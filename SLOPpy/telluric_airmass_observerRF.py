from __future__ import print_function, division
from SLOPpy.subroutines.common import *
from SLOPpy.subroutines.spectral_subroutines import *
from SLOPpy.subroutines.io_subroutines import *
from SLOPpy.subroutines.fit_subroutines import *
from SLOPpy.subroutines.plot_subroutines import *
from SLOPpy.subroutines.shortcuts import *

__all__ = ["compute_telluric_airmass_berv_observerRF",
           "plot_telluric_airmass_berv_observerRF",
           "compute_telluric_airmass_observerRF",
           "plot_telluric_airmass_observerRF",
           "compute_telluric_airmass_reference_observerRF",
           "plot_telluric_airmass_reference_observerRF",
           "compute_telluric_airmass_berv_reference_observerRF",
           "plot_telluric_airmass_berv_reference_observerRF"
           ]


def compute_telluric_airmass_berv_observerRF(config_in):

    compute_telluric_observerRF(config_in,
                                n_iterations=1,
                                use_berv=True,
                                use_reference_airmass=False,
                                subroutine_name='telluric_airmass_berv_observerRF')


def compute_telluric_airmass_observerRF(config_in):

    compute_telluric_observerRF(config_in,
                                n_iterations=1,
                                use_berv=False,
                                use_reference_airmass=False,
                                subroutine_name='telluric_airmass_observerRF')


def compute_telluric_airmass_reference_observerRF(config_in):

    compute_telluric_observerRF(config_in,
                                n_iterations=1,
                                use_berv=False,
                                use_reference_airmass=True,
                                subroutine_name='telluric_airmass_reference_observerRF')

def compute_telluric_airmass_berv_reference_observerRF(config_in):
    compute_telluric_observerRF(config_in,
                                n_iterations=1,
                                use_berv=True,
                                use_reference_airmass=True,
                                subroutine_name='telluric_airmass_berv_reference_observerRF')


def plot_telluric_airmass_berv_observerRF(config_in, night_input):
    """ Alias to simplify the configuration file"""
    plot_telluric_airmass_observerRF(config_in, night_input)


def plot_telluric_airmass_reference_observerRF(config_in, night_input):
    """ Alias to simplify the configuration file"""
    plot_telluric_airmass_observerRF(config_in, night_input)


def plot_telluric_airmass_berv_reference_observerRF(config_in, night_input):
    """ Alias to simplify the configuration file"""
    plot_telluric_airmass_observerRF(config_in, night_input)


def compute_telluric_observerRF(config_in, **kwargs):

    night_dict = from_config_get_nights(config_in)

    for night in night_dict:

        print()
        print("compute_telluric_airmass_observerRF        Night: ", night)

        try:
            telluric = load_from_cpickle('telluric', config_in['output'], night)
            continue
        except:
            print("No telluric correction file found, computing now ")
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

        # There must be a more elegant way to do this, but I'm, not aware of it
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

            if processed['n_orders'] == 0:
                processed['n_orders'] = input_data[obs]['orders']
                processed['n_pixels'] = input_data[obs]['wave_size']

        """ Reference airmass for iterative correction of airmass"""
        if kwargs['use_reference_airmass']:
            airmass_temp = np.zeros(lists['n_transit_in'])
            for n_obs, obs in enumerate(lists['transit_in']):
                # This is to ensure that airmass, berv and rvc are associated to the correct spectra
                airmass_temp[n_obs] = input_data[obs]['AIRMASS']
            processed['airmass_ref'] = np.average(airmass_temp)
        else:
            processed['airmass_ref'] = 0.000

        for obs in lists['observations']:
            processed[obs]['e2ds_precorrected'] = processed[obs]['e2ds_rescaled'][:]
            processed[obs]['e2ds_precorrected_err'] = input_data[obs]['e2ds_err'] / processed[obs]['e2ds_rescaling']

        for niter in range(0, kwargs['n_iterations']):
            if kwargs['n_iterations'] > 1:
                print("NITER:            ", niter)

            for obs in lists['telluric']:
                processed[obs]['logI'] = np.log(processed[obs]['e2ds_precorrected'])
                processed[obs]['logI_err'] = processed[obs]['e2ds_precorrected_err']/processed[obs]['e2ds_precorrected']

            processed['telluric'] = {}

            abs_slope = np.ones([processed['n_orders'], processed['n_pixels']], dtype=np.double)
            line_shift = np.ones([processed['n_orders'], processed['n_pixels']], dtype=np.double)
            zero_point = np.ones([processed['n_orders'], processed['n_pixels']], dtype=np.double)
            pearson_r = np.zeros([processed['n_orders'], processed['n_pixels']], dtype=np.double)
            pearson_p = np.zeros([processed['n_orders'], processed['n_pixels']], dtype=np.double)

            airmass = np.zeros(lists['n_tellurics'], dtype=np.double)
            berv = np.zeros(lists['n_tellurics'], dtype=np.double)
            rvc = np.zeros(lists['n_tellurics'], dtype=np.double)

            for n_obs, obs in enumerate(lists['telluric']):
                # This is to ensure that airmass, berv and rvc are associated to the correct spectra
                processed['telluric'][obs] = {'n_obs': n_obs}
                airmass[n_obs] = input_data[obs]['AIRMASS']
                berv[n_obs] = input_data[obs]['BERV']
                rvc[n_obs] = input_data[obs]['RVC']

            for order in range(0, processed['n_orders']):

                logi_array = np.empty([lists['n_tellurics'], processed['n_pixels']], dtype=np.double)
                sigi_array = np.empty([lists['n_tellurics'], processed['n_pixels']], dtype=np.double)

                for obs in lists['telluric']:
                    n_obs = processed['telluric'][obs]['n_obs']
                    logi_array[n_obs, :] = processed[obs]['logI'][order, :]
                    sigi_array[n_obs, :] = processed[obs]['logI_err'][order, :]

                """ The user has the option to select between different approaches to
                    extract the telluric absorption spectrum
                    To-Do: move this section to a subroutine for cythonization"""

                if kwargs['use_berv']:
                    if observational_pams['linear_fit_method'] == 'linear_curve_fit':

                        abs_slope[order, :], line_shift[order, :], zero_point[order, :] = \
                            berv_linear_curve_fit_modified(airmass, berv, logi_array, sigi_array, processed['n_pixels'])

                    else:
                        abs_slope[order, :], line_shift[order, :], zero_point[order, :] = \
                            berv_linear_lstsq(airmass, berv, logi_array)

                else:
                    if observational_pams['linear_fit_method'] == 'linear_curve_fit':

                        abs_slope[order, :], zero_point[order, :] = \
                            airmass_linear_curve_fit(airmass, logi_array, sigi_array, processed['n_pixels'])

                        #abs_slope[order, :], zero_point[order, :] = \
                        #    airmass_linear_curve_fit_ransac(airmass, logi_array, sigi_array, processed['n_pixels'])

                        #obs_ref = lists['observations'][0]
                        #plt.plot(processed[obs_ref]['wave'][order,:], processed[obs_ref]['e2ds_rescaled'][order,:])
                        #for iii in range(0, processed['n_pixels']):
                        #    if iii < 3700 or iii > 3720: continue
                        #    plt.axvline(processed[obs_ref]['wave'][order,iii])
                        #plt.show()
                        #
                        #ik=0
                        #air_arr = np.arange(1.2, 2.5, 0.1)
                        #for iii in range(0, processed['n_pixels']):
                        #
                        #    if iii < 3700 or iii > 3720: continue
                        #    plt.errorbar(airmass, logi_array[:, iii]+ik, yerr=sigi_array[:, iii], fmt='o')
                        #    print(np.exp(abs_slope[order, iii]))
                        #    plt.plot(air_arr, air_arr*abs_slope[order, iii] + zero_point[order, iii]+ik)
                        #    ik -= 0.20
                        #plt.show()
                    else:

                        abs_slope[order, :], zero_point[order, :] = \
                            airmass_linear_lstsq(airmass, logi_array)

                plt.show()
                """ Saving the outcome to dictionary """
                processed['telluric']['order_'+repr(order)] = {'logi_array': logi_array, 'sigi_array': sigi_array}

            if kwargs.get('use_template', False):
                telluric_template_data = np.genfromtxt(night_dict[night]['telluric_template'])

                spectrum_noairmass = np.exp(abs_slope)

                obs_reference = lists['observations'][0]
                telluric['template'] = {
                    'input':{
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

                plt.plot(telluric['template']['input']['wave'], telluric['template']['input']['flux'], zorder=1, c='C0')
                plt.scatter(telluric['template']['rebinned']['wave'], telluric['template']['rebinned']['flux'],zorder=2, s=2)
                plt.scatter(telluric['template']['rebinned']['wave'],spectrum_noairmass, alpha=0.5, s=1, zorder=3)

                factor_list = []
                slope_list = []

                for order in range(0, processed['n_orders']):
                    fit_selection = (telluric['template']['rebinned']['flux'][order, :] < 1.0)
                    # Check if there are telluric lines in this wavelength range
                    if np.sum(fit_selection) > 30:
                        #telluric_factor, telluric_slope, success_flag = find_telluric_rescaling_factor(
                        #    spectrum_noairmass[order, :],
                        #    telluric['template']['rebinned']['flux'][order, :]
                        #)
                        telluric_factor, telluric_slope, success_flag = find_telluric_rescaling_factor_2steps(
                            spectrum_noairmass[order, :],
                            telluric['template']['rebinned']['flux'][order, :]
                        )

                        if success_flag:
                            factor_list.extend([telluric_factor])
                            slope_list.extend([telluric_slope])

                if len(factor_list)>1:
                    telluric_slope = np.median(slope_list)
                    telluric_factor = np.median(factor_list)
                elif len(factor_list) == 1:
                    telluric_slope = slope_list[0]
                    telluric_factor = factor_list[0]
                else:
                    telluric_slope = 0.00
                    telluric_factor = 0.00

                print('   telluric factor: {0:7f} (correction slope: {1:7f}'.format(telluric_factor,telluric_slope))
                print()
                #print(telluric_factor, success_flag)
                #quit()

                processed['telluric']['spectrum_noairmass'] = \
                    (telluric['template']['rebinned']['flux']-1.)*telluric_factor + 1.0

                plt.plot(telluric['template']['input']['wave'],
                         (telluric['template']['input']['flux']-1.)*telluric_factor + 1.0, zorder=1, c='C1')

                plt.show()


            else:
                processed['telluric']['spectrum_noairmass'] = np.exp(abs_slope)

            for obs in lists['observations']:
                """ Correction of telluric lines for the average airmass value, following Wyttenbach et al. 2015 """
                processed[obs]['e2ds_corrected'] = processed[obs]['e2ds_precorrected'] / \
                                                        np.power(processed['telluric']['spectrum_noairmass'],
                                                                 input_data[obs]['AIRMASS'] -
                                                                 processed['airmass_ref'])
                processed[obs]['e2ds_corrected_err'] = processed[obs]['e2ds_precorrected_err'] / \
                                                        np.power(processed['telluric']['spectrum_noairmass'],
                                                                 input_data[obs]['AIRMASS'] -
                                                                 processed['airmass_ref'])

        for obs in lists['observations']:
            # Correction of telluric lines

            telluric[obs] = {}

            telluric[obs]['spectrum_noairmass'] = processed['telluric']['spectrum_noairmass']

            telluric[obs]['airmass'] = input_data[obs]['AIRMASS']
            telluric[obs]['airmass_ref'] = processed['airmass_ref']
            telluric[obs]['null'] = telluric[obs]['spectrum_noairmass'] < 0.001
            telluric[obs]['spectrum_noairmass'][ telluric[obs]['null']] = 1.0

            telluric[obs]['spectrum'] = np.power(processed['telluric']['spectrum_noairmass'],
                                       input_data[obs]['AIRMASS'] - processed['airmass_ref'])

            telluric[obs]['spline_noairmass'] = np.ones([input_data[obs]['n_orders'],
                                                         input_data[obs]['n_pixels']],
                                                        dtype=np.double)

            for order in range(0, processed['n_orders']):
                telluric[obs]['spline_noairmass'][order, :], _, _ = \
                    compute_spline(input_data[obs]['wave'][order, :],
                                   telluric[obs]['spectrum_noairmass'][order, :],
                                   0.05)

            telluric[obs]['spline'] = np.power(telluric[obs]['spline_noairmass'],
                                       input_data[obs]['AIRMASS'] - processed['airmass_ref'])

            telluric[obs]['telluric_corrected'] = processed[obs]['e2ds_corrected']
            telluric[obs]['telluric_corrected_err'] = processed[obs]['e2ds_corrected_err']

        save_to_cpickle('telluric', telluric, config_in['output'], night)
        save_to_cpickle('telluric_processed', processed, config_in['output'], night)

        print()
        print("Night ", night, " completed")


def plot_telluric_airmass_observerRF(config_in, night_input=''):
    import matplotlib.pyplot as plt

    night_dict = from_config_get_nights(config_in)
    instrument_dict = from_config_get_instrument(config_in)
    system_dict = from_config_get_system(config_in)

    if night_input == '':
        night_list = night_dict
    else:
        night_list = np.atleast_1d(night_input)

    for night in night_list:

        print("plot_telluric_airmass_observerRF           Night: ", night)

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