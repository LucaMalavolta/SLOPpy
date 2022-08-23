from __future__ import print_function, division
from SLOPpy.subroutines.common import *
from SLOPpy.subroutines.spectral_subroutines import *
from SLOPpy.subroutines.io_subroutines import *
from SLOPpy.subroutines.fit_subroutines import *
from SLOPpy.subroutines.plot_subroutines import *
from SLOPpy.subroutines.shortcuts import *
from SLOPpy.subroutines.eso_skycalc_cli import get_eso_sckycalc_harps

__all__ = ["compute_telluric_observerRF_skycalc", "plot_telluric_observerRF_skycalc"]


def compute_telluric_observerRF_skycalc(config_in):

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
            'subroutine': 'telluric_observerRF_skycalc',
            'n_orders': 0,
            'n_pixels': 0
        }

        telluric = {
            'subroutine': 'telluric_observerRF_skycalc',
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

        telluric['sky_template'] = {}
        instrument = night_dict[night]['instrument']

        for n_obs, obs in enumerate(lists['transit_in']):
            if n_obs >= lists['n_transit_in']/2.:
                obs_ref = obs
                break

        telluric['sky_template']['ref_observation'] = obs_ref
        if instrument == 'HARPS':
            telluric['sky_template']['use_eso_skycalc'] = True
            telluric['sky_template']['ref_airmass'] = input_data[obs_ref]['AIRMASS']
            telluric['sky_template']['data'] = np.ones([processed['n_orders'], processed['n_pixels']], dtype=np.double)

        else:
            telluric_model_file = config_in['instruments'][instrument]['telluric_model']
            telluric['sky_template']['use_eso_skycalc'] = False
            telluric['sky_template']['ref_airmass'] = 1.00000
            telluric['sky_template']['data'] = fits.open(telluric_model_file)

        for obs in lists['observations']:
            processed[obs]['e2ds_precorrected'] = processed[obs]['e2ds_rescaled'][:]
            processed[obs]['e2ds_precorrected_err'] = input_data[obs]['e2ds_err'] / processed[obs]['e2ds_rescaling']

        for niter in xrange(0, 1):
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

                if observational_pams['linear_fit_method']== 'linear_curve_fit':

                    abs_slope[order, :], line_shift[order, :], zero_point[order, :] = \
                        berv_linear_curve_fit_modified(airmass, berv, logi_array, sigi_array, processed['n_pixels'])

                else:
                    abs_slope[order, :], line_shift[order, :], zero_point[order, :] = \
                        berv_linear_lstsq(airmass, berv, logi_array)

                """ Saving the outcome to dictionary """
                processed['telluric']['order_'+repr(order)] = {'logi_array': logi_array, 'sigi_array': sigi_array}

                if telluric['sky_template']['use_eso_skycalc']:
                    wave_model, step_model, tran_model, terr_model = \
                        get_eso_sckycalc_harps(obs_ref,
                                               [input_data[obs_ref]['wave'][0], input_data[obs_ref]['wave'][-1]],
                                               input_data[obs_ref]['RA'],
                                               input_data[obs_ref]['DEC'],
                                               night, config_in['output'])

                    tran_model_rebinned = \
                        rebin_1d_to_1d(wave_model,
                                       step_model,
                                       tran_model,
                                       input_data[obs_ref]['wave'],
                                       input_data[obs_ref]['step'],
                                       preserve_flux=False)

                    telluric['sky_template']['data'][order, :] = \
                        np.power(tran_model_rebinned,
                                 1./telluric['sky_template']['ref_airmass'])

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

            telluric[obs]['airmass'] = input_data[obs]['AIRMASS']
            telluric[obs]['airmass_ref'] = processed['airmass_ref']
            telluric[obs]['null'] = telluric[obs]['spectrum_noairmass'] < 0.001
            telluric[obs]['spectrum_noairmass'][ telluric[obs]['null']] = 1.0

            telluric[obs]['telluric_corrected'] = processed[obs]['e2ds_corrected']
            telluric[obs]['telluric_corrected_err'] = processed[obs]['e2ds_corrected_err']

        save_to_cpickle('telluric', telluric, config_in['output'], night)
        save_to_cpickle('telluric_processed', processed, config_in['output'], night)

        print()
        print("Night ", night, " completed")


def plot_telluric_observerRF_skycalc(config_in, night_input=''):
    import matplotlib.pyplot as plt

    night_dict = from_config_get_nights(config_in)
    instrument_dict = from_config_get_instrument(config_in)
    system_dict = from_config_get_system(config_in)

    if night_input=='':
        night_list = night_dict
    else:
        night_list = np.atleast_1d(night_input)

    for night in night_list:

        print("plot_telluric_airmass_stellarRF            Night: ", night)

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

            _, e2ds_rescaled , _ = \
                perform_rescaling(processed[obs]['wave'],
                                  processed[obs]['e2ds'],
                                  processed[obs]['e2ds_err'],
                                  observational_pams['wavelength_rescaling'])

            e2ds_rescaled_corrected_spectrum = e2ds_rescaled / telluric[obs]['spectrum']

            for order in range(0, processed[obs]['n_orders']):

                ax1.plot(processed[obs]['wave'][order, :],
                            e2ds_rescaled[order, :],
                            c=line_colors[i], lw=1, alpha=0.5)
                ax1.scatter(processed[obs]['wave'][order, :],
                            e2ds_rescaled_corrected_spectrum[order, :],
                            s=1, c=line_colors[i])

                ax2.plot(processed[obs]['wave'][order, :],
                         telluric[obs]['spectrum'][order, :],
                         c=line_colors[i])
                ax2.axhline(1.00, c='k')


        ax1.legend(loc=3)
        ax1.set_title('Night: ' + night)

        ax2.set_xlabel('$\lambda$ [$\AA$]')

        obs_ref = telluric['sky_template']['ref_observation']

        for order in range(0, processed[obs]['n_orders']):
            ax2.scatter(processed[obs_ref]['wave'][order, :], telluric[obs_ref]['spectrum_noairmass'][order, :], s=2, c='C0', zorder=1000)
            ax2.plot(processed[obs_ref]['wave'][order, :], telluric['sky_template']['data'][order, :], c='C0', zorder=1000)
            ax1.plot(processed[obs_ref]['wave'][order, :], telluric['sky_template']['data'][order, :], c='C0', zorder=1000)

        sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=colors[0], vmax=colors[-1]))
        sm.set_array([])  # You have to set a dummy-array for this to work...
        cbar = plt.colorbar(sm, cax=cbax1)
        cbar.set_label('BJD - 2450000.0')
        fig.subplots_adjust(wspace=0.05, hspace=0.4)
        plt.show()