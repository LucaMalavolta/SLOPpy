from __future__ import print_function, division
from SLOPpy.subroutines.common import *
from SLOPpy.subroutines.spectral_subroutines import *
from SLOPpy.subroutines.io_subroutines import *
from SLOPpy.subroutines.fit_subroutines import *
from SLOPpy.subroutines.shortcuts import *

__all__ = ["compute_telluric_airmass_observerRF_chunks"]


def compute_telluric_airmass_observerRF_chunks(config_in):

    compute_telluric_observerRF_chunks(config_in,
                                       n_iterations=1,
                                       use_berv=False,
                                       use_reference_airmass=False,
                                       subroutine_name='compute_telluric_airmass_observerRF_chunks')


def compute_telluric_observerRF_chunks(config_in, **kwargs):

    night_dict = from_config_get_nights(config_in)

    for night in night_dict:

        print()
        print("compute_telluric_airmass_observerRF_chunks         Night: ", night)

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
        input_data = retrieve_observations(config_in['output'], night, lists['observations'])
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

            processed[obs] = {}

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

            """ for plotting purpose only"""
            processed[obs]['wave'] = input_data[obs]['wave']
            processed[obs]['e2ds'] = input_data[obs]['e2ds']
            processed[obs]['e2ds_err'] = input_data[obs]['e2ds_err']

        for niter in xrange(0, kwargs['n_iterations']):
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

            for order in xrange(0, processed['n_orders']):
                print("                                 - order ", repr(order))

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
                    if observational_pams['linear_fit_method']== 'linear_curve_fit':

                        abs_slope[order, :], zero_point[order, :] = \
                            airmass_linear_curve_fit(airmass, logi_array, sigi_array, processed['n_pixels'])

                    else:

                        abs_slope[order, :], zero_point[order, :] = \
                            airmass_linear_lstsq(airmass, logi_array)

                """ Saving the outcome to dictionary """
                processed['telluric']['order_'+repr(order)] = {'logi_array': logi_array, 'sigi_array': sigi_array}

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

            for order in xrange(0, processed['n_orders']):
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
