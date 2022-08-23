from __future__ import print_function, division
from SLOPpy.subroutines.common import *
from SLOPpy.subroutines.spectral_subroutines import *
from SLOPpy.subroutines.fit_subroutines import *
from SLOPpy.subroutines.io_subroutines import *
from SLOPpy.subroutines.shortcuts import *
from SLOPpy.subroutines.plot_subroutines import *

__all__ = ["compute_differential_refraction_preparation"]

def compute_differential_refraction_preparation(config_in, append_name=None):

    if append_name:
        subroutine_name = 'differential_refraction_preparation_' + append_name
        filename = 'refraction_'+append_name
    else:
        subroutine_name = 'differential_refraction_preparation'
        filename = 'refraction'

    night_dict = from_config_get_nights(config_in)
    instrument_dict = from_config_get_instrument(config_in)

    shared_data = load_from_cpickle('shared', config_in['output'])
    reference_flux = np.empty([len(night_dict), shared_data['coadd']['size']], dtype=np.double)
    reference_wght = np.zeros([len(night_dict), shared_data['coadd']['size']], dtype=np.double)
    reference_mask = np.ones([len(night_dict), shared_data['coadd']['size']], dtype=bool)
    compute_reference = False

    """ make sure that all the spectra are computed in the same reference system is cross-calibrations is used"""
    absolute_SRF = False
    for n_night, night in enumerate(night_dict):
        if night_dict[night]['refraction'].get('reference_night', False) \
                or night_dict[night]['refraction'].get('reference_instrument', False):
            absolute_SRF = True

    for n_night, night in enumerate(night_dict):

        try:
            preparation = load_from_cpickle(filename +'_preparation', config_in['output'], night)
            print("{0:45s} Night:{1:15s}   {2:s}".format(subroutine_name, night, 'Retrieved'))
            continue
        except:
            print("{0:45s} Night:{1:15s}   {2:s}".format(subroutine_name, night, 'Computing'))
            print()

        instrument = night_dict[night]['instrument']

        """ Retrieving the list of observations"""
        lists = load_from_cpickle('lists', config_in['output'], night)

        calib_data = load_from_cpickle('calibration_fibA', config_in['output'], night)
        input_data = retrieve_observations(config_in['output'], night, lists['observations'],
                                           use_refraction=False)
        observational_pams = load_from_cpickle('observational_pams', config_in['output'], night)

        defined_reference = night_dict[night]['refraction'].get('reference', False)

        if defined_reference:
            preparation = {
                'subroutine': 'differential_refraction_preparation',
                'coadd': {
                    'wave': shared_data['coadd']['wave'],
                    'step': shared_data['coadd']['step'],
                    'size': shared_data['coadd']['size']
                },
                'absolute_SRF': True,
                'reference_coadd': True
            }
        else:
            preparation = {
                'subroutine': 'differential_refraction_preparation',
                'coadd': {
                    'wave': input_data['coadd']['wave'],
                    'step': input_data['coadd']['step'],
                    'size': input_data['coadd']['size'],
                },
                'absolute_SRF': absolute_SRF,
                'reference_coadd': False
            }

        total_flux = np.empty([len(lists['observations']), preparation['coadd']['size']], dtype=np.double)
        total_wght = np.zeros([len(lists['observations']), preparation['coadd']['size']], dtype=np.double)
        total_mask = np.ones([len(lists['observations']), preparation['coadd']['size']], dtype=bool)

        """ Rebinning of all the spectra """
        for n_obs, obs in enumerate(lists['observations']):

            print("  Spectral rebinning - Processing: ", obs)
            preparation[obs] = {}

            """ Rebinning of the spectra in the SRF, except for a fixed constant in order to minimize
                the difference between  
            """
            if preparation['absolute_SRF']:
                rv_shift = observational_pams[obs]['rv_shift_ORF2SRF']
            else:
                rv_shift = observational_pams[obs]['rv_shift_ORF2SRF_mod']

            preparation[obs]['flux_rebinned_stellarRF'] = \
                rebin_2d_to_1d(input_data[obs]['wave'],
                               input_data[obs]['step'],
                               input_data[obs]['e2ds'],
                               calib_data['blaze'],
                               preparation['coadd']['wave'],
                               preparation['coadd']['step'],
                               rv_shift=rv_shift)

            err_flux_rebinned_SRF = \
                rebin_2d_to_1d(input_data[obs]['wave'],
                               input_data[obs]['step'],
                               input_data[obs]['e2ds_err'],
                               calib_data['blaze'],
                               preparation['coadd']['wave'],
                               preparation['coadd']['step'],
                               rv_shift=rv_shift,
                               is_error=True)

            """ Zero or negative values are identified, flagged and substituted with another value """
            #processed[obs]['flux_rebinned_stellarRF'], \
            #processed[obs]['err_flux_rebinned_SRF'], \
            #processed[obs]['flux_rebinned_SRF_null'] = \
            #    replace_values_errors_with_interpolation_1d(processed[obs]['flux_rebinned_stellarRF'],
            #                                                processed[obs]['err_flux_rebinned_SRF'],
            #                                                force_positive=True)

            #processed[obs]['rescaling'], processed[obs]['rescaled'], processed[obs]['rescaled_err'] = \
            #    perform_rescaling(processed['coadd']['wave'],
            #                      processed[obs]['flux_rebinned_stellarRF'],
            #                      processed[obs]['err_flux_rebinned_SRF'],
            #                      observational_pams['wavelength_rescaling'])
            #

            """ Zero or negative values are identified, flagged and substituted with another value """
            preparation[obs]['flux_rebinned_stellarRF'], \
            err_flux_rebinned_SRF, \
            flux_rebinned_SRF_null = \
                replace_values_errors_with_interpolation_1d(preparation[obs]['flux_rebinned_stellarRF'],
                                                            err_flux_rebinned_SRF,
                                                            force_positive=True)

            preparation[obs]['rescaling'], preparation[obs]['rescaled'], preparation[obs]['rescaled_err'] = \
                perform_rescaling(preparation['coadd']['wave'],
                                  preparation[obs]['flux_rebinned_stellarRF'],
                                  err_flux_rebinned_SRF,
                                  observational_pams['wavelength_rescaling'])

            #if instrument_dict[instrument]['refraction'].get('reference_night', False):
            #    if night != instrument_dict[instrument]['refraction']['reference_night']:
            #        continue
            #
            #if instrument_dict[instrument]['refraction'].get('reference_instrument', False):
            #    if instrument != instrument_dict[instrument]['refraction']['reference_instrument']:
            #        continue

            if night_dict[night]['refraction'].get('use_all_observations', False) or obs in lists['telluric']:
                total_flux[n_obs, :] = preparation[obs]['rescaled']
                total_mask[n_obs, :] = flux_rebinned_SRF_null
                total_wght[n_obs, :] = 1. / (preparation[obs]['rescaled_err'] ** 2)
                print("      Observation added to reference spectrum")

        masked_array = np.ma.array(total_flux, mask=total_mask)
        rescaled_mask, sum_weights = np.ma.average(masked_array,
                                                   weights=total_wght,
                                                   axis=0,
                                                   returned=True)

        preparation['coadd']['rescaled'] = rescaled_mask.filled(0.00)
        sum_weights[sum_weights <= 0.0] = 1.0
        preparation['coadd']['rescaled_err'] = 1. / np.sqrt(sum_weights)

        preparation['coadd']['rescaled'], preparation['coadd']['rescaled_err'], preparation['coadd']['null'] = \
            replace_values_errors_with_interpolation_1d(preparation['coadd']['rescaled'],
                                                        preparation['coadd']['rescaled_err'],
                                                        force_positive=True)

        save_to_cpickle(filename + '_preparation', preparation, config_in['output'], night)

        if defined_reference == night \
                or defined_reference == instrument \
                or defined_reference == 'all' :

                compute_reference = True
                reference_flux[n_night, :] = preparation['coadd']['rescaled']
                reference_mask[n_night, :] = preparation['coadd']['null']
                reference_wght[n_night, :] = 1. / (preparation['coadd']['rescaled_err'] ** 2)

    if compute_reference:
        reference = {
            'wave': shared_data['coadd']['wave'],
            'step': shared_data['coadd']['step'],
            'size': shared_data['coadd']['size'],
        }

        masked_array = np.ma.array(reference_flux, mask=reference_mask)
        rescaled_mask, sum_weights = np.ma.average(masked_array,
                                                   weights=reference_wght,
                                                   axis=0,
                                                   returned=True)

        reference['rescaled'] = rescaled_mask.filled(0.00)
        sum_weights[sum_weights <= 0.0] = 1.0
        reference['rescaled_err'] = 1. / np.sqrt(sum_weights)

        reference['rescaled'], reference['rescaled_err'], reference['null'] = \
            replace_values_errors_with_interpolation_1d(reference['rescaled'],
                                                        reference['rescaled_err'],
                                                        force_positive=True)

        save_to_cpickle(filename + '_reference', reference, config_in['output'])


    print()

