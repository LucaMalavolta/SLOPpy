from __future__ import print_function, division
from SLOPpy.subroutines.common import *
from SLOPpy.subroutines.spectral_subroutines import *
from SLOPpy.subroutines.io_subroutines import *
from SLOPpy.subroutines.fit_subroutines import *
from SLOPpy.subroutines.shortcuts import *
from SLOPpy.subroutines.rebin_subroutines import *
from SLOPpy.subroutines.clv_rm_subroutines import *

__all__ = ['write_output_spectra']

def write_output_spectra(config_in):

    subroutine_name = 'write_output_spectra'
    clv_rm_correction = True

    night_dict = from_config_get_nights(config_in)
    instrument_dict = from_config_get_instrument(config_in)
    system_dict = from_config_get_system(config_in)
    planet_dict = from_config_get_planet(config_in)

    shared_data = load_from_cpickle('shared', config_in['output'])

    lightcurve_dict = from_config_get_transmission_lightcurve(config_in)

    for night in night_dict:

        clv_rm_correction = True
        try:
            clv_rm_modelling = load_from_cpickle('clv_rm_modelling', config_in['output'], night)
        except:
            clv_rm_correction = False

        message = 'Computing'
        try:
            output_spectra = load_from_cpickle(subroutine_name, config_in['output'], night)

            if clv_rm_correction and not output_spectra['clv_rm_correction']:
                message = 'Updating with CLV-corrected spectra'
                raise ValueError()

            print("{0:45s} Night:{1:15s}   {2:s}".format(subroutine_name, night, 'Retrieved'))
            continue
        except:
            print("{0:45s} Night:{1:15s}   {2:s}".format(subroutine_name, night, message))
            print()

        """ Retrieving the list of observations"""
        lists = load_from_cpickle('lists', config_in['output'], night)

        calib_data = load_from_cpickle('calibration_fibA', config_in['output'], night)
        input_data = retrieve_observations( config_in['output'], night, lists['observations'])
        observational_pams = load_from_cpickle('observational_pams', config_in['output'], night)

        output_spectra = {
            'subroutine': subroutine_name,
            'clv_rm_correction': clv_rm_correction
        }
        """ Adding the C-bands arrays to the dictionary"""

        if clv_rm_correction:
            clv_rm_modelling = load_from_cpickle('clv_rm_modelling', config_in['output'], night)

        for n_obs, obs in enumerate( lists['observations']):

            output_spectra[obs] = {}

            output_spectra[obs]['BJD'] = input_data[obs]['BJD']

            output_spectra[obs]['SRF_rebinned'] = \
                rebin_2d_to_1d(input_data[obs]['wave'],
                               input_data[obs]['step'],
                               input_data[obs]['e2ds'],
                               calib_data['blaze'],
                               shared_data['coadd']['wave'],
                               shared_data['coadd']['step'],
                               rv_shift=observational_pams[obs]['rv_shift_ORF2SRF'])

            output_spectra[obs]['SRF_rebinned_err'] = \
                rebin_2d_to_1d(input_data[obs]['wave'],
                               input_data[obs]['step'],
                               input_data[obs]['e2ds_err'],
                               calib_data['blaze'],
                               shared_data['coadd']['wave'],
                               shared_data['coadd']['step'],
                               rv_shift=observational_pams[obs]['rv_shift_ORF2SRF'])

            output_spectra[obs]['SRF_rescaling'], \
            output_spectra[obs]['SRF_rescaled'], \
            output_spectra[obs]['SRF_rescaled_err'] = perform_rescaling(
                shared_data['coadd']['wave'], output_spectra[obs]['SRF_rebinned'], output_spectra[obs]['SRF_rebinned_err'],
                observational_pams['wavelength_rescaling'])

            if clv_rm_correction:

                rv_shift = 0.0 # we always stay in SRF
                correction, _ = clv_rm_correction_factor_computation(
                    clv_rm_modelling, shared_data['coadd']['wave'], shared_data['coadd']['step'], rv_shift, obs)

                output_spectra[obs]['SRF_clv_rm_correction'] = correction
                output_spectra[obs]['SRF_clv_rm_rebinned'] = output_spectra[obs]['SRF_rebinned'] / correction
                output_spectra[obs]['SRF_clv_rm_rebinned_err'] = output_spectra[obs]['SRF_rebinned_err'] / correction

                output_spectra[obs]['SRF_clv_rm_rescaled'] = output_spectra[obs]['SRF_rescaled'] / correction
                output_spectra[obs]['SRF_clv_rm_rescaled_err'] = output_spectra[obs]['SRF_rescaled_err'] / correction

            try:
                output_spectra[obs]['phase'] = \
                    (observational_pams[obs]['BJD'] - night_dict[night]['time_of_transit'][0])/planet_dict['period'][0]
            except:
                output_spectra[obs]['phase'] = \
                    (observational_pams[obs]['BJD'] - night_dict[night]['time_of_transit'])/planet_dict['period'][0]

        save_to_cpickle(subroutine_name, output_spectra, config_in['output'], night)

