import numpy as np
from SLOPpy.subroutines.rebin_subroutines import *


def clv_rm_correction_factor_computation(clv_rm_modelling,  wave, step, rv_shift, obs):
    ancillary = {}

    ancillary['norm_convolved_shifted'] = \
        rebin_1d_to_1d(clv_rm_modelling['common']['wave'],
                       clv_rm_modelling['common']['step'],
                       clv_rm_modelling['common']['norm_convolved'],
                       wave,
                       step,
                       rv_shift=rv_shift,
                       preserve_flux=False)

    ancillary['stellar_spectra_convolved_shifted'] = \
        rebin_1d_to_1d(clv_rm_modelling['common']['wave'],
                       clv_rm_modelling['common']['step'],
                       clv_rm_modelling[obs]['stellar_spectra_convolved'],
                       wave,
                       step,
                       rv_shift=rv_shift,
                       preserve_flux=False)

    wavelength_exclusion = \
        (wave <= clv_rm_modelling['common']['wave'][0] + 1) | \
        (wave >= clv_rm_modelling['common']['wave'][-1] - 1)

    wavelength_selection = \
        (wave > clv_rm_modelling['common']['wave'][0] + 1) & \
        (wave > clv_rm_modelling['common']['wave'][-1] - 1)

    ancillary['norm_convolved_shifted'][wavelength_exclusion] = \
        np.amax(ancillary['norm_convolved_shifted'][wavelength_selection])
    ancillary['stellar_spectra_convolved_shifted'][wavelength_exclusion] = \
        np.amax(ancillary['stellar_spectra_convolved_shifted'][wavelength_selection])

    ancillary['correction'] = ancillary['stellar_spectra_convolved_shifted'] / ancillary['norm_convolved_shifted']

    return ancillary['correction'], ancillary
