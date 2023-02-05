from __future__ import print_function, division

from SLOPpy.instruments.common_DRSv3 import *

def HARPNv3_get_instrument_keywords():
    """ These definitions applt to DRS version 3.x """
    keywords = {
        'header_rvc': 'HIERARCH TNG DRS CCF RVC',
        'header_berv': 'HIERARCH TNG DRS BERV',
        'header_bjd': 'HIERARCH TNG DRS BJD',
        'header_mjd': 'MJD-OBS', # MJD in days. This parameter is required for the retrieval of GDAS data

        'header_blaze': 'HIERARCH TNG DRS BLAZE FILE',
        'header_ccd': 'HIERARCH TNG DRS CCD SIGDET',
        'header_conad': 'HIERARCH TNG DRS CCD CONAD',

        'header_dpr_catg': 'HIERARCH TNG DPR CATG',
        'header_dpr_type': 'HIERARCH TNG DPR TYPE',

        'header_deg_ll': 'HIERARCH TNG DRS CAL TH DEG LL',
        'header_coeff_ll': 'HIERARCH TNG DRS CAL TH COEFF LL',

        'airmass_alt_start': 'HIERARCH TNG TEL AIRM START',
        'airmass_alt_end': 'HIERARCH TNG TEL AIRM END',

        ## Telescope altitude is computed using the middle values obtained from airmass
        'humidity':'HIERARCH TNG METEO HUMIDITY', # Relative humidity in % for GEOELEV.
        'pressure':'HIERARCH TNG METEO PRESSURE',
        'temperature_env': 'HIERARCH TNG METEO TEMP10M', #Ambient temperature in C for GEOELEV
        'temperature_m1': 'HIERARCH TNG M1 CH1TEMP', # Temperature of primary mirror M1 in C (for emission spectra only)
    }

    properties = {
        # DRS-specific keywords
        'time_stamp': 'mid_exposure',
        'time_standard': 'UTC',

        # Observatory-specific keywords
        'geoelev': 2387.2, # meters
        'longitude' : -17.889, # Tel geo longitude (+=East) (deg)
        'latitude' : 28.754,  # Tel geo latitute (+=North) (deg)
        # Instrument-specific keyword
        'n_orders_A': 69,
        'n_orders_B': 69,
        'orders_BtoA':
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
             10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
             20, 21, 22, 23, 24, 25, 26, 27, 28, 29,
             30, 31, 32, 33, 34, 35, 36, 37, 38, 39,
             40, 41, 42, 43, 44, 45, 46, 47, 48, 49,
             50, 51, 52, 53, 54, 55, 56, 57, 58, 59,
             60, 61, 62, 63, 64, 65, 66, 67, 68],
            # after many experiments, I found out the easiest and more robust way to define
            # the order correspondence between fiber A anf B is just to write it down
        'red_ccd':
            [        42, 43, 44, 45, 46, 47, 48, 49,
             50, 51, 52, 53, 54, 55, 56, 57, 58, 59,
             60, 61, 62, 63, 64, 65, 66, 67, 68],
        'blue_ccd':
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
             10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
             20, 21, 22, 23, 24, 25, 26, 27, 28, 29,
             30, 31, 32, 33, 34, 35, 36, 37, 38, 39,
             40, 41],
        'full_ccd':
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
             10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
             20, 21, 22, 23, 24, 25, 26, 27, 28, 29,
             30, 31, 32, 33, 34, 35, 36, 37, 38, 39,
             40, 41, 42, 43, 44, 45, 46, 47, 48, 49,
             50, 51, 52, 53, 54, 55, 56, 57, 58, 59,
             60, 61, 62, 63, 64, 65, 66, 67, 68],
        # The following are the input values used by Molecfit, taken from Allart+2017
        # for convenience, all the default values are listed here instead of being scattered into the code
        'molecfit': {
            'default_wstep': 0.01000, # default wavelength step size for the input stellar spectra
            'molecules': ['H2O', 'O2'],
            'ftol': "1e-9",
            'xtol': "1e-9",

            'cont_const': 1.0, # a0,  This value differs from Allart+2017 since we are using normalized spectra
            'cont_n': 3, # n_cont, Degree of coefficients for continuum fit
            'wlc_n': 2, # n_lambda, Polynomial degree of the refined wavelength solution
            'wlc_const': 0.0, # b0, Initial constant term for wavelength correction (shift relative to half  wavelength range)
            'res_gauss': 4.8, # omega_Gaussian,  Initial value for FWHM of Gaussian in pixels
            'kernfac': 15, #kernel_size, Size of Gaussian/Lorentzian/Voigtian kernel in FWHM

            'slitwidth': 1.00, # in arcseconds
            'pixelscale': 0.16,
        }
    }

    return keywords, properties

# Shortcut from DRS-geenral to instrument-specific subroutine
def HARPN_DRSv3_get_calib_data(archive, file_rad, fiber='A', order_selection=None):

    return DRSv3_get_calib_data(archive, file_rad, fiber=fiber, order_selection=order_selection)

# Shortcut from DRS-geenral to instrument-specific subroutine
def HARPN_DRSv3_get_input_data(archive, file_rad, mask, fiber='A', skip_ccf=None, skip_s1d=True, order_selection=None):

    return DRSv3_get_input_data(archive, file_rad, mask, fiber=fiber, skip_ccf=skip_ccf, skip_s1d=skip_s1d, order_selection=order_selection)
