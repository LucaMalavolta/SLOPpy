from __future__ import print_function, division
import numpy as np
from astropy.io import fits
from SLOPpy.subroutines.rebin_subroutines import *
from SLOPpy.subroutines.constants import *


def difference_utc2tdb(jd):
    ### jd        date                sec        leap diff
    #  2441317.5 1972-01-01T00:00:00 2272060800 10   42.184
    #  2441499.5 1972-07-01T00:00:00 2287785600 11   43.184
    #  2441683.5 1973-01-01T00:00:00 2303683200 12   44.184
    #  2442048.5 1974-01-01T00:00:00 2335219200 13   45.184
    #  2442413.5 1975-01-01T00:00:00 2366755200 14   46.184
    #  2442778.5 1976-01-01T00:00:00 2398291200 15   47.184
    #  2443144.5 1977-01-01T00:00:00 2429913600 16   48.184
    #  2443509.5 1978-01-01T00:00:00 2461449600 17   49.184
    #  2443874.5 1979-01-01T00:00:00 2492985600 18   50.184
    #  2444239.5 1980-01-01T00:00:00 2524521600 19   51.184
    #  2444786.5 1981-07-01T00:00:00 2571782400 20   52.184
    #  2445151.5 1982-07-01T00:00:00 2603318400 21   53.184
    #  2445516.5 1983-07-01T00:00:00 2634854400 22   54.184
    #  2446247.5 1985-07-01T00:00:00 2698012800 23   55.184
    #  2447161.5 1988-01-01T00:00:00 2776982400 24   56.184
    #  2447892.5 1990-01-01T00:00:00 2840140800 25   57.184
    #  2448257.5 1991-01-01T00:00:00 2871676800 26   58.184
    #  2448804.5 1992-07-01T00:00:00 2918937600 27   59.184
    #  2449169.5 1993-07-01T00:00:00 2950473600 28   60.184
    #  2449534.5 1994-07-01T00:00:00 2982009600 29   61.184
    #  2450083.5 1996-01-01T00:00:00 3029443200 30   62.184
    #  2450630.5 1997-07-01T00:00:00 3076704000 31   63.184
    #  2451179.5 1999-01-01T00:00:00 3124137600 32   64.184
    #  2453736.5 2006-01-01T00:00:00 3345062400 33   65.184
    #  2454832.5 2009-01-01T00:00:00 3439756800 34   66.184
    #  2456109.5 2012-07-01T00:00:00 3550089600 35   67.184
    #  2457204.5 2015-07-01T00:00:00 3644697600 36   68.184
    #  2457754.5 2017-01-01T00:00:00 3692217600 37   69.184

    jd_table = np.asarray([2441317.5, 2441499.5, 2441683.5, 2442048.5, 2442413.5, 2442778.5, 2443144.5, 2443509.5, 2443874.5, 2444239.5,
                           2444786.5, 2445151.5, 2445516.5, 2446247.5, 2447161.5, 2447892.5, 2448257.5, 2448804.5, 2449169.5, 2449534.5,
                           2450083.5, 2450630.5, 2451179.5, 2453736.5, 2454832.5, 2456109.5, 2457204.5, 2457754.5])

    df_table = np.asarray([42.184, 43.184, 44.184, 45.184, 46.184, 47.184, 48.184, 49.184, 50.184, 51.184,
                52.184, 53.184, 54.184, 55.184, 56.184, 57.184, 58.184, 59.184, 60.184, 61.184,
                62.184, 63.184, 64.184, 65.184, 66.184, 67.184, 68.184, 69.184])/86400.

    return df_table[(jd_table-jd<0)][-1]

def get_instrument_keywords(instrument):

    if instrument == 'HARPS-N':
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



    elif instrument == 'HARPS':
        keywords = {
            'header_rvc': 'HIERARCH ESO DRS CCF RVC',
            'header_berv': 'HIERARCH ESO DRS BERV',
            'header_bjd': 'HIERARCH ESO DRS BJD',
            'header_mjd': 'MJD-OBS', # MJD in days. This parameter is required for the retrieval of GDAS data

            'header_blaze': 'HIERARCH ESO DRS BLAZE FILE',
            'header_ccd': 'HIERARCH ESO DRS CCD SIGDET',
            'header_conad': 'HIERARCH ESO DRS CCD CONAD',

            'header_dpr_catg': 'HIERARCH ESO DPR CATG',
            'header_dpr_type': 'HIERARCH ESO DPR TYPE',

            'header_deg_ll': 'HIERARCH ESO DRS CAL TH DEG LL',
            'header_coeff_ll': 'HIERARCH ESO DRS CAL TH COEFF LL',

            'airmass_alt_start': 'HIERARCH ESO TEL AIRM START',
            'airmass_alt_end': 'HIERARCH ESO TEL AIRM END',

            ## Telescope altitude is computed using the middle values obtained from airmass
            'humidity':'HIERARCH ESO TEL AMBI RHUM', # Relative humidity in % for GEOELEV.
            #'pressure_start' : 'HIERARCH ESO TEL AMBI PRES START',
            #'pressure_end': 'HIERARCH ESO TEL AMBI PRES END',
            'pressure':'HIERARCH ESO TEL AMBI PRES END',
            'temperature_env': 'HIERARCH ESO TEL AMBI TEMP', #Ambient temperature in C for GEOELEV
            'temperature_m1': 'HIERARCH ESO TEL TH M1 TEMP', # Temperature of primary mirror M1 in C (for emission spectra only)

        }

        properties = {
            # DRS-specific keywords
            'time_stamp': 'mid_exposure',
            'time_standard': 'UTC',

            # Observatory-specific keywords
            'geoelev': 2400.0, # meters
            'longitude' : -70.7345, # Tel geo longitude (+=East) (deg)
            'latitude' : -29.2584,  # Tel geo latitute (+=North) (deg)
            # Instrument-specific keyword
            'n_orders_A': 72,
            'n_orders_B': 71,
            'orders_BtoA':
                [ 0,  1,  2,  3,  4,  5,  6,  7,  8,  9,
                 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
                 20, 21, 22, 23, 24, 25, 26, 27, 28, 29,
                 30, 31, 32, 33, 34, 35, 36, 37, 38, 39,
                 40, 41, 42, 43, 44, -1, 45, 46, 47, 48,
                 49, 50, 51, 52, 53, 54, 55, 56, 57, 58,
                 59, 60, 61, 62, 63, 64, 65, 66, 67, 68,
                 69, 70],
            # after many experiments, I found out the easiest and more robust way to define
            # the order correspondence between fiber A anf B is just to write it down
            'red_ccd':
                [                    45, 46, 47, 48, 49,
                 50, 51, 52, 53, 54, 55, 56, 57, 58, 59,
                 60, 61, 62, 63, 64, 65, 66, 67, 68, 69,
                 70, 71],
            'blue_ccd':
                [0,  1,  2,  3,  4,  5,  6,  7,  8,  9,
                 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
                 20, 21, 22, 23, 24, 25, 26, 27, 28, 29,
                 30, 31, 32, 33, 34, 35, 36, 37, 38, 39,
                 40, 41, 42, 43, 44],
            'full_ccd':
                [0,  1,  2,  3,  4,  5,  6,  7,  8,  9,
                 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
                 20, 21, 22, 23, 24, 25, 26, 27, 28, 29,
                 30, 31, 32, 33, 34, 35, 36, 37, 38, 39,
                 40, 41, 42, 43, 44, 45, 46, 47, 48, 49,
                 50, 51, 52, 53, 54, 55, 56, 57, 58, 59,
                 60, 61, 62, 63, 64, 65, 66, 67, 68, 69,
                 70, 71],
            # The following are the input values used by Molecfit, taken from Allart+2017
            # for convenience, all the default values are listed here instead of being scattered into the code
            'molecfit': {
                'default_wstep': 0.01000, # default wavelength step size for the input stellar spectra
                'molecules': ['H2O', 'O2'],
                'ftol': 1e-9,
                'xtol': 1e-9,

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

    else:
        raise ValueError("Instrument not supported")

    return keywords, properties


def map_orders_AB(properties, order_selection):

    n_orders_A = properties['n_orders_A']
    map_orders_A_full = np.arange(0, n_orders_A, dtype=np.int16)
    map_orders_BtoA = np.asarray(properties['orders_BtoA'])

    map_orders_A = []
    map_orders_B = []

    #if order_selection == 'red_ccd':
    #    working_A_full = map_orders_A_full[properties['red_ccd']]
    #    working_BtoA = map_orders_BtoA[properties['red_ccd']]
    if order_selection in ['red_ccd', 'blue_ccd', 'full_ccd']:
        working_A_full = map_orders_A_full[properties[order_selection]]
        working_BtoA = map_orders_BtoA[properties[order_selection]]

    elif order_selection:
        working_A_full = map_orders_A_full[order_selection]
        working_BtoA = map_orders_BtoA[order_selection]
    else:
        working_A_full = map_orders_A_full
        working_BtoA = map_orders_BtoA

    for order_A, order_B in zip(working_A_full, working_BtoA):
        if order_B < 0: continue
        map_orders_A.extend([order_A])
        map_orders_B.extend([order_B])

    return map_orders_A, map_orders_B


def give_back_selected_orders(properties, fiber, order_selection):

    map_orders_A, map_orders_B = map_orders_AB(properties, order_selection)
    if fiber != 'A':
        return map_orders_B
    else:
        return map_orders_A

def get_calib_data(instrument, archive, file_rad, fiber='A', order_selection=None):

    calib_dict = {}

    keywords, properties = get_instrument_keywords(instrument)
    selected_orders = give_back_selected_orders(properties, fiber, order_selection)

    map_orders_A, map_orders_B = map_orders_AB(properties, order_selection)

    if fiber=='A':
        calib_dict['fibAB_orders_match'] = map_orders_A - np.min(selected_orders)
        """ The first selected order could not have a match in fiber B, so we need to renumber from the first order of
            the input selection, not from the first order that had a match """
    else:
        calib_dict['fibAB_orders_match'] = map_orders_B - np.min(map_orders_B)
        """ Only the orders with a match in fiber A are read in the first place, so we can safely rescale with respect
            to the number of the first order in the matched list """

    e2ds_fits = fits.open(archive+'/'+file_rad+'_e2ds_'+fiber+'.fits')

    if e2ds_fits[0].header[keywords['header_dpr_catg']] != 'SCIENCE':
        return

    try:
        blaze_file = e2ds_fits[0].header[keywords['header_blaze']]
        blaze_fits = fits.open(archive + '/' + blaze_file)
    except:
        blaze_file = e2ds_fits[0].header[keywords['header_blaze']].replace(':', '-')
        blaze_fits = fits.open(archive + '/' + blaze_file)

    # getting blaze file
    calib_dict['blaze'] = blaze_fits[0].data[selected_orders, :]

    # getting lamp file
    try:
        lamp_fits = fits.open(archive + '/' + blaze_file[:29] + '_lamp_' + fiber + '.fits')
        calib_dict['lamp'] = lamp_fits[0].data[selected_orders, :]
        lamp_fits.close()
    except:
        print("lamp files not available, sky correction will not be performed")

    calib_dict['n_pixels'] = blaze_fits[0].header['NAXIS1']
    calib_dict['n_orders'] = len(selected_orders)

    blaze_fits.close()

    return calib_dict


def get_input_data(instrument, archive, file_rad, mask, fiber='A', skip_ccf=None, skip_s1d=True, order_selection=None):

    input_dict = {'mask': mask, 'header':{}}
    input_s1d = {'header':{}}

    keywords, properties = get_instrument_keywords(instrument)
    selected_orders = give_back_selected_orders(properties, fiber, order_selection)

    if mask is None:
        skip_ccf = True

    e2ds_fits = fits.open(archive+'/'+file_rad+'_e2ds_'+fiber+'.fits')
    input_dict['header']['e2ds'] = e2ds_fits[0].header

    input_dict['n_pixels'] = e2ds_fits[0].header['NAXIS1']
    input_dict['n_orders'] = len(selected_orders)

    input_dict['DPR_CATG'] = e2ds_fits[0].header[keywords['header_dpr_catg']]

    if input_dict['DPR_CATG'] != 'SCIENCE':
        return
    input_dict['DPR_TYPE'] = e2ds_fits[0].header[keywords['header_dpr_type']]

    if not skip_s1d:
        s1d_fits = fits.open(archive + '/' + file_rad + '_s1d_'+fiber+'.fits')
        input_dict['header']['s1d'] = s1d_fits[0].header
        input_s1d['header']['s1d'] = s1d_fits[0].header

        temp_wave, temp_step = get_s1d_wave(s1d_fits)
        sel_wave = (temp_wave >= 3879.99990) & (temp_wave <= 6900.0001)

        input_s1d['flux'] = s1d_fits[0].data[sel_wave]
        input_s1d['wave'] = temp_wave[sel_wave]
        input_s1d['step'] = temp_step[sel_wave]
        input_s1d['size'] = np.size(input_s1d['wave'])

        s1d_fits.close()

    if not skip_ccf:
        ccf_fits = fits.open(archive+'/'+file_rad+'_ccf_'+mask+'_'+fiber+'.fits')
        input_dict['RVC'] = ccf_fits[0].header[keywords['header_rvc']]
        input_dict['header']['ccf'] = ccf_fits[0].header

    try:
        input_dict['BERV'] = e2ds_fits[0].header[keywords['header_berv']]
        input_dict['EXPTIME'] = e2ds_fits[0].header['EXPTIME']

        if properties['time_stamp'] == 'start_exposure':
            input_dict['BJD'] = e2ds_fits[0].header[keywords['header_bjd']] + input_dict['EXPTIME']/86400.
            input_dict['MJD'] = e2ds_fits[0].header[keywords['header_mjd']] + input_dict['EXPTIME']/86400.
        elif properties['time_stamp'] == 'mid_exposure':
            input_dict['BJD'] = e2ds_fits[0].header[keywords['header_bjd']]
            input_dict['MJD'] = e2ds_fits[0].header[keywords['header_mjd']]
        elif properties['time_stamp'] == 'end_exposure':
            input_dict['BJD'] = e2ds_fits[0].header[keywords['header_bjd']] - input_dict['EXPTIME']/86400.
            input_dict['MJD'] = e2ds_fits[0].header[keywords['header_mjd']] - input_dict['EXPTIME']/86400.
        else:
            print('*** please specify the relationship between epoch and exposure time - assuming mid-exposure epochs')
            input_dict['BJD'] = e2ds_fits[0].header[keywords['header_bjd']]
            input_dict['MJD'] = e2ds_fits[0].header[keywords['header_mjd']]

        if properties['time_standard'] == 'UTC':
            input_dict['BJD']+=  difference_utc2tdb(input_dict['MJD']+2400000.5)


        input_dict['LST'] = e2ds_fits[0].header['LST']

        try:
	        input_dict['AIRMASS'] = e2ds_fits[0].header['AIRMASS']
        except:
            input_dict['AIRMASS'] = (e2ds_fits[0].header[keywords['airmass_alt_start']]
										+ e2ds_fits[0].header[keywords['airmass_alt_end']])/2.

        input_dict['UTC'] = (input_dict['MJD'] - int(input_dict['MJD'])) * 86400.
        input_dict['HUMIDITY'] = e2ds_fits[0].header[keywords['humidity']]
        input_dict['PRESSURE'] = e2ds_fits[0].header[keywords['pressure']]
        input_dict['TEMPERATURE_EN'] = e2ds_fits[0].header[keywords['temperature_env']]
        input_dict['TEMPERATURE_M1'] = e2ds_fits[0].header[keywords['temperature_m1']]
        input_dict['ELEVATION'] = np.arcsin(1./input_dict['AIRMASS']) * (180./np.pi)

        input_dict['GEOELEV'] = properties['geoelev']
        input_dict['GEOLONG'] = properties['longitude']
        input_dict['GEOLAT'] = properties['latitude']

        input_dict['molecfit'] = properties['molecfit']

        try:
            try:
                input_dict['RA'] = e2ds_fits[0].header['RA-DEG']
                input_dict['DEC'] = e2ds_fits[0].header['DEC-DEG']
            except:
                input_dict['RA'] = e2ds_fits[0].header['RA-RAD'] * 180.00 / np.pi
                input_dict['DEC'] = e2ds_fits[0].header['DEC-RAD'] * 180.00 / np.pi  # weird choice of using DEC in hours
        except:
                input_dict['RA'] = e2ds_fits[0].header['RA']
                input_dict['DEC'] = e2ds_fits[0].header['DEC']

    except:
        print('Keyword error in prepare_dataset - check the FITS header of your files')
        quit()
        pass

    input_dict['BLAZE_file'] = e2ds_fits[0].header[keywords['header_blaze']]
    input_dict['CCD_SIGDET'] = e2ds_fits[0].header[keywords['header_ccd']]
    input_dict['CCD_GAIN'] = e2ds_fits[0].header[keywords['header_conad']]

    # getting data
    input_dict['e2ds'] = e2ds_fits[0].data[selected_orders, :]
    temp_wave, temp_step = get_e2ds_wave(e2ds_fits, keywords['header_deg_ll'], keywords['header_coeff_ll'])

    input_dict['wave'] = temp_wave[selected_orders, :]
    input_dict['step'] = temp_step[selected_orders, :]
    input_dict['orders'] = len(selected_orders)

    input_dict['wave_size'] = e2ds_fits[0].header['NAXIS1']

    e2ds_fits.close()
    if not skip_ccf:
        ccf_fits.close()

    return input_dict,input_s1d


def get_s1d_wave(s1d_fits):
    return np.arange(0, s1d_fits[0].header['NAXIS1'], 1.)*s1d_fits[0].header['CDELT1'] + s1d_fits[0].header['CRVAL1'], \
           np.ones(s1d_fits[0].header['NAXIS1'])*s1d_fits[0].header['CDELT1']


def get_e2ds_wave(e2ds_fits, header_deg_ll, header_coeff_ll, order=None):

    e2ds_o = e2ds_fits[0].header['NAXIS2']
    e2ds_w = e2ds_fits[0].header['NAXIS1']

    e2ds_wave = np.zeros([e2ds_o, e2ds_w], dtype=np.double)
    e2ds_step = np.zeros([e2ds_o, e2ds_w], dtype=np.double)

    d = e2ds_fits[0].header[header_deg_ll]
    x = np.arange(0, e2ds_w, 1.)

    for n in range(0, e2ds_o):
        for i in range(d, -1, -1):
            a_sel = i + n*(1+d)
            a_coeff = e2ds_fits[0].header[header_coeff_ll+repr(a_sel)]
            if i == d:
                y_w = a_coeff
                y_s = i*a_coeff
            else:
                y_w = y_w*x + a_coeff
                if i > 0: y_s = y_s*x + i*a_coeff
        e2ds_wave[n, :] = y_w
        e2ds_step[n, :] = y_s

    if order is None:
        return e2ds_wave, e2ds_step
    else:
        return e2ds_wave[order, :], e2ds_step[order, :]


#def shift_wavelength(wave, step, rv_shift):
#    wave_shift = rv_shift/(speed_of_light/1000.000) + 1.00000
#    return wave*wave_shift, step*wave_shift
#
#
#def shift_wavelength_to_rest(wave, step, rv_shift):
#    inverse_wave_shift = (-rv_shift)/(speed_of_light/1000.000) + 1.00000
#    return wave/inverse_wave_shift, step/inverse_wave_shift
