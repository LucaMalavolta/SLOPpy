from __future__ import print_function, division
from SLOPpy.subroutines.common import *
from SLOPpy.subroutines.spectral_subroutines import *
from SLOPpy.subroutines.io_subroutines import *
from SLOPpy.subroutines.fit_subroutines import *
from SLOPpy.subroutines.plot_subroutines import *
from SLOPpy.subroutines.shortcuts import *

__all__ = ["check_differential_refraction", "plot_check_differential_refraction", "write_differential_refraction"]

subroutine_name = 'check_differential_refraction'


def check_differential_refraction(config_in):

    night_dict = from_config_get_nights(config_in)

    for night in night_dict:

        try:
            processed = load_from_cpickle('check_differential_refraction_processed', config_in['output'], night)
            check_drc = load_from_cpickle('check_differential_refraction', config_in['output'], night)
            print("{0:45s} Night:{1:15s}   {2:s}".format(subroutine_name, night, 'Retrieved'))
            continue
        except:
            print("{0:45s} Night:{1:15s}   {2:s}".format(subroutine_name, night, 'Computing'))
            print()

        """ Retrieving the list of observations"""
        lists = load_from_cpickle('lists', config_in['output'], night)

        calib_data = load_from_cpickle('calibration_fibA', config_in['output'], night)

        input_data = retrieve_observations(config_in['output'], night, lists['observations'],
                                           use_refraction=False, use_telluric=False)

        input_data_corrected = retrieve_observations(config_in['output'], night, lists['observations'],
                                           use_refraction=True, use_telluric=False)

        input_data_s1d = load_from_cpickle('input_dataset_s1d_fibA', config_in['output'], night)

        observational_pams = load_from_cpickle('observational_pams', config_in['output'], night)

        check_drc = {
            'subroutine': subroutine_name,
            'wave': input_data['coadd']['wave']
        }
        processed = {
            'subroutine': subroutine_name
        }

        for obs in lists['observations']:

            processed[obs] = {
                'n_orders': input_data[obs]['n_orders'],
                'n_pixels': input_data[obs]['n_pixels']
            }

            """ for plotting purpose only"""
            #processed[obs]['e2ds'] = input_data[obs]['e2ds']
            #processed[obs]['e2ds_err'] = input_data[obs]['e2ds_err']

            #processed[obs]['flux'] = input_data[obs]['e2ds']/calib_data['blaze']/input_data[obs]['step']
            #processed[obs]['flux_err'] = np.sqrt(input_data[obs]['e2ds'])/calib_data['blaze']/input_data[obs]['step']
            preserve_flux = input_data[obs].get('absolute_flux', True)

            processed[obs]['flux_s1d'] = \
                rebin_2d_to_1d(input_data[obs]['wave'], input_data[obs]['step'], input_data[obs]['e2ds'],
                               calib_data['blaze'], input_data['coadd']['wave'], input_data['coadd']['step'],
                               preserve_flux=preserve_flux,
                               rv_shift=observational_pams[obs]['rv_shift_ORF2BRF'])

            """
            processed[obs]['flux_s1d_err'] = \
                rebin_2d_to_1d(input_data[obs]['wave'], input_data[obs]['step'],  input_data[obs]['e2ds_err'],
                               calib_data['blaze'], input_data['coadd']['wave'], input_data['coadd']['step'],
                               preserve_flux=preserve_flux,
                               rv_shift=0.00, is_error=True)
            """
            processed[obs]['flux_s1d_err'] = processed[obs]['flux_s1d']


            processed[obs]['s1d_rescaling'], processed[obs]['s1d_rescaled'], processed[obs]['s1d_rescaled_err'] = \
                perform_rescaling(input_data['coadd']['wave'],
                                  processed[obs]['flux_s1d'],
                                  processed[obs]['flux_s1d_err'],
                                  [5450.0, 5550.0])
                                  #observational_pams['wavelength_rescaling'])

            """ for plotting purpose only"""
            #processed[obs]['e2ds_corr'] = input_data_corrected[obs]['e2ds']
            #processed[obs]['e2ds_err_corr'] = input_data_corrected[obs]['e2ds_err']

            #processed[obs]['flux_corr'] = input_data_corrected[obs]['e2ds']/calib_data['blaze']/input_data_corrected[obs]['step']
            #processed[obs]['flux_err_corr'] = np.sqrt(input_data_corrected[obs]['e2ds'])/calib_data['blaze']/input_data_corrected[obs]['step']


            processed[obs]['flux_s1d_corr'] = \
                rebin_2d_to_1d(input_data_corrected[obs]['wave'],
                               input_data_corrected[obs]['step'],
                               input_data_corrected[obs]['e2ds'],
                               calib_data['blaze'],
                               input_data_corrected['coadd']['wave'],
                               input_data_corrected['coadd']['step'],
                               preserve_flux=preserve_flux,
                               rv_shift=observational_pams[obs]['rv_shift_ORF2BRF'])

            """  
            processed[obs]['flux_s1d_corr_err'] = \
                rebin_2d_to_1d(input_data_corrected[obs]['wave'],
                               input_data_corrected[obs]['step'],
                               input_data_corrected[obs]['e2ds_err'],
                               calib_data['blaze'],
                               input_data_corrected['coadd']['wave'],
                               input_data_corrected['coadd']['step'],
                               preserve_flux=preserve_flux,
                               rv_shift=0.00,
                               is_error=True)
            """
            processed[obs]['flux_s1d_corr_err'] = processed[obs]['flux_s1d_corr']

            processed[obs]['s1d_corr_rescaling'], processed[obs]['s1d_corr_rescaled'], processed[obs]['s1d_corr_rescaled_err'] = \
                perform_rescaling(input_data['coadd']['wave'],
                                  processed[obs]['flux_s1d_corr'],
                                  processed[obs]['flux_s1d_corr_err'],
                                  [5450.0, 5550.0])

            processed[obs]['dr_correction'] = processed[obs]['s1d_corr_rescaled']/processed[obs]['s1d_rescaled']

            processed[obs]['s1d_DRS_rescaling'], processed[obs]['s1d_DRS_rescaled'], processed[obs]['s1d_DRS_rescaled_err'] = \
                perform_rescaling(input_data_s1d[obs]['wave'],
                                  input_data_s1d[obs]['flux'],
                                  np.sqrt(np.abs(input_data_s1d[obs]['flux'])),
                                  [5450.0, 5550.0])
                                  #observational_pams['wavelength_rescaling'])

            processed[obs]['DRS_coeff_flux'] = []
            processed[obs]['DRS_corr'] = np.zeros(input_data_s1d[obs]['size'], dtype=np.double)

            for coeff_index in np.arange(0, 10, 1, dtype=np.int16):
                try:
                    keyword = 'HIERARCH TNG DRS FLUX CORR COEFF' + repr(coeff_index)
                    processed[obs]['DRS_coeff_flux'].extend([input_data[obs]['header']['ccf'][keyword]])
                    processed[obs]['DRS_corr'] += \
                        input_data[obs]['header']['ccf'][keyword] * \
                        np.power(input_data_s1d[obs]['wave'], coeff_index)
                except:
                    continue

            processed[obs]['DRS_corr_rescaling'], processed[obs]['DRS_corr_rescaled'], _ = \
                perform_rescaling(input_data_s1d[obs]['wave'],
                                  processed[obs]['DRS_corr'],
                                  processed[obs]['DRS_corr'],
                                  observational_pams['wavelength_rescaling'])

            check_drc[obs] = {
                's1d': {
                    'wave': input_data['coadd']['wave'],
                    'flux': processed[obs]['flux_s1d'],
                    #'flux_err': processed[obs]['flux_s1d_err'],
                    'rescaled': processed[obs]['s1d_rescaled'],
                    #'rescaled_err': processed[obs]['s1d_rescaled_err']
                    },
                's1d_corr': {
                    'correction': processed[obs]['dr_correction'],
                    'correction_rescaled': processed[obs]['dr_correction'],
                    'flux': processed[obs]['flux_s1d_corr'],
                    #'flux_err': processed[obs]['flux_s1d_corr_err'],
                    'rescaled': processed[obs]['s1d_corr_rescaled'],
                    #'rescaled_err': processed[obs]['s1d_corr_rescaled_err']
                    },
                'DRS_s1d':{
                    'wave': input_data_s1d[obs]['wave'],
                    'flux': input_data_s1d[obs]['flux'],
                    #'flux_err': np.sqrt(input_data_s1d[obs]['flux']+0.1),
                    'rescaled': processed[obs]['s1d_DRS_rescaled'],
                    #'rescaled_err': processed[obs]['s1d_DRS_rescaled_err']
                    },
                'DRS_s1d_corr': {
                    'correction': processed[obs]['DRS_corr'],
                    'correction_rescaled': processed[obs]['DRS_corr_rescaled'],
                    'flux': input_data_s1d[obs]['flux']/processed[obs]['DRS_corr'],
                    #'flux_err': np.sqrt(np.abs(input_data_s1d[obs]['flux']))/processed[obs]['DRS_corr'],
                    'rescaled': processed[obs]['s1d_DRS_rescaled']/processed[obs]['DRS_corr_rescaled'],
                    #'rescaled_err': processed[obs]['s1d_DRS_rescaled_err']/processed[obs]['DRS_corr_rescaled'],
                    },
            }

        save_to_cpickle('check_differential_refraction_processed', processed, config_in['output'], night)
        save_to_cpickle('check_differential_refraction', check_drc, config_in['output'], night)

        print('Night ', night, ' completed')
        print()


def plot_check_differential_refraction(config_in, night_input=''):
    night_dict = from_config_get_nights(config_in)

    if night_input == '':
        night_list = night_dict
    else:
        night_list = np.atleast_1d(night_input)

    for night in night_list:

        """ Retrieving the list of observations"""
        lists = load_from_cpickle('lists', config_in['output'], night)

        try:
            """ Retrieving the analysis"""
            check_drc = load_from_cpickle('check_differential_refraction', config_in['output'], night)
        except:
            print("  Failed in retrieving the data")
            return

        observational_pams = load_from_cpickle('observational_pams', config_in['output'], night)

        """ Creation of the color array, based on the BJD of the observations
        """
        bjd = []
        am = []

        for obs in lists['observations']:
            bjd.append(observational_pams[obs]['BJD'] - 2450000.0)
            am.append(observational_pams[obs]['AIRMASS'])

        n_obs = len(lists['observations']) * 1.0 + 1.

        colors = np.asarray(bjd)
        cmap = plt.cm.viridis
        line_colors = cmap(np.linspace(0, 1, len(lists['observations'])))

        fig = plt.figure(figsize=(12, 6))
        gs = GridSpec(2, 2, width_ratios=[50, 1])
        gs.update(hspace=0.1)

        ax1 = plt.subplot(gs[0, 0])
        ax2 = plt.subplot(gs[1, 0], sharex=ax1, sharey=ax1)
        cbax1 = plt.subplot(gs[:, 1])

        for i, obs in enumerate(lists['observations']):

            sel = (check_drc[obs]['s1d']['rescaled'] > -1000000.05)

            ax1.plot(check_drc['wave'][sel],
                     check_drc[obs]['s1d']['rescaled'][sel]+i/5.,
                     c = line_colors[i], lw = 1, alpha = 1)

            ax2.plot(check_drc[obs]['DRS_s1d']['wave'],
                     check_drc[obs]['DRS_s1d']['rescaled']+i/5.,
                     c = line_colors[i], lw = 1, alpha = 1)

            i_max = 1.5 + i/5.

        ax1.set_xlim(check_drc['wave'][0], check_drc['wave'][-1])

        ax1.set_ylim(0.00, i_max)
        ax1.legend(loc=3)
        ax1.set_title('Night: {0:s} \n SLOPpy input s1d'.format(night))
        ax2.set_title('DRS input s1d')

        ax2.set_xlabel('$\lambda$ [$\AA$]')

        sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=colors[0], vmax=colors[-1]))
        sm.set_array([])  # You have to set a dummy-array for this to work...
        cbar = plt.colorbar(sm, cax=cbax1)
        cbar.set_label('BJD - 2450000.0')
        fig.subplots_adjust(wspace=0.05, hspace=0.4)
        plt.show()


        fig = plt.figure(figsize=(12, 6))
        gs = GridSpec(2, 2, width_ratios=[50, 1])
        #gs.update(wspace=0.025, hspace=0.05)
        gs.update(hspace=0.1)

        ax1 = plt.subplot(gs[0, 0])
        ax2 = plt.subplot(gs[1, 0], sharex=ax1, sharey=ax1)
        cbax1 = plt.subplot(gs[:, 1])

        for i, obs in enumerate(lists['observations']):

            sel = (check_drc[obs]['s1d']['flux']> 0.05)

            #ax1.plot(check_drc['wave'][sel],
            #         check_drc[obs]['s1d'][sel]+i/5.,
            #         c=line_colors[i], lw=1, alpha=1.0, zorder=0)

            ax1.plot(check_drc['wave'],
                     check_drc[obs]['s1d_corr']['correction_rescaled']+i/5.,
                     c=line_colors[i], lw=1, alpha=1)

            #ax2.plot(check_drc[obs]['wave_DRS'],
            #         check_drc[obs]['s1d_DRS']+i/5.,
            #         c=line_colors[i], lw=1, alpha=1)

            ax2.plot(check_drc[obs]['DRS_s1d']['wave'],
                     1./check_drc[obs]['DRS_s1d_corr']['correction_rescaled']+i/5.,
                     c=line_colors[i], lw=1, alpha=1)

            i_max = 1.5 + i/5.

        #ax1.plot(processed['coadd']['wave'], processed['coadd']['rescaled'], c='k', lw=1)
        #ax2.plot(processed['coadd']['wave'], processed['coadd']['rescaled'], c='k', lw=1)
        ax1.set_xlim(check_drc['wave'][0], check_drc['wave'][-1])

        ax1.set_ylim(0.00, i_max)
        ax1.legend(loc=3)
        ax1.set_title('Night: {0:s} \n SLOPpy correction function'.format(night))
        ax2.set_title('DRS correction function')

        ax1.set_xlabel('$\lambda$ [$\AA$]')

        sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=colors[0], vmax=colors[-1]))
        sm.set_array([])  # You have to set a dummy-array for this to work...
        cbar = plt.colorbar(sm, cax=cbax1)
        cbar.set_label('BJD - 2450000.0')
        fig.subplots_adjust(wspace=0.05, hspace=0.4)
        plt.show()

        fig = plt.figure(figsize=(12, 6))
        gs = GridSpec(2, 2, width_ratios=[50, 1])
        #gs.update(wspace=0.025, hspace=0.05)
        gs.update(hspace=0.1)

        ax1 = plt.subplot(gs[0, 0])
        ax2 = plt.subplot(gs[1, 0], sharex=ax1, sharey=ax1)
        cbax1 = plt.subplot(gs[:, 1])

        for i, obs in enumerate(lists['observations']):

            sel = (check_drc[obs]['s1d_corr']['rescaled']> 0.05)

            ax1.plot(check_drc['wave'][sel],
                     check_drc[obs]['s1d_corr']['rescaled'][sel]+i/5.,
                     c=line_colors[i], lw=1, alpha=0.5)

            ax2.plot(check_drc[obs]['DRS_s1d']['wave'],
                     check_drc[obs]['DRS_s1d_corr']['rescaled']+i/5.,
                     c=line_colors[i], lw=1, alpha=0.5)

            i_max = 1.5 + i/5.

        #ax1.plot(processed['coadd']['wave'], processed['coadd']['rescaled'], c='k', lw=1)
        #ax2.plot(processed['coadd']['wave'], processed['coadd']['rescaled'], c='k', lw=1)
        ax1.set_xlim(check_drc['wave'][0], check_drc['wave'][-1])

        ax1.set_ylim(0.00, i_max)
        ax1.legend(loc=3)
        ax1.set_title('Night: {0:s} \n SLOPpy corrected spectra'.format(night))
        ax2.set_title('DRS corrected spectra')
        ax2.set_xlabel('$\lambda$ [$\AA$]')

        sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=colors[0], vmax=colors[-1]))
        sm.set_array([])  # You have to set a dummy-array for this to work...
        cbar = plt.colorbar(sm, cax=cbax1)
        cbar.set_label('BJD - 2450000.0')
        fig.subplots_adjust(wspace=0.05, hspace=0.4)
        plt.show()


def write_differential_refraction(config_in):

    night_dict = from_config_get_nights(config_in)

    for night in night_dict:

        print()
        print('write_differential_refraction              Night: ', night)
        print()

        """ Retrieving the list of observations"""
        lists = load_from_cpickle('lists', config_in['output'], night)

        calib_data = load_from_cpickle('calibration_fibA', config_in['output'], night)

        input_data = retrieve_observations(config_in['output'], night, lists['observations'],
                                           use_refraction=True)

        input_data_s1d = load_from_cpickle('input_dataset_s1d_fibA', config_in['output'], night)

        observational_pams = load_from_cpickle('observational_pams', config_in['output'], night)

        # Let's keep it simple to save memory

        dir_SLOPpy_drc = night + '_SLOPpy_drc/'
        dir_DRS_drc = night + '_DRS_drc/'

        os.system('mkdir -p ' + dir_SLOPpy_drc)
        os.system('mkdir -p ' + dir_DRS_drc)

        for obs in lists['observations']:

            processed = dict(n_orders=input_data[obs]['n_orders'], n_pixels=input_data[obs]['n_pixels'])

            preserve_flux = input_data[obs].get('absolute_flux', True)

            processed['flux_s1d_BRF_corr'] = \
                rebin_2d_to_1d(input_data[obs]['wave'], input_data[obs]['step'],input_data[obs]['e2ds'],
                               calib_data['blaze'], input_data['coadd']['wave'], input_data['coadd']['step'],
                               preserve_flux=preserve_flux,
                               rv_shift=observational_pams[obs]['rv_shift_ORF2BRF'])

            processed['flux_s1d_BRF_corr_err'] = \
                rebin_2d_to_1d(input_data[obs]['wave'], input_data[obs]['step'],input_data[obs]['e2ds_err'],
                               calib_data['blaze'], input_data['coadd']['wave'], input_data['coadd']['step'],
                               preserve_flux=preserve_flux,
                               rv_shift=observational_pams[obs]['rv_shift_ORF2BRF'], is_error=True)

            processed['flux_s1d_SRF_corr'] = \
                rebin_2d_to_1d(input_data[obs]['wave'], input_data[obs]['step'],input_data[obs]['e2ds'],
                               calib_data['blaze'], input_data['coadd']['wave'], input_data['coadd']['step'],
                               preserve_flux=preserve_flux,
                               rv_shift=observational_pams[obs]['rv_shift_ORF2SRF'])

            processed['flux_s1d_SRF_corr_err'] = \
                rebin_2d_to_1d(input_data[obs]['wave'], input_data[obs]['step'],input_data[obs]['e2ds_err'],
                               calib_data['blaze'], input_data['coadd']['wave'], input_data['coadd']['step'],
                               preserve_flux=preserve_flux,
                               rv_shift=observational_pams[obs]['rv_shift_ORF2SRF'], is_error=True)


            processed['DRS_coeff_flux'] = []
            processed['DRS_s1d_corr'] = np.zeros(input_data_s1d[obs]['size'], dtype=np.double)
            processed['DRS_e2ds_corr'] = np.zeros(np.shape(input_data[obs]['wave']), dtype=np.double)

            for coeff_index in np.arange(0, 10, 1, dtype=np.int16):
                try:
                    keyword = 'HIERARCH TNG DRS FLUX CORR COEFF' + repr(coeff_index)
                    processed['DRS_coeff_flux'].extend([input_data[obs]['header']['ccf'][keyword]])
                    processed['DRS_s1d_corr'] += \
                        input_data[obs]['header']['ccf'][keyword] * \
                        np.power(input_data_s1d[obs]['wave'], coeff_index)
                    processed['DRS_e2ds_corr'] += \
                        input_data[obs]['header']['ccf'][keyword] * \
                        np.power(input_data[obs]['wave'], coeff_index)
                except:
                    continue

            processed['DRS_s1d_corr'] = input_data_s1d[obs]['flux']/processed['DRS_s1d_corr']
            processed['DRS_e2ds_corr'] = input_data[obs]['e2ds']/processed['DRS_e2ds_corr']

            """Saving the e2ds files"""
            hdu_e2ds_SLOPpy = fits.PrimaryHDU()
            hdu_e2ds_DRS = fits.PrimaryHDU()

            hdu_e2ds_SLOPpy.data = np.asarray(input_data[obs]['e2ds'], dtype=np.float32)
            hdu_e2ds_DRS.data = np.asarray(processed['DRS_e2ds_corr'], dtype=np.float32)

            for key_name, key_val in input_data[obs]['header']['e2ds'].items():

                if key_name == 'SIMPLE' or \
                        key_name=='BITPIX' or \
                        key_name == 'NAXIS' or \
                        key_name == 'NAXIS1' or \
                        key_name == 'NAXIS2':
                    continue

                if len(key_name) > 8:
                    hdu_e2ds_SLOPpy.header['HIERARCH '+ key_name] = key_val
                    hdu_e2ds_DRS.header['HIERARCH '+ key_name] = key_val

                else:
                    hdu_e2ds_SLOPpy.header[key_name] = key_val
                    hdu_e2ds_DRS.header[key_name] = key_val

            hdu_e2ds_SLOPpy.writeto(dir_SLOPpy_drc + obs + '_e2ds_A.fits', overwrite=True)
            hdu_e2ds_DRS.writeto(dir_DRS_drc + obs + '_e2ds_A.fits', overwrite=True)

            """Saving the s1d files"""
            hdu_s1d_SLOPpy_BRF = fits.PrimaryHDU()
            hdu_s1d_SLOPpy_SRF = fits.PrimaryHDU()
            hdu_s1d_DRS = fits.PrimaryHDU()

            hdu_s1d_SLOPpy_BRF.data = np.asarray(processed['flux_s1d_BRF_corr'], dtype=np.float32)
            hdu_s1d_SLOPpy_SRF.data = np.asarray(processed['flux_s1d_SRF_corr'], dtype=np.float32)
            hdu_s1d_DRS.data = np.asarray(processed['DRS_s1d_corr'], dtype=np.float32)

            for key_name, key_val in input_data[obs]['header']['s1d'].items():

                if key_name == 'SIMPLE' or \
                        key_name=='BITPIX' or \
                        key_name == 'NAXIS' or \
                        key_name == 'NAXIS1' or \
                        key_name == 'NAXIS2':
                    continue

                if len(key_name) > 8:
                    hdu_s1d_SLOPpy_BRF.header['HIERARCH '+ key_name] = key_val
                    hdu_s1d_SLOPpy_SRF.header['HIERARCH '+ key_name] = key_val
                    hdu_s1d_DRS.header['HIERARCH '+ key_name] = key_val

                else:
                    hdu_s1d_SLOPpy_BRF.header[key_name] = key_val
                    hdu_s1d_SLOPpy_SRF.header[key_name] = key_val
                    hdu_s1d_DRS.header[key_name] = key_val

            """ Fixing SLOPpy s1d keywords """
            hdu_s1d_SLOPpy_BRF.header['CRVAL1'] = input_data['coadd']['wave'][0]
            hdu_s1d_SLOPpy_BRF.header['CDELT1'] = input_data['coadd']['step'][0]

            hdu_s1d_SLOPpy_SRF.header['CRVAL1'] = input_data['coadd']['wave'][0]
            hdu_s1d_SLOPpy_SRF.header['CDELT1'] = input_data['coadd']['step'][0]

            """ Fixing DRS s1d keywords """
            hdu_s1d_DRS.header['CRVAL1'] = input_data_s1d[obs]['wave'][0]
            hdu_s1d_DRS.header['CDELT1'] = input_data_s1d[obs]['step'][0]

            hdu_s1d_SLOPpy_BRF.writeto(dir_SLOPpy_drc + obs + '_s1d_A.fits', overwrite=True)
            hdu_s1d_SLOPpy_SRF.writeto(dir_SLOPpy_drc + obs + '_s1d_A_SRF.fits', overwrite=True)
            hdu_s1d_DRS.writeto(dir_DRS_drc + obs + '_s1d_A.fits', overwrite=True)

        print()
        print('Night ', night, ' completed')