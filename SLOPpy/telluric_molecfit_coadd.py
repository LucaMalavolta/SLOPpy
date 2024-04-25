from __future__ import print_function, division
from SLOPpy.subroutines.common import *
from SLOPpy.subroutines.spectral_subroutines import *
from SLOPpy.subroutines.io_subroutines import *
from SLOPpy.subroutines.fit_subroutines import *
from SLOPpy.subroutines.plot_subroutines import *
from SLOPpy.subroutines.shortcuts import *
from SLOPpy.telluric_molecfit_preparation import compute_telluric_molecfit_preparation

__all__ = ["compute_telluric_molecfit_coadd",
           "plot_telluric_molecfit_coadd",
           "compute_telluric_molecfit",
           "plot_telluric_molecfit"]

subroutine_name = 'telluric_molecfit_coadd'


def compute_telluric_molecfit_coadd(config_in, no_coadding=False):

    night_dict = from_config_get_nights(config_in)
    instrument_dict = from_config_get_instrument(config_in)
    molecfit_dict = from_config_get_molecfit(config_in)

    compute_telluric_molecfit_preparation(config_in)

    aer_version = molecfit_dict.get('aer_version', '3.8.1.2')

    if no_coadding==True:
        print()
        print("--- WARNING ---: MOLECFIT model will be computed on individual spectra rather than on coadded spectra ")
        print("                 This task is accomplished by the same routine that computes the model on coadded spectra")
        print("                 With the only difference that the coadded spectra is made by a single observations")
        print("                 This choice proves a lot the maintenance of the code. As a consequence you will see a lot")
        print("                 of 'coadded' messages while running the code. ")
        print()
        molecfit_dict['exptime_coadd'] = 0.001
        subroutine_name = 'telluric_molecfit'


    for night in night_dict:

        instrument_name = night_dict[night]['instrument']
        template_dict = instrument_dict[instrument_name]['telluric_template']

        try:
            telluric = load_from_cpickle('telluric', config_in['output'], night)
            print("{0:45s} Night:{1:15s}   {2:s}".format(subroutine_name, night, 'Retrieved'))
            continue
        except:
            print("{0:45s} Night:{1:15s}   {2:s}".format(subroutine_name, night, 'Computing'))
            print()

        print('    instrument :', instrument_name)
        print()

        tellprep = load_from_cpickle('telluric_molecfit_preparation', config_in['output'], night)

        """ Retrieving the list of observations"""
        lists = load_from_cpickle('lists', config_in['output'], night)

        """ Retrieving the observations"""
        calib_data = load_from_cpickle('calibration_fibA', config_in['output'], night)
        input_data = retrieve_observations(config_in['output'], night, lists['observations'], use_telluric=False)
        observational_pams = load_from_cpickle('observational_pams', config_in['output'], night)

        processed = {
            'subroutine': 'telluric_molecfit',
            'n_orders': 0,
            'n_pixels': 0,
        }

        telluric = {
            'subroutine': 'telluric_molecfit',
            'reference_frame': 'observer'
        }

        processed['airmass_ref'] = 0.000
        processed['telluric'] = {}
        processed['rebin'] = {}
        processed['work_dir'] = tellprep['work_dir']

        """
        Molecfit works on pixel grid, so we must ensure that the spectra are rebinned always on the same wavelength
        scale and same wavelength step. We use local arrays for this purpose
        """

        processed['rebin']['wave'] = np.arange(input_data['coadd']['wavelength_range'][0],
                                               input_data['coadd']['wavelength_range'][1],
                                               molecfit_dict['rebinning_step'],
                                               dtype=np.double)

        processed['rebin']['size'] = np.size(processed['rebin']['wave'])
        processed['rebin']['step'] = np.ones(processed['rebin']['size'],
                                             dtype=np.double) * molecfit_dict['rebinning_step']

        processed['rebin'] = {
            'wave': input_data['coadd']['wave'],
            'size': input_data['coadd']['size'],
            'step': input_data['coadd']['step'],
        }

        # TODO: fix the wave:include files
        wave_include = '"'
        for wli_s, wli_e in zip(tellprep['include']['vacuum'][:, 0], tellprep['include']['vacuum'][:, 1]):
            wave_include = wave_include+str(wli_s)+','+str(wli_e)+','
        wave_include = wave_include[:-1]+'"'

        n_coadd = 0
        n_reference = 0
        texp_cumulated = 0.00
        texp_total = 0.000
        coadd_list = []

        # Computing the total integration time
        for n_obs, obs in enumerate(lists['observations']):
            texp_total += input_data[obs]['EXPTIME']

        print('  Writing data and configuration files for molecfit+calctrans')
        print()

        # There must be a more elegant way to do this, but I'm, not aware of it
        for n_obs, obs in enumerate(lists['observations']):

            input_data[obs]['molecfit']['aer_version'] = aer_version

            processed[obs] = {
                'n_orders': input_data[obs]['n_orders'],
                'n_pixels': input_data[obs]['n_pixels']
            }

            """ e2ds spectra are rescaled and then rebinned while keeping them in the Observer Reference Frame"""

            processed[obs]['e2ds_rescaling'], processed[obs]['e2ds_rescaled'], processed[obs]['e2ds_rescaled_err'] = \
                perform_rescaling(input_data[obs]['wave'],
                                  input_data[obs]['e2ds'],
                                  input_data[obs]['e2ds_err'],
                                  observational_pams['wavelength_rescaling'])

            preserve_flux = input_data[obs].get('absolute_flux', True)

            processed[obs]['rebin_ORF'] = \
                rebin_2d_to_1d(input_data[obs]['wave'],
                               input_data[obs]['step'],
                               processed[obs]['e2ds_rescaled'],
                               calib_data['blaze'],
                               processed['rebin']['wave'],
                               processed['rebin']['step'],
                               preserve_flux=preserve_flux,
                               rv_shift=0.00)

            """ This part is relative to the coadded spectrum, must be placed here because
                some variables such as directory names must be defined before the next step
                spectra are coadded to increase the SNR of the spectrum analyzed by molecfit
            """
            if n_coadd == 0:

                reference_name = 'coadded_{0:03d}'.format(n_reference)
                reference_dirname = './' + processed['work_dir'] + '/' + reference_name + '/'
                os.system('mkdir -p ' + reference_dirname)

                rebin_coadd = processed[obs]['rebin_ORF'].copy()

                molecfit_pams = {
                    'MJD': input_data[obs]['MJD'],
                    'UTC': input_data[obs]['UTC'],
                    'ELEVATION': input_data[obs]['ELEVATION'],
                    'HUMIDITY': input_data[obs]['HUMIDITY'],
                    'PRESSURE': input_data[obs]['PRESSURE'],
                    'TEMPERATURE_EN': input_data[obs]['TEMPERATURE_EN'],
                    'TEMPERATURE_M1': input_data[obs]['TEMPERATURE_M1'],
                    'AIRM_START': input_data[obs]['AIRM_START'],
                    'AIRM_END': input_data[obs]['AIRM_END']}

                coadded_files = open(reference_dirname + reference_name + '_files.list', 'w')

                coadd_list.append(reference_name)
                observations_dirlist = []
                observations_exelist = []

            else:
                rebin_coadd += processed[obs]['rebin_ORF']

                molecfit_pams['MJD'] += input_data[obs]['MJD']
                molecfit_pams['UTC'] += input_data[obs]['UTC']
                molecfit_pams['ELEVATION'] += input_data[obs]['ELEVATION']
                molecfit_pams['HUMIDITY'] += input_data[obs]['HUMIDITY']
                molecfit_pams['PRESSURE'] += input_data[obs]['PRESSURE']
                molecfit_pams['TEMPERATURE_EN'] += input_data[obs]['TEMPERATURE_EN']
                molecfit_pams['TEMPERATURE_M1'] += input_data[obs]['TEMPERATURE_M1']
                molecfit_pams['AIRM_END'] = input_data[obs]['AIRM_END']

            n_coadd += 1
            coadded_files.write(obs + '\n')

            texp_cumulated += input_data[obs]['EXPTIME']

            # """ Molecfit analysis is skipped if the telluric correction has been computed already"""
            # if os.path.isfile('./molecfit_'+night +'/output/'+obs+'_ORF_s1d_TAC.dat'):
            # print('    molecfit+calctrans results for ' + obs + ' already available')
            # continue

            """
                This is the directory for MOLECFIT_CALCTRANS and MOLECFIT_CORRECT,
                which is different from the one where the coadded spectrum is saved
            """
            observation_dirname = './' + processed['work_dir'] + '/' + 'obs_{0:03d}'.format(n_obs) + '/'
            os.system('mkdir -p ' + observation_dirname)

            """ the spectrum is saved as a BinTable Fits file in a format suitable for molecfit
                this is the spectrum for MOLECFIT_CALCTRANS and MOLECFIT_CORRECT, so it is saved inside
                the folder with the observation name
            """
            observation_name = obs
            observation_tabname = obs + '_ORF_s1d.fits'

            """ HDU keywords required starting from Molecfit 4.3"""

            hdu_keywords = {
                'MJD-OBS': (input_data[obs]['MJD'], 'MJD-OBS'),
                'ESO OBS EXECTIME': (input_data[obs]['EXPTIME'], 'EXPTIME'), #not sure about this...
                'ESO TEL AIRM START': (input_data[obs]['AIRM_START'], 'AIRMASS START'),
                'ESO TEL AIRM END': (input_data[obs]['AIRM_END'], 'AIRMASS END')
            }

            write_molecfit_input_spectrum(processed['rebin']['wave'],
                                          processed[obs]['rebin_ORF'],
                                          observation_dirname + observation_tabname,
                                          keywords = hdu_keywords)

            observation_calctrans_parname = observation_name + '_calctrans.rc'
            write_calctrans_par(observation_dirname + observation_calctrans_parname)

            """ Writing the SOF files for MOLECFIT_CALCTRANS and MOLECFIT_CORRECT
                For the observed spectrum
            """
            observation_calctrans_sofname = obs + '_calctrans.sof'

            observation_calctrans_soffile = open(observation_dirname + observation_calctrans_sofname, 'w')
            observation_calctrans_soffile.write(observation_tabname+' SCIENCE\n')
            observation_calctrans_soffile.write('../' + reference_name + '/MODEL_MOLECULES.fits MODEL_MOLECULES\n')
            observation_calctrans_soffile.write('../' + reference_name + '/ATM_PARAMETERS.fits ATM_PARAMETERS\n')
            observation_calctrans_soffile.write(
                '../' + reference_name + '/BEST_FIT_PARAMETERS.fits BEST_FIT_PARAMETERS\n')
            observation_calctrans_soffile.close()

            """ Writing the bash script to execute MOLECFIT_CALCTRANS in the directory containing the science fits
            """
            bash_file = './' + processed['work_dir'] + '/calctrans_exec_' + obs + '.source'
            bash_script = open(bash_file, 'w')
            bash_script.write('#!/bin/bash \n')

            bash_script.write('export TMPDIR=$PWD\n')
            bash_script.write('echo  "   " executing calctrans on ' + obs + ' \n')
            bash_script.write('cd ' + observation_dirname + ' \n')

            bash_script.write(molecfit_dict['esorex_exec'] + ' --recipe-config=' + observation_calctrans_parname
                              + ' molecfit_calctrans ' + observation_calctrans_sofname + '> ' + obs + '_calctrans.log\n')
            bash_script.write('cd $TMPDIR \n')
            bash_script.close()

            observations_dirlist.append(observation_dirname)
            observations_exelist.append(bash_file)

            processed[obs]['dir_name'] = observation_dirname
            processed[obs]['tab_name'] = observation_tabname

            if (texp_cumulated >= molecfit_dict['exptime_coadd'] and
                    texp_total-texp_cumulated >= molecfit_dict['exptime_coadd']) \
                    or n_obs == len(lists['observations'])-1:

                coadded_files.close()
                print('   Coadded spectrum: ', n_reference)

                if os.path.exists(reference_dirname + 'TELLURIC_CORR.fits'):
                    print('      molecfit for ' + reference_name + ' previously completed')
                    print()
                else:

                    rebin_coadd /= n_coadd

                    """ the spectra is saved as a FITS file in a format suitable for molecfit """
                    reference_tabname = reference_name + '_ORF_s1d.fits'

                    """ New FITS keywords required by Molecfit 4.3"""
                    hdu_keywords = {
                        'MJD-OBS': (molecfit_pams['MJD'], 'MJD-OBS'),
                        'ESO OBS EXECTIME': (texp_cumulated, 'EXPTIME'),
                        'ESO TEL AIRM START': (molecfit_pams['AIRM_START'], 'AIRMASS START'),
                        'ESO TEL AIRM END': (molecfit_pams['AIRM_END'], 'AIRMASS END')
                    }

                    write_molecfit_input_spectrum(processed['rebin']['wave'],
                                                    rebin_coadd,
                                                    reference_dirname + reference_tabname,
                                                    keywords=hdu_keywords)

                    """ Average of the observational parameters """
                    for key in molecfit_pams:
                        molecfit_pams[key] /= n_coadd

                    molecfit_pams['GEOELEV'] = input_data[obs]['GEOELEV']
                    molecfit_pams['GEOLONG'] = input_data[obs]['GEOLONG']
                    molecfit_pams['GEOLAT'] = input_data[obs]['GEOLAT']

                    reference_molecfit_parname = reference_name + '_molecfit.rc'
                    write_molecfit_par(reference_dirname + reference_molecfit_parname,
                                        wave_include,
                                        input_data[obs]['molecfit'],
                                        molecfit_pams)

                    reference_calctrans_parname = reference_name + '_calctrans.rc'
                    write_calctrans_par(reference_dirname + reference_calctrans_parname)

                    reference_molecfit_sofname = reference_name + '_molecfit.sof'

                    reference_molecfit_soffile = open(reference_dirname + reference_molecfit_sofname, 'w')
                    reference_molecfit_soffile.write(reference_tabname + ' SCIENCE\n')
                    reference_molecfit_soffile.close()

                    """ Writing the SOF files for MOLECFIT_CALCTRANS and MOLECFIT_CORRECT
                        For the observed spectrum
                    """
                    reference_calctrans_sofname = obs + '_calctrans.sof'

                    reference_calctrans_soffile = open(reference_dirname + reference_calctrans_sofname, 'w')
                    reference_calctrans_soffile.write(reference_tabname+' SCIENCE\n')
                    reference_calctrans_soffile.write('MODEL_MOLECULES.fits MODEL_MOLECULES\n')
                    reference_calctrans_soffile.write('ATM_PARAMETERS.fits ATM_PARAMETERS\n')
                    reference_calctrans_soffile.write('BEST_FIT_PARAMETERS.fits BEST_FIT_PARAMETERS\n')
                    reference_calctrans_soffile.close()

                    """ Writing the bash script to execute MOLECFIT_MODEL and MOLECFIT_CALCTRANS in the directory containing the coadded fits
                    """

                    bash_file = './' + processed['work_dir'] + '/molecfit_exec_' + reference_name + '.source'
                    bash_script = open(bash_file, 'w')
                    bash_script.write('#!/bin/bash \n')

                    bash_script.write('export TMPDIR=$PWD\n')
                    bash_script.write('echo  "   " executing molecfit on ' + reference_name + ' \n')
                    bash_script.write('cd ' + reference_dirname + ' \n')

                    bash_script.write(molecfit_dict['esorex_exec'] + ' --recipe-config=' + reference_molecfit_parname
                                      + ' molecfit_model ' + reference_molecfit_sofname + '> ' + obs + '_molecfit.log\n')
                    bash_script.write(molecfit_dict['esorex_exec'] + ' --recipe-config=' + reference_calctrans_parname
                                      + ' molecfit_calctrans ' + reference_calctrans_sofname + '> ' + obs + '_calctrans.log\n')
                    bash_script.write('cd $TMPDIR \n')
                    bash_script.close()

                    os.system('. ' + bash_file)

                for dirname, exename in zip(observations_dirlist, observations_exelist):
                    if os.path.exists(dirname + 'TELLURIC_CORR.fits'):
                        print('      molecfit for ' + dirname + ' previously completed')
                        print()
                    else:
                        os.system('. ' + exename)

                n_coadd = 0
                n_reference += 1
                texp_total -= texp_cumulated
                texp_cumulated = 0.0

        print()

        for n_obs, obs in enumerate(lists['observations']):

            telluric[obs] = {}
            observation_dirname = processed[obs]['dir_name']

            print('  Telluric correction for ', obs, 'retrieved from ', observation_dirname + 'TELLURIC_CORR.fits')

            """ Loading the telluric spectrum from the output directory of molecfit """
            corr_fits = fits.open(observation_dirname + 'TELLURIC_CORR.fits')
            # orig_fits = fits.open(observation_dirname + observation_tabname)
            telluric_molecfit = corr_fits[1].data
            """ rebinning onto the e2ds wave scale"""

            if molecfit_dict.get('fix_telluric', True):
                print('      fix_telluric applied -  temporary workaround for line at 5885.97 A [ORF]')
                line_boundaries = [5885.74, 5886.21]
                sel = (processed['rebin']['wave'] > line_boundaries[0]) \
                    & (processed['rebin']['wave'] < line_boundaries[1])
                tell_cont = np.amax(telluric_molecfit[sel])

                telluric_molecfit[sel] = (telluric_molecfit[sel] - tell_cont) / 2.0 + tell_cont

            telluric[obs]['spectrum'] = \
                rebin_1d_to_2d(processed['rebin']['wave'],
                               processed['rebin']['step'],
                               telluric_molecfit,
                               input_data[obs]['wave'],
                               input_data[obs]['step'],
                               preserve_flux=False)

            try:
                telluric[obs]['spectrum'] = np.nan_to_num(nan=1.0, posinf=1.0, neginf=1.0)
            except:
                temp = ~(np.isfinite(telluric[obs]['spectrum']))
                telluric[obs]['spectrum'][temp] = 1.0
            sel = telluric[obs]['spectrum'] < 0.0001
            telluric[obs]['spectrum'][sel] = 1.0

            telluric[obs]['airmass'] = input_data[obs]['AIRMASS']

            " for compatibilty to some plots, even if it doesn't make any sense"
            telluric[obs]['airmass_ref'] = 0.000
            telluric[obs]['spectrum_noairmass'] = np.power(telluric[obs]['spectrum'],
                                                           telluric[obs]['airmass_ref'] - input_data[obs]['AIRMASS'])
            telluric[obs]['null'] = telluric[obs]['spectrum_noairmass'] < 0.001
            telluric[obs]['spectrum_noairmass'][telluric[obs]['null']] = 1.0
            # we just copy the spectrum file, it's it's a model itself
            telluric[obs]['spline'] = telluric[obs]['spectrum'].copy()

            processed[obs]['e2ds_corrected'] = processed[obs]['e2ds_rescaled'] / telluric[obs]['spectrum']
            processed[obs]['e2ds_corrected_err'] = processed[obs]['e2ds_rescaled_err'] / telluric[obs]['spectrum']

        save_to_cpickle('telluric', telluric, config_in['output'], night)
        save_to_cpickle('telluric_processed', processed, config_in['output'], night)

        print()
        print("Night ", night, " completed")


def plot_telluric_molecfit_coadd(config_in, night_input='', no_coadding=False):
    import matplotlib.pyplot as plt

    night_dict = from_config_get_nights(config_in)
    instrument_dict = from_config_get_instrument(config_in)
    system_dict = from_config_get_system(config_in)

    if no_coadding==True:
        print()
        print("--- WARNING ---: MOLECFIT model will be computed on individual spectra rather than on coadded spectra ")
        print("                 This task is accomplished by the same routine that computes the model on coadded spectra")
        print("                 With the only difference that the coadded spectra is made by a single observations")
        print("                 This choice proves a lot the maintenance of the code. As a consequence you will see a lot")
        print("                 of 'coadded' messages while running the code. ")
        print()
        subroutine_name = 'telluric_molecfit'


    if night_input == '':
        night_list = night_dict
    else:
        night_list = np.atleast_1d(night_input)

    for night in night_list:

        # plt.scatter(rescaling_array, computed_std, c='C0', zorder=1)
        # plt.scatter(sel_factor, sel_stdev, c='C1', zorder=2)
        # plt.plot(rescaling_array, np.polyval(coeff, rescaling_array))
        # plt.plot(rescaling_array, 2*rescaling_array*coeff[0] + coeff[1] )
        # plt.plot()

        print("{0:45s} Night:{1:15s}   {2:s}".format(subroutine_name, night, 'Plotting'))


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

        input_data = retrieve_observations(config_in['output'], night, lists['observations'], use_telluric=False)

        colors, cmap, line_colors = make_color_array(lists, observational_pams)

        fig = plt.figure(figsize=(12, 6))
        gs = GridSpec(2, 2, width_ratios=[50, 1])

        ax1 = plt.subplot(gs[0, 0])
        ax2 = plt.subplot(gs[1, 0], sharex=ax1)
        cbax1 = plt.subplot(gs[:, 1])

        lift_spectrum = 0.25

        for i, obs in enumerate(lists['observations']):
            color_array = cmap(i / len(lists['observations']))

            for order in range(0, processed[obs]['n_orders']):

                if order == 0 and i == 0:
                    ax1.plot(input_data[obs]['wave'][order, :],
                             processed[obs]['e2ds_rescaled'][order, :],
                             c=color_array, lw=1, alpha=0.5, label='uncorrected')
                    ax1.scatter(input_data[obs]['wave'][order, :],
                                processed[obs]['e2ds_corrected'][order, :],
                                s=1, c=np.atleast_2d(color_array), label='corrected')
                else:
                    ax1.plot(input_data[obs]['wave'][order, :],
                             processed[obs]['e2ds_rescaled'][order, :],
                             c=color_array, lw=1, alpha=0.5)
                    ax1.scatter(input_data[obs]['wave'][order, :],
                                processed[obs]['e2ds_corrected'][order, :],
                                s=1, c=np.atleast_2d(color_array))

                # ax1.plot(processed[obs]['wave'][order, :],
                #            e2ds_rescaled[order, :]+lift_spectrum,
                #            c=color_array, lw=1, alpha=0.5)
                # ax1.scatter(processed[obs]['wave'][order, :],
                #            e2ds_rescaled_corrected_spline[order, :]+lift_spectrum,
                #            s=1, c=np.atleast_2d(color_array))

                ax2.plot(input_data[obs]['wave'][order, :],
                         telluric[obs]['spectrum'][order, :],
                         c=color_array)
                ax2.axhline(1.00, c='k')

                # ax2.plot(processed[obs]['wave'][order, :],
                #         telluric[obs]['spline'][order, :]+lift_spectrum,
                #         c=color_array)
                # ax2.axhline(1.00+lift_spectrum, c='k')

        # ax2.plot(input_data['coadd']['wave'],telluric['stellarRF']['spline_eval']+0.1,c='k')
        # ax2.scatter(input_data['coadd']['wave'],telluric['stellarRF']['spectrum']+0.1,c='r', s=2)

        ax1.legend(loc=3)
        ax1.set_title('Night: ' + night)

        ax2.set_xlabel('$\lambda$ [$\AA$]')

        try:
            instrument = night_dict[night]['instrument']
            comparison_file = config_in['instruments'][instrument]['telluric_comparison']
            comparison_data = np.genfromtxt(comparison_file, skip_header=1)
            if comparison_data[0, 0] < 1000.0:
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


def compute_telluric_molecfit(config_in):
    compute_telluric_molecfit_coadd(config_in, no_coadding=True)

def plot_telluric_molecfit(config_in, night_input=''):
    plot_telluric_molecfit_coadd(config_in, night_input=night_input, no_coadding=True)
