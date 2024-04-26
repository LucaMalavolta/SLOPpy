from __future__ import print_function, division
from SLOPpy.subroutines.common import *
from SLOPpy.subroutines.spectral_subroutines import *
from SLOPpy.subroutines.io_subroutines import *
from SLOPpy.subroutines.fit_subroutines import *
from SLOPpy.subroutines.plot_subroutines import *
from SLOPpy.subroutines.shortcuts import *
from SLOPpy.telluric_molecfit_v1_preparation import compute_telluric_molecfit_v1_preparation

__all__ = ["compute_telluric_molecfit_v1_coadd",
           "plot_telluric_molecfit_v1_coadd",
           "compute_telluric_molecfit_v1",
           "plot_telluric_molecfit_v1"]

subroutine_name = 'telluric_molecfit_v1_coadd'


def compute_telluric_molecfit_v1_coadd(config_in, no_coadding=False):
    """
    Lazy workaround
    :param config_in:
    :param kwargs:
    :return:
    """

    night_dict = from_config_get_nights(config_in)
    instrument_dict = from_config_get_instrument(config_in)
    molecfit_dict = from_config_get_molecfit(config_in)

    compute_telluric_molecfit_v1_preparation(config_in)

    if no_coadding==True:
        print()
        print("--- WARNING ---: MOLECFIT model will be computed on individual spectra rather than on coadded spectra ")
        print("                 This task is accomplished by the same routine that computes the model on coadded spectra")
        print("                 With the only difference that the coadded spectra is made by a single observations")
        print("                 This choice proves a lot the maintenance of the code. As a consequence you will see a lot")
        print("                 of 'coadded' messages while running the code. ")
        print()
        molecfit_dict['exptime_coadd'] = 0.001
        subroutine_name = 'telluric_molecfit_v1'
    else:
        subroutine_name = 'telluric_molecfit_coadd_v1'

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

            """ Molecfit analysis is skipped if the telluric correction has been computed already"""
            # if os.path.isfile('./molecfit_'+night +'/output/'+obs+'_ORF_s1d_TAC.dat'):
            # print('    molecfit+calctrans results for ' + obs + ' already available')
            # continue

            """ the spectra is saved as an ASCII file in a format suitable for molecfit """
            fileout = open('./' + processed['work_dir'] + '/' + obs + '_ORF_s1d.dat', 'w')
            for w, f in zip(processed['rebin']['wave'], processed[obs]['rebin_ORF']):
                fileout.write('{0:12.6f} {1:12.6f} \n'.format(w, f))
            fileout.close()

            """
            processed[obs]['rebin_SRF'] = \
            rebin_2d_to_1d(input_data[obs]['wave'],
                           input_data[obs]['step'],
                           processed[obs]['e2ds_rescaled'],
                           calib_data['blaze'],
                           processed['rebin']['wave'],
                           processed['rebin']['step'],
                           preserve_flux=preserve_flux,
                           rv_shift = observational_pams[obs]['rv_shift_ORF2SRF'])

            fileout = open('./' + processed['work_dir'] + '/' + obs + '_SRF_s1d.dat','w')
            for w, f in zip(processed['rebin']['wave'], processed[obs]['rebin_SRF']):
                fileout.write('{0:12.6f} {1:12.6f} \n'.format(w, f))
            fileout.close()
            """

            """ spectra is coadded to increase the SNR of the spectrum analyzed by molecfit """
            if n_coadd == 0:

                reference_name = 'coadded_{0:03d}'.format(n_reference)

                rebin_coadd = processed[obs]['rebin_ORF'].copy()

                molecfit_pams = {
                    'MJD': input_data[obs]['MJD'],
                    'UTC': input_data[obs]['UTC'],
                    'ELEVATION': input_data[obs]['ELEVATION'],
                    'HUMIDITY': input_data[obs]['HUMIDITY'],
                    'PRESSURE': input_data[obs]['PRESSURE'],
                    'TEMPERATURE_EN': input_data[obs]['TEMPERATURE_EN'],
                    'TEMPERATURE_M1': input_data[obs]['TEMPERATURE_M1']}

                coadded_files = open('./' + processed['work_dir'] + '/' + reference_name + '_files.list', 'w')

                coadd_list.append(reference_name)

            else:
                rebin_coadd += processed[obs]['rebin_ORF']

                molecfit_pams['MJD'] += input_data[obs]['MJD']
                molecfit_pams['UTC'] += input_data[obs]['UTC']
                molecfit_pams['ELEVATION'] += input_data[obs]['ELEVATION']
                molecfit_pams['HUMIDITY'] += input_data[obs]['HUMIDITY']
                molecfit_pams['PRESSURE'] += input_data[obs]['PRESSURE']
                molecfit_pams['TEMPERATURE_EN'] += input_data[obs]['TEMPERATURE_EN']
                molecfit_pams['TEMPERATURE_M1'] += input_data[obs]['TEMPERATURE_M1']

            n_coadd += 1
            coadded_files.write(obs + '\n')

            texp_cumulated += input_data[obs]['EXPTIME']

            # TODO: input from configuration file for molecfit installation path
            bash_script = open('./' + processed['work_dir'] + '/molecfit_exec_' + obs + '.source', 'w')
            bash_script.write('#!/bin/bash \n')

            bash_script.write('export TMPDIR=$PWD\n')
            bash_script.write('echo  "   " executing calctrans on ' + obs + ' \n')
            bash_script.write(molecfit_dict['installation_path'] + 'calctrans ' +
                              obs + '.par > ' + obs + '_calctrans.log\n')
            bash_script.close()

            write_molecfit_v1_par('./' + processed['work_dir'] + '/' + obs + '.par',
                                  obs + '_ORF_s1d.dat',
                                  reference_name,
                                  'include_' + night + '.dat',
                                  input_data[obs]['molecfit'],
                                  input_data[obs])

            if (texp_cumulated >= molecfit_dict['exptime_coadd'] and
                    texp_total-texp_cumulated >= molecfit_dict['exptime_coadd']) \
                    or n_obs == len(lists['observations'])-1:

                coadded_files.close()
                print('   Coadded spectrum: ', n_reference)

                rebin_coadd /= n_coadd

                """ the spectra is saved as an ASCII file in a format suitable for molecfit """
                fileout = open('./' + processed['work_dir'] + '/' + reference_name + '_ORF_s1d.dat', 'w')
                for w, f in zip(processed['rebin']['wave'], rebin_coadd):
                    fileout.write('{0:12.6f} {1:12.6f} \n'.format(w, f))
                fileout.close()

                """ Average of the observational parameters """
                for key in molecfit_pams:
                    molecfit_pams[key] /= n_coadd

                molecfit_pams['GEOELEV'] = input_data[obs]['GEOELEV']
                molecfit_pams['GEOLONG'] = input_data[obs]['GEOLONG']
                molecfit_pams['GEOLAT'] = input_data[obs]['GEOLAT']

                # TODO: input from configuration file for molecfit installation path
                bash_script = open('./' + processed['work_dir'] + '/molecfit_exec_' + reference_name + '.source', 'w')
                bash_script.write('#!/bin/bash \n')
                bash_script.write('export TMPDIR=$PWD\n')

                bash_script.write('echo  "   " executing molecfit+calctrans on ' + reference_name + ' \n')
                bash_script.write(molecfit_dict['installation_path'] + 'molecfit ' +
                                  reference_name + '.par > ' + reference_name + '_molecfit.log\n')
                bash_script.write(molecfit_dict['installation_path'] + 'calctrans ' +
                                  reference_name + '.par > ' + reference_name + '_calctrans.log\n')
                bash_script.close()

                # TODO: cycle with variation in UTC until molecfit exits succesfully
                #
                # while True:
                #
                #   # write parameter file
                #   # execute molecfit
                #   # check if file _tac.fits has been written (= successful run)
                #   if cond:
                #       break
                #

                utc_reference = molecfit_pams['UTC'] * 1.
                utc_incremental = True
                utc_increase = 500.

                while True:

                    if os.path.exists('./' + processed['work_dir'] + '/output/' + reference_name + '_tac.asc'):
                        print('  molecfit for ' + reference_name + ' previously completed')
                        print()
                        break

                    write_molecfit_v1_par('./' + processed['work_dir'] + '/' + reference_name + '.par',
                                          reference_name + '_ORF_s1d.dat',
                                          reference_name,
                                          'include_' + night + '.dat',
                                          input_data[obs]['molecfit'],
                                          molecfit_pams)

                    os.system('cd ' + processed['work_dir'] + '/ && . ./molecfit_exec_' + reference_name + '.source')

                    if os.path.exists('./' + processed['work_dir'] + '/output/' + reference_name + '_tac.asc'):
                        print('  molecfit for ' + reference_name + ' successfully completed')
                        print()
                        break

                    if molecfit_pams['UTC'] > 86400 - utc_increase:
                        utc_incremental = False
                        molecfit_pams['UTC'] = utc_reference

                    if utc_incremental:
                        molecfit_pams['UTC'] += utc_increase
                        print('  molecfit for {0:s} crashed, UTC increased from {1:6.0f} to {2:6.0f} '.format(
                            reference_name, utc_reference, molecfit_pams['UTC']))
                    else:
                        molecfit_pams['UTC'] -= utc_increase
                        print('  molecfit for {0:s} crashed, UTC decreased from {1:6.0f} to {2:6.0f} '.format(
                            reference_name, utc_reference, molecfit_pams['UTC']))

                n_coadd = 0
                n_reference += 1
                texp_total -= texp_cumulated
                texp_cumulated = 0.0

        """
        Execute molecfit runs on all the coadded spectra
        """

        # for reference_name in coadd_list:
        #    os.system('cd  molecfit_' + night + '/ && . ./molecfit_exec_' + reference_name + '.source')
        #
        print()
        print('  molecfit completed')

        for obs in lists['observations']:

            if os.path.exists('./' + processed['work_dir'] + '/output/' + obs + '_ORF_s1d_TAC.dat'):
                print('  skipping calctrans execution for observation ' + obs)
            else:
                print('  calctrans execution for observation ' + obs)
                os.system('cd ' + processed['work_dir'] + '/ && . ./molecfit_exec_' + obs + '.source')

        print()
        print('  calctrans completed')

        for n_obs, obs in enumerate(lists['observations']):

            telluric[obs] = {}

            """ Loading the telluric spectrum from the output directory of molecfit """
            telluric_molecfit = np.genfromtxt(
                './' + processed['work_dir'] + '/output/'+obs+'_ORF_s1d_TAC.dat', usecols=2)
            """ rebinning onto the e2ds wave scale"""

            if molecfit_dict.get('fix_telluric', True):
                print('  fix_telluric applied -  temporary workaround for line at 5885.97 A [ORF]')
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


def plot_telluric_molecfit_v1_coadd(config_in, night_input='', no_coadding=False):
    import matplotlib.pyplot as plt

    night_dict = from_config_get_nights(config_in)
    instrument_dict = from_config_get_instrument(config_in)
    system_dict = from_config_get_system(config_in)

    if night_input == '':
        night_list = night_dict
    else:
        night_list = np.atleast_1d(night_input)

    if no_coadding==True:
        print()
        print("--- WARNING ---: MOLECFIT model will be computed on individual spectra rather than on coadded spectra ")
        print("                 This task is accomplished by the same routine that computes the model on coadded spectra")
        print("                 With the only difference that the coadded spectra is made by a single observations")
        print("                 This choice proves a lot the maintenance of the code. As a consequence you will see a lot")
        print("                 of 'coadded' messages while running the code. ")
        print()
        subroutine_name = 'telluric_molecfit'
    else:
        subroutine_name = 'telluric_molecfit_coadd_v1'

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


def compute_telluric_molecfit_v1(config_in):
    compute_telluric_molecfit_v1_coadd(config_in, no_coadding=True)

def plot_telluric_molecfit_v1(config_in, night_input=''):
    plot_telluric_molecfit_v1_coadd(config_in, night_input=night_input, no_coadding=True)
