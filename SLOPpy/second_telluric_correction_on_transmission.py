from __future__ import print_function, division
from SLOPpy.subroutines.common import *
from SLOPpy.subroutines.spectral_subroutines import *
from SLOPpy.subroutines.fit_subroutines import *
from SLOPpy.subroutines.io_subroutines import *
from SLOPpy.subroutines.shortcuts import *

from scipy.stats import linregress


__all__ = ["compute_second_telluric_correction_on_transmission", "plot_second_telluric_correction_on_transmission"]

subroutine_name = 'second_telluric_correction_on_transmission'


def compute_second_telluric_correction_on_transmission(config_in):
    night_dict = from_config_get_nights(config_in)
    shared_data = load_from_cpickle('shared', config_in['output'])

    for night in night_dict:

        try:
            transmission = load_from_cpickle('transmission_planetRF_second_correction', config_in['output'], night)
            continue
        except:
            print("No transmission spectra with second correction found, computing now ")
            print()

        print()
        print("compute_second_telluric_correction_on_transmission     Night: ", night)

        """ Retrieving the list of observations"""
        lists = load_from_cpickle('lists', config_in['output'], night)

        transmission = load_from_cpickle('transmission_planetRF', config_in['output'], night)

        telluric = load_from_cpickle('telluric', config_in['output'], night)

        calib_data = load_from_cpickle('calibration_fibA', config_in['output'], night)
        input_data = retrieve_observations(config_in['output'], night, lists['observations'], use_telluric=False)

        observational_pams = load_from_cpickle('observational_pams', config_in['output'], night)

        processed = {
            'subroutine': subroutine_name
        }

        for obs in lists['observations']:

            processed[obs] = {}


            processed[obs]['telluric_shifted'] = \
                rebin_2d_to_1d(input_data[obs]['wave'],
                               input_data[obs]['step'],
                               telluric[obs]['spectrum'],
                               telluric[obs]['spectrum'],
                               transmission['wave'],
                               transmission['step'],
                               rv_shift=observational_pams[obs]['rv_shift_ORF2PRF'],
                               preserve_flux=False,
                               skip_blaze_correction=True)


            processed[obs]['selection'] = (np.abs(1.0000-transmission[obs]['rescaled']) < 2*np.std(transmission[obs]['rescaled'])) \
                                          & (np.abs(1.0000-processed[obs]['telluric_shifted']) > 0.02)

            transmission[obs]['slope'], \
            transmission[obs]['intercept'], \
            transmission[obs]['rvalue'], \
            transmission[obs]['pvalue'], \
            transmission[obs]['stderr'], = linregress(processed[obs]['telluric_shifted'][processed[obs]['selection']],
                                                      transmission[obs]['rescaled'][processed[obs]['selection']])

            transmission[obs]['rescaled'] += 1.000 - (transmission[obs]['intercept'] +
                                                               transmission[obs]['slope']*processed[obs]['telluric_shifted'])

        array_average = np.zeros([len(lists['transit_in']), transmission['size']])
        weights_average = np.zeros([len(lists['transit_in']), transmission['size']])
        for i, obs in enumerate(lists['transit_in']):
            array_average[i, :] = transmission[obs]['rescaled'][:]
            weights_average[i, :] = 1./(transmission[obs]['rescaled_err']**2.)


        transmission['average'], transmission['sum_weights'] = np.average(
            array_average, axis=0, weights=weights_average, returned=True)
        transmission['average_err'] = 1./np.sqrt(transmission['sum_weights'])

        array_average = np.zeros([len(lists['transit_out']), transmission['size']])
        weights_average = np.zeros([len(lists['transit_out']), transmission['size']])

        for i, obs in enumerate(lists['transit_out']):
            array_average[i, :] = transmission[obs]['rescaled'][:]
            weights_average[i, :] = 1./(transmission[obs]['rescaled_err']**2.)


        transmission['average_out'], transmission['sum_weights_out'] = np.average(
            array_average, axis=0, weights=weights_average, returned=True)
        transmission['average_out_err'] = 1./np.sqrt(transmission['sum_weights_out'])

        save_to_cpickle('transmission_planetRF_second_correction_processed', processed, config_in['output'], night)
        save_to_cpickle('transmission_planetRF_second_correction', transmission, config_in['output'], night)


def plot_second_telluric_correction_on_transmission(config_in, night_input=''):
    night_dict = from_config_get_nights(config_in)

    if night_input == '':
        night_list = night_dict
    else:
        night_list = np.atleast_1d(night_input)

    for night in night_list:

        """ Retrieving the list of observations"""
        lists = load_from_cpickle('lists', config_in['output'], night)

        previous_transmission = load_from_cpickle('transmission_planetRF', config_in['output'], night)

        transmission = load_from_cpickle('transmission_planetRF_second_correction', config_in['output'], night)
        processed = load_from_cpickle('transmission_planetRF_second_correction_processed', config_in['output'], night)

        f, (ax1, ax2) = plt.subplots(2, sharex=True, sharey=True)

        cmap = plt.cm.Spectral
        line_colors = cmap(np.linspace(0, 1, len(lists['observations'])))

        for i, obs in enumerate(lists['observations']):

            ax1.axhline(1.+i/10., c='k', zorder=0)
            ax2.axhline(1.+i/10., c='k', zorder=0)

            ax1.scatter(processed[obs]['telluric_shifted'], previous_transmission[obs]['rescaled']+i / 10., s=2,
                        c=np.atleast_2d(line_colors[i]), zorder=1)

            ax2.scatter(processed[obs]['telluric_shifted'], transmission[obs]['rescaled'] + i / 10., s=2,
                        #c=np.atleast_2d(line_colors[i]), zorder=2)
                        c = 'r', zorder = 2)

        ax1.set_xlim(0.80, 1.02)
        plt.show()



        """ Creation of the color array, based on the BJD of the observations
        """
        bjd = []
        am = []

        for obs in lists['observations']:
            bjd.append(transmission[obs]['BJD'] - 2450000.0)
            am.append(transmission[obs]['AIRMASS'])

        colors = np.asarray(bjd)
        cmap = plt.cm.Spectral
        line_colors = cmap(np.linspace(0, 1, len(lists['observations'])))

        fig = plt.figure(figsize=(12, 6))
        gs = GridSpec(1, 2, width_ratios=[50, 1])
        ax = plt.subplot(gs[0, 0])
        cbax1 = plt.subplot(gs[:, 1])

        i_shift = 0.00
        for i, obs in enumerate(lists['observations']):
            ax.errorbar(transmission['wave'],
                        transmission[obs]['rescaled']-i_shift,
                        yerr=transmission[obs]['rescaled_err'],
                        marker = 'o', c=line_colors[i], ms=1, alpha=0.5)

            i_shift += 0.05

        ax.set_ylim(0.0-i_shift, 1.2)
        ax.set_xlabel('$\lambda$ [$\AA$]')
        ax.legend(loc=3)
        ax.set_title('Night: ' + night)

        sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=colors[0], vmax=colors[-1]))
        sm.set_array([])  # You have to set a dummy-array for this to work...
        cbar = plt.colorbar(sm, cax=cbax1)
        cbar.set_label('BJD - 2450000.0')
        fig.subplots_adjust(wspace=0.05, hspace=0.4)
        plt.show()

