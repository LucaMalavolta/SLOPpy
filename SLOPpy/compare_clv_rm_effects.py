from __future__ import print_function, division
from SLOPpy.subroutines.common import *
from SLOPpy.subroutines.spectral_subroutines import *
from SLOPpy.subroutines.io_subroutines import *
from SLOPpy.subroutines.fit_subroutines import *
from SLOPpy.subroutines.shortcuts import *

__all__ = ['plot_compare_clv_rm_effects_planetRF',
           'plot_compare_clv_rm_effects_observerRF',
           'plot_compare_clv_rm_effects_stellarRF',
           'plot_compare_clv_rm_effects']


def plot_compare_clv_rm_effects_planetRF(config_in, night_input=''):
    plot_compare_clv_rm_effects(config_in, night_input, reference='planetRF')


def plot_compare_clv_rm_effects_observerRF(config_in, night_input=''):
    plot_compare_clv_rm_effects(config_in, night_input, reference='observerRF')


def plot_compare_clv_rm_effects_stellarRF(config_in, night_input=''):
    plot_compare_clv_rm_effects(config_in, night_input, reference='stellarRF')


def plot_compare_clv_rm_effects(config_in, night_input='', reference='planetRF'):
    
    transmission_average = load_from_cpickle('transmission_average_'+reference, config_in['output'])
    transmission_clv_rm_average = load_from_cpickle('transmission_clv_rm_average_'+reference, config_in['output'])
    night_dict = from_config_get_nights(config_in)

    fig = plt.figure(figsize=(12, 9))

    gs = GridSpec(2, 1)
    ax1 = plt.subplot(gs[0, 0])
    ax2 = plt.subplot(gs[1, 0], sharex=ax1, sharey=ax1)

    ax1.errorbar(transmission_clv_rm_average['wave'],
                 transmission_clv_rm_average['average'],
                 yerr=transmission_clv_rm_average['average_err'],
                 fmt='ko', ms=1, zorder=10, alpha=0.10, label='average, with CLV RM')

    ax1.errorbar(transmission_clv_rm_average['binned_wave'],
                 transmission_clv_rm_average['binned'],
                 yerr=transmission_clv_rm_average['binned_err'],
                 fmt='ko', ms=3, zorder=20, label='binned, with CLV RM')

    #ax1.errorbar(transmission_average['binned_wave'],
    #             transmission_average['binned'],
    #             yerr=transmission_average['binned_err'],
    #             fmt='mo', ms=3, zorder=15, label='binned, no CLV RM')


    ax2.errorbar(transmission_average['wave'],
                 transmission_average['average'],
                 yerr=transmission_average['average_err'],
                 fmt='ko', ms=1, zorder=10, alpha=0.10, label='average, no CLV RM')

    ax2.errorbar(transmission_average['binned_wave'],
                 transmission_average['binned'],
                 yerr=transmission_average['binned_err'],
                 fmt='ko', ms=3, zorder=20, label='binned, no CLV RM')

    #ax2.errorbar(transmission_clv_rm_average['binned_wave'],
    #             transmission_clv_rm_average['binned'],
    #             yerr=transmission_clv_rm_average['binned_err'],
    #             fmt='mo', ms=3, zorder=15, alpha=0.5, label='binned, with CLV RM')

    total_night = len(night_dict)
    side_step = config_in['master-out']['wavelength_step'] * config_in['master-out']['binning_factor'] / 10
    for n_night, night in enumerate(night_dict):
        ax1.errorbar(transmission_clv_rm_average['binned_wave'] + (n_night-total_night/2) * side_step,
                     transmission_clv_rm_average[night]['binned'],
                     yerr=transmission_clv_rm_average[night]['binned_err'],
                     color='C'+repr(n_night), label='{0:s} with CLV RM'.format(night),
                     fmt='o', ms=1, zorder=17, alpha=0.75)

        ax2.errorbar(transmission_average['binned_wave'] + (n_night-total_night/2) * side_step,
                     transmission_average[night]['binned'],
                     yerr=transmission_average[night]['binned_err'],
                     color='C' + repr(n_night), label='{0:s} no CLV RM'.format(night),
                     fmt='o', ms=1, zorder=17, alpha=0.75)


    #ax1.set_ylim(0.95-spec_offset*(1.+n_night), 1.05)
    ax1.set_xlim(config_in['master-out']['wavelength_range'][0], config_in['master-out']['wavelength_range'][1])
    ax1.set_ylim(0.985, 1.01)
    ax2.set_xlabel('$\lambda$ [$\AA$]')
    ax1.legend(loc=3)
    ax1.set_title('Average transmission spectrum with CLV and RM correction, in {0:s}'.format(reference))
    ax2.set_title('Average transmission spectrum without CLV and RM correction, in {0:s}'.format(reference))

    plt.show()

