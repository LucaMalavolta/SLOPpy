from __future__ import print_function, division
from SLOPpy.subroutines.common import *
from SLOPpy.subroutines.spectral_subroutines import *
from SLOPpy.subroutines.io_subroutines import *
from SLOPpy.subroutines.fit_subroutines import *
from SLOPpy.subroutines.shortcuts import *
from SLOPpy.subroutines.rebin_subroutines import *
from astropy.convolution import convolve, Box1DKernel

__all__ = ['compute_transmission_lightcurve_average_planetRF',
           'plot_transmission_lightcurve_average_planetRF',
           'compute_transmission_lightcurve_average_stellarRF',
           'plot_transmission_lightcurve_average_stellarRF',
           'compute_transmission_lightcurve_average_observerRF',
           'plot_transmission_lightcurve_average_observerRF',
           'compute_transmission_lightcurve_average',
           'plot_transmission_lightcurve_average'
           ]

def compute_transmission_lightcurve_average_planetRF(config_in, lines_label):
    compute_transmission_lightcurve_average(config_in, lines_label, reference='planetRF')


def plot_transmission_lightcurve_average_planetRF(config_in, night_input):
    plot_transmission_lightcurve_average(config_in, night_input, reference='planetRF')


def compute_transmission_lightcurve_average_stellarRF(config_in, lines_label):
    compute_transmission_lightcurve_average(config_in, lines_label, reference='stellarRF')


def plot_transmission_lightcurve_average_stellarRF(config_in, night_input):
    plot_transmission_lightcurve_average(config_in, night_input, reference='stellarRF')


def compute_transmission_lightcurve_average_observerRF(config_in, lines_label):
    compute_transmission_lightcurve_average(config_in, lines_label, reference='observerRF')


def plot_transmission_lightcurve_average_observerRF(config_in, night_input):
    plot_transmission_lightcurve_average(config_in, night_input, reference='observerRF')


subroutine_name = 'transmission_lightcurve_average'
pick_files = 'transmission_lightcurve'

def compute_transmission_lightcurve_average(config_in, lines_label, reference='planetRF'):


    night_dict = from_config_get_nights(config_in)
    #instrument_dict = from_config_get_instrument(config_in)
    #system_dict = from_config_get_system(config_in)
    planet_dict = from_config_get_planet(config_in)
    spectral_lines = from_config_get_spectral_lines(config_in)
    lines_dict = spectral_lines[lines_label] # from_config_get_transmission_lightcurve(config_in)

    output_list = ['user',
                    'mcmc_night_MED',
                    'mcmc_night_MAP',
                    'mcmc_global_MED',
                    'mcmc_global_MAP'
                    'user_uncorrected']

    append_list = ['', '_uncorrected', '_clv_model']

    shared_data = load_from_cpickle('shared', config_in['output'])

    """ Using the line-specific range to define the transmission spectrum region """

    shared_selection = (shared_data['coadd']['wave'] >= lines_dict['range'][0]) \
        & (shared_data['coadd']['wave'] < lines_dict['range'][1])
    binned_selection = (shared_data['binned']['wave'] >= lines_dict['range'][0]) \
        & (shared_data['binned']['wave'] < lines_dict['range'][1])

    lightcurve_average_template = {
        'subroutine': subroutine_name,
        'range': lines_dict['range'],
        'wave': shared_data['coadd']['wave'][shared_selection],
        'step': shared_data['coadd']['step'][shared_selection],
        'size': int(np.sum(shared_selection)),
        'binned_wave': shared_data['binned']['wave'][binned_selection],
        'binned_step': shared_data['binned']['step'][binned_selection],
        'binned_size': int(np.sum(binned_selection))
    }

    for output_selection in output_list:
        skip_iteration = False

        try:
            lightcurve_average = load_from_cpickle(subroutine_name+'_'+reference+ '_' + output_selection, config_in['output'], lines=lines_label)
            continue
        except (FileNotFoundError, IOError):
            print("  No average transmission lightcurve found for case:{0:s}, computing now ".format(output_selection))
            print()


        lightcurve_average = lightcurve_average_template.copy()

        # doublet sodium in the lab reference frame

        """
        C stabds for central
        """
        C_bands = {}
        for passband_key, passband_val in lines_dict['passbands'].items():
            C_bands[passband_key] = {}
            for line_key, line_val in lines_dict['lines'].items():
                C_bands[passband_key][line_key] = (np.abs(lightcurve_average['wave'] - line_val)*2. <passband_val)

        """
        S stand for side
        """
        S_bands = {}
        for band_key, band_val in lines_dict['continuum'].items():
            S_bands[band_key] = (lightcurve_average['wave'] >= band_val[0]) & (lightcurve_average['wave'] <= band_val[1])


        if 'full_transit_duration' in planet_dict:
            full_transit_duration = planet_dict['total_transit_duration'][0]
        else:
            full_transit_duration = planet_dict['transit_duration'][0]

        if 'total_transit_duration' in planet_dict:
            total_transit_duration = planet_dict['total_transit_duration'][0]
        else:
            total_transit_duration = planet_dict['transit_duration'][0]

        transit_in_bins = np.linspace(
            -total_transit_duration/2./planet_dict['period'][0],
            total_transit_duration/2./planet_dict['period'][0],
            6
        )
        transit_full_bins = np.linspace(
            -full_transit_duration/2./planet_dict['period'][0],
            full_transit_duration/2./planet_dict['period'][0],
            6
        )

        transit_in_step = np.average(transit_in_bins[1:]-transit_in_bins[:-1])
        transit_full_step = np.average(transit_full_bins[1:]-transit_full_bins[:-1])


        lightcurve_average['transit_in_flag'] = []
        lightcurve_average['transit_full_flag'] = []
        lightcurve_average['transit_out_flag'] = []
        lightcurve_average['transit_in'] = {}
        lightcurve_average['transit_full'] = {}
        lightcurve_average['transit_out'] = {}
        lightcurve_average['observations'] = {'phase': []}
        lightcurve_average['bands_list'] = []
        lightcurve_average['C_bands'] = C_bands
        lightcurve_average['S_bands'] = S_bands
        lightcurve_average['bins'] = {
                'transit_in_bins': transit_in_bins,
                'transit_in_step': transit_in_step,
                'transit_full_bins': transit_full_bins,
                'transit_full_step': transit_full_step
            }

        for band_key in C_bands:
            for name_append in append_list:
                lightcurve_average['observations']['delta_' + band_key + name_append] = []
            lightcurve_average['bands_list'].extend([band_key])

        for night in night_dict:

            try:
                lightcurve = load_from_cpickle(pick_files+'_'+reference+ '_' + output_selection, config_in['output'], night, lines_label)
            except:
                skip_iteration = True
                continue

            print("compute_transmission_lightcurve            Night: ", night)


            lightcurve_average['observations']['phase'].extend(lightcurve['arrays']['observations']['phase'].tolist())

            lightcurve_average['transit_in_flag'].extend(
                lightcurve['arrays']['observations']['transit_in_flag'].tolist())
            lightcurve_average['transit_full_flag'].extend(
                lightcurve['arrays']['observations']['transit_full_flag'].tolist())
            lightcurve_average['transit_out_flag'].extend(
                lightcurve['arrays']['observations']['transit_out_flag'].tolist())

            for band_key in lightcurve_average['bands_list']:
                for name_append in append_list:
                    lightcurve_average['observations']['delta_' + band_key + name_append].extend(
                        lightcurve['arrays']['observations']['delta_' + band_key + name_append].tolist())


        if skip_iteration: continue

        sorting_index = np.argsort(lightcurve_average['observations']['phase'])
        lightcurve_average['observations']['phase'] = np.asarray(lightcurve_average['observations']['phase'])[sorting_index]
        lightcurve_average['transit_in_flag'] = np.asarray(lightcurve_average['transit_in_flag'])[sorting_index]
        lightcurve_average['transit_full_flag'] = np.asarray(lightcurve_average['transit_full_flag'])[sorting_index]
        lightcurve_average['transit_out_flag'] = np.asarray(lightcurve_average['transit_out_flag'])[sorting_index]

        lightcurve_average['transit_in']['phase'] = \
            lightcurve_average['observations']['phase'][lightcurve_average['transit_in_flag']]
        lightcurve_average['transit_full']['phase'] = \
            lightcurve_average['observations']['phase'][lightcurve_average['transit_full_flag']]
        lightcurve_average['transit_out']['phase'] = \
            lightcurve_average['observations']['phase'][lightcurve_average['transit_out_flag']]

        for band_key in lightcurve_average['bands_list']:
            for name_append in append_list:

                lightcurve_average['observations']['delta_' + band_key + name_append] = \
                    np.asarray(lightcurve_average['observations']['delta_' + band_key + name_append])[sorting_index]
                lightcurve_average['transit_in']['delta_' + band_key + name_append] = \
                    lightcurve_average['observations']['delta_' + band_key + name_append][lightcurve_average['transit_in_flag']]
                lightcurve_average['transit_full']['delta_' + band_key + name_append] = \
                    lightcurve_average['observations']['delta_' + band_key + name_append][lightcurve_average['transit_full_flag']]
                lightcurve_average['transit_out']['delta_' + band_key + name_append] = \
                    lightcurve_average['observations']['delta_' + band_key + name_append][lightcurve_average['transit_out_flag']]


        pre_duration = transit_full_bins[0] - lightcurve_average['transit_out']['phase'][0]
        if pre_duration > 0:
            nsteps_pre = int(pre_duration / transit_full_step)
            if pre_duration % transit_full_step > 0.0:
                nsteps_pre += 1
        else:
            nsteps_pre = 0

        post_duration = lightcurve_average['transit_out']['phase'][-1] - transit_full_bins[-1]
        if post_duration > 0:
            nsteps_post = int(post_duration / transit_full_step)
            if post_duration % transit_full_step > 0.0:
                nsteps_post += 1
        else:
            nsteps_post = 0

        transit_bins = np.arange(transit_full_bins[0] - nsteps_pre * transit_full_step,
                                transit_full_bins[-1] + (nsteps_post + 1.1) * transit_full_step,
                                transit_full_step)

        lightcurve_average['binned'] = {
            'observations': {
                'phase': np.zeros(len(transit_bins)),
            },
            'transit_in': {},
            'transit_full': {},
            'transit_step': {},
            'transit_out': {},
        }
        for band_key in C_bands:
            for name_append in append_list:
                lightcurve_average['binned']['observations']['delta_' + band_key + name_append] = np.zeros([len(transit_bins), 2])

        transit_out_flag = np.zeros(len(transit_bins), dtype=bool)
        transit_in_flag = np.zeros(len(transit_bins), dtype=bool)
        transit_full_flag = np.zeros(len(transit_bins), dtype=bool)

        n_a = 0
        for nb in range(0, len(transit_bins) - 1):
            sel = (lightcurve_average['observations']['phase'] >= transit_bins[nb]) \
                & (lightcurve_average['observations']['phase'] < transit_bins[nb + 1])

            if np.sum(sel) <= 0: continue
            lightcurve_average['binned']['observations']['phase'][n_a] = np.average(
                lightcurve_average['observations']['phase'][sel])

            for band_key in C_bands:
                for name_append in append_list:

                    lightcurve_average['binned']['observations']['delta_' + band_key + name_append][n_a, 0], sum_weights = np.average(
                        lightcurve_average['observations']['delta_' + band_key + name_append][sel, 0],
                        weights=1. / lightcurve_average['observations']['delta_' + band_key + name_append][sel, 1] ** 2,
                        returned=True)

                    lightcurve_average['binned']['observations']['delta_' + band_key + name_append][n_a, 1] = np.sqrt(1. / sum_weights)


            if np.abs(lightcurve_average['binned']['observations']['phase'][n_a]) >= \
                total_transit_duration/2./planet_dict['period'][0]:
                    transit_out_flag[n_a] = True
            elif np.abs(lightcurve_average['binned']['observations']['phase'][n_a]) >= \
                full_transit_duration/2./planet_dict['period'][0]:
                    transit_in_flag[n_a] = True
            else:
                transit_full_flag[n_a] = True

            n_a += 1  # bins actually computed

        lightcurve_average['binned']['transit_in']['phase'] = lightcurve_average['binned']['observations']['phase'][transit_in_flag]
        lightcurve_average['binned']['transit_full']['phase'] = lightcurve_average['binned']['observations']['phase'][transit_full_flag]
        lightcurve_average['binned']['transit_out']['phase'] = lightcurve_average['binned']['observations']['phase'][transit_out_flag]
        lightcurve_average['binned']['observations']['phase'] = lightcurve_average['binned']['observations']['phase'][:n_a]

        for band_key in C_bands:
            for name_append in append_list:
                lightcurve_average['binned']['transit_in']['delta_' + band_key + name_append] = \
                    lightcurve_average['binned']['observations']['delta_' + band_key + name_append][transit_in_flag, :]
                lightcurve_average['binned']['transit_full']['delta_' + band_key + name_append] = \
                    lightcurve_average['binned']['observations']['delta_' + band_key + name_append][transit_full_flag, :]
                lightcurve_average['binned']['transit_out']['delta_' + band_key + name_append] = \
                    lightcurve_average['binned']['observations']['delta_' + band_key + name_append][transit_out_flag, :]
                lightcurve_average['binned']['observations']['delta_' + band_key + name_append] = \
                    lightcurve_average['binned']['observations']['delta_' + band_key + name_append][:n_a, :]


        save_to_cpickle(subroutine_name + '_'+reference+ '_' + output_selection, lightcurve_average, config_in['output'], lines=lines_label)


def plot_transmission_lightcurve_average(config_in, night_input='', reference='planetRF'):

    import matplotlib.pyplot as plt

    lightcurve_average = load_from_cpickle('transmission_lightcurve_average_'+reference, config_in['output'])

    for band_key in lightcurve_average['C_bands']:
        plt.figure(figsize=(12, 6))
        plt.title('Average transmission lightcurve\n {0:s}'.format(band_key))
        plt.errorbar(lightcurve_average['observations']['phase'],
                    lightcurve_average['observations']['delta_' + band_key][:,0],
                     yerr= lightcurve_average['observations']['delta_' + band_key][:,1] ,
                     fmt='.', c='k', alpha=0.25, label='observations')
        plt.errorbar(lightcurve_average['binned']['observations']['phase'],
                    lightcurve_average['binned']['observations']['delta_' + band_key][:, 0],
                     yerr= lightcurve_average['binned']['observations']['delta_' + band_key][:,1],
                     fmt='.', c='k', alpha=1.0, label='observations')

        plt.axvspan(-1, lightcurve_average['bins']['transit_in_bins'][0], alpha=0.25, color='green')
        plt.axvspan(lightcurve_average['bins']['transit_in_bins'][-1], 1., alpha=0.25, color='green')

        plt.axhline(0, c='C1')
        plt.xlim(lightcurve_average['observations']['phase'][0]-0.01,
                 lightcurve_average['observations']['phase'][-1]+0.01)
        plt.xlabel('$\lambda$ [$\AA$]')
        plt.ylabel('$\mathcal{R}$ - 1.')
        plt.legend()
        plt.show()

