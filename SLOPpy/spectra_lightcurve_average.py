from __future__ import print_function, division
from SLOPpy.subroutines.common import *
from SLOPpy.subroutines.spectral_subroutines import *
from SLOPpy.subroutines.io_subroutines import *
from SLOPpy.subroutines.fit_subroutines import *
from SLOPpy.subroutines.shortcuts import *
from SLOPpy.subroutines.rebin_subroutines import *
from astropy.convolution import convolve, Box1DKernel

__all__ = ['compute_spectra_lightcurve_average', 'plot_spectra_lightcurve_average',
           'compute_spectra_lightcurve_average_clv_rm_correction',
           'plot_spectra_lightcurve_average_clv_rm_correction']


def compute_spectra_lightcurve_average_clv_rm_correction(config_in, lines_label):
    compute_spectra_lightcurve_average(config_in, lines_label)


def plot_spectra_lightcurve_average_clv_rm_correction(config_in, night_input=''):
    plot_spectra_lightcurve_average(config_in, night_input)

subroutine_name = 'spectra_lightcurve_average'
pick_files = 'spectra_lightcurve'
sampler_name = 'emcee'

def compute_spectra_lightcurve_average(config_in, lines_label):

    night_dict = from_config_get_nights(config_in)
    instrument_dict = from_config_get_instrument(config_in)
    system_dict = from_config_get_system(config_in)
    planet_dict = from_config_get_planet(config_in)
    spectral_lines = from_config_get_spectral_lines(config_in)

    lines_dict = spectral_lines[lines_label] # from_config_get_transmission_lightcurve(config_in)


    shared_data = load_from_cpickle('shared', config_in['output'])

    results_list = ['user',
                    'mcmc_night_MED',
                    'mcmc_night_MAP',
                    'mcmc_global_MED',
                    'mcmc_global_MAP'
                    'user_uncorrected']

    append_list = ['', '_uncorrected', '_clv_model']

    for results_selection in results_list:
        skip_iteration = False

        try:
            lightcurve_average = load_from_cpickle(subroutine_name+ '_' + results_selection, config_in['output'], lines=lines_label)
            print("{0:45s}                         {1:s}".format(subroutine_name, 'Retrieved'))
            return
        except (FileNotFoundError, IOError):
            print("  No average transmission lightcurve found for case:{0:s}, computing now ".format(results_selection))
            print("{0:45s}                         {1:s}".format(subroutine_name, 'Computing'))
            print()


        """
        C stabds for central
        """
        C_bands = {}
        for passband_key, passband_val in lines_dict['passbands'].items():
            C_bands[passband_key] = {}
            for line_key, line_val in lines_dict['lines'].items():
                C_bands[passband_key][line_key] = (np.abs(shared_data['coadd']['wave'] - line_val)*2. <passband_val)

        """
        S stand for side
        """
        S_bands = {}
        for band_key, band_val in lines_dict['continuum'].items():
            S_bands[band_key] = (shared_data['coadd']['wave'] >= band_val[0]) & (shared_data['coadd']['wave'] <= band_val[1])


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

        lightcurve_average = {
            'subroutine': subroutine_name,
            'transit_in_flag': [],
            'transit_full_flag': [],
            'transit_out_flag': [],
            'transit_in': {},
            'transit_full': {},
            'transit_out': {},
            'observations': {
                'phase': []
            },
            'bands_list': [],
            'C_bands': C_bands,
            'S_bands': S_bands,
            'average': {},
            'bins': {
                'transit_in_bins': transit_in_bins,
                'transit_in_step': transit_in_step,
                'transit_full_bins': transit_full_bins,
                'transit_full_step': transit_full_step
            }
        }

        for band_key in C_bands:
            for name_append in append_list:
                lightcurve_average['observations']['ratio_' + band_key + name_append] = []
            lightcurve_average['bands_list'].extend([band_key])

        for night in night_dict:

            try:
                lightcurve = load_from_cpickle(pick_files + '_' + results_selection, config_in['output'], night, lines_label)
            except:
                print("  No night spectra lightcurve found for case:{0:s}, skipping ".format(results_selection))
                skip_iteration = True
                continue

            #lightcurve = load_from_cpickle(pick_files, config_in['output'], night)

            lightcurve_average['observations']['phase'].extend(lightcurve['arrays']['observations']['phase'].tolist())

            lightcurve_average['transit_in_flag'].extend(
                lightcurve['arrays']['observations']['transit_in_flag'].tolist())
            lightcurve_average['transit_full_flag'].extend(
                lightcurve['arrays']['observations']['transit_full_flag'].tolist())
            lightcurve_average['transit_out_flag'].extend(
                lightcurve['arrays']['observations']['transit_out_flag'].tolist())

            for band_key in lightcurve_average['bands_list']:
                for name_append in append_list:
                    lightcurve_average['observations']['ratio_' + band_key+ name_append].extend(
                        lightcurve['arrays']['observations']['ratio_' + band_key+ name_append].tolist())

        if skip_iteration: continue


        sorting_index = np.argsort(lightcurve_average['observations']['phase'])
        lightcurve_average['observations']['phase'] = np.asarray(lightcurve_average['observations']['phase'])[sorting_index]
        lightcurve_average['transit_in_flag'] = np.asarray(lightcurve_average['transit_in_flag'], dtype=np.int16)[sorting_index]
        lightcurve_average['transit_full_flag'] = np.asarray(lightcurve_average['transit_full_flag'], dtype=np.int16)[sorting_index]
        lightcurve_average['transit_out_flag'] = np.asarray(lightcurve_average['transit_out_flag'], dtype=np.int16)[sorting_index]

        lightcurve_average['transit_in']['phase'] = \
            lightcurve_average['observations']['phase'][lightcurve_average['transit_in_flag']]
        lightcurve_average['transit_full']['phase'] = \
            lightcurve_average['observations']['phase'][lightcurve_average['transit_full_flag']]
        lightcurve_average['transit_out']['phase'] = \
            lightcurve_average['observations']['phase'][lightcurve_average['transit_out_flag']]

        #for band_key in lightcurve_average['bands_list']:
        for band_key in C_bands:
            for name_append in append_list:
                lightcurve_average['observations']['ratio_' + band_key + name_append] = \
                    np.asarray(lightcurve_average['observations']['ratio_' + band_key + name_append])[sorting_index]
                lightcurve_average['transit_in']['ratio_' + band_key + name_append] = \
                    lightcurve_average['observations']['ratio_' + band_key + name_append][lightcurve_average['transit_in_flag']]
                lightcurve_average['transit_full']['ratio_' + band_key + name_append] = \
                    lightcurve_average['observations']['ratio_' + band_key + name_append][lightcurve_average['transit_full_flag']]
                lightcurve_average['transit_out']['ratio_' + band_key + name_append] = \
                    lightcurve_average['observations']['ratio_' + band_key + name_append][lightcurve_average['transit_out_flag']]

                avg_out, avg_out_sq = \
                    np.average(lightcurve_average['transit_out']['ratio_' + band_key + name_append][:, 0],
                            weights=1. / (lightcurve_average['transit_out']['ratio_' + band_key + name_append][:, 1]) ** 2,
                            returned=True)
                avg_in, avg_in_sq = \
                    np.average(lightcurve_average['transit_in']['ratio_' + band_key + name_append][:, 0],
                            weights=1. / (lightcurve_average['transit_in']['ratio_' + band_key + name_append][:, 1]) ** 2,
                            returned=True)
                avg_full, avg_full_sq = \
                    np.average(lightcurve_average['transit_full']['ratio_' + band_key + name_append][:, 0],
                            weights=1. / (lightcurve_average['transit_full']['ratio_' + band_key + name_append][:, 1]) ** 2,
                            returned=True)

                avg_out = \
                    np.average(lightcurve_average['transit_out']['ratio_' + band_key + name_append][:, 0])
                avg_in = \
                    np.average(lightcurve_average['transit_in']['ratio_' + band_key + name_append][:, 0])
                avg_full = \
                    np.average(lightcurve_average['transit_full']['ratio_' + band_key + name_append][:, 0])


                lightcurve_average['average'][band_key + name_append] = {
                    'average_out': np.asarray([avg_out, 1. / np.power(avg_out_sq, 0.5)]),
                    'average_in': np.asarray([avg_in, 1. / np.power(avg_in_sq, 0.5)]),
                    'average_full': np.asarray([avg_full, 1. / np.power(avg_full_sq, 0.5)]),
                }

                delta_fac = \
                    lightcurve_average['average'][band_key + name_append]['average_full'][0] / lightcurve_average['average'][band_key + name_append]['average_out'][0]
                delta_err = delta_fac * np.sqrt(
                    (lightcurve_average['average'][band_key + name_append]['average_out'][1]
                    / lightcurve_average['average'][band_key + name_append]['average_out'][0]) ** 2
                    + (lightcurve_average['average'][band_key + name_append]['average_full'][1]
                    / lightcurve_average['average'][band_key + name_append]['average_full'][0]) ** 2)

                lightcurve_average['average'][band_key + name_append]['delta'] = np.asarray([(1. - delta_fac) * 100., delta_err * 100.])


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
            'transit_out': {},
        }
        for band_key in C_bands:
            for name_append in append_list:
                lightcurve_average['binned']['observations']['ratio_' + band_key + name_append] = np.zeros([len(transit_bins), 2])

        transit_out_flag = np.zeros(len(transit_bins), dtype=bool)
        transit_full_flag = np.zeros(len(transit_bins), dtype=bool)
        transit_in_flag = np.zeros(len(transit_bins), dtype=bool)

        n_a = 0
        for nb in range(0, len(transit_bins) - 1):
            sel = (lightcurve_average['observations']['phase'] >= transit_bins[nb]) \
                & (lightcurve_average['observations']['phase'] < transit_bins[nb + 1])

            if np.sum(sel) <= 0: continue
            lightcurve_average['binned']['observations']['phase'][n_a] = np.average(
                lightcurve_average['observations']['phase'][sel])

            for band_key in C_bands:
                for name_append in append_list:
                    lightcurve_average['binned']['observations']['ratio_' + band_key + name_append][n_a, 0], sum_weights = np.average(
                        lightcurve_average['observations']['ratio_' + band_key + name_append][sel, 0],
                        weights=1. / lightcurve_average['observations']['ratio_' + band_key + name_append][sel, 1] ** 2,
                        returned=True)

                    lightcurve_average['binned']['observations']['ratio_' + band_key + name_append][n_a, 1] = np.sqrt(1. / sum_weights)

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
        lightcurve_average['binned']['transit_out']['phase'] = lightcurve_average['binned']['observations']['phase'][transit_out_flag]
        lightcurve_average['binned']['observations']['phase'] = lightcurve_average['binned']['observations']['phase'][:n_a]

        for band_key in C_bands:
            for name_append in append_list:
                lightcurve_average['binned']['transit_in']['ratio_' + band_key + name_append] = \
                    lightcurve_average['binned']['observations']['ratio_' + band_key + name_append][transit_in_flag, :]
                lightcurve_average['binned']['transit_out']['ratio_' + band_key + name_append] = \
                    lightcurve_average['binned']['observations']['ratio_' + band_key + name_append][transit_out_flag, :]
                lightcurve_average['binned']['observations']['ratio_' + band_key + name_append] = \
                    lightcurve_average['binned']['observations']['ratio_' + band_key + name_append][:n_a, :]

        save_to_cpickle(subroutine_name+ '_' + results_selection, lightcurve_average, config_in['output'], lines=lines_label)


def plot_spectra_lightcurve_average(config_in, night_input='', clv_rm_correction=False):

    import matplotlib.pyplot as plt

    if clv_rm_correction:
        subroutine_name = 'spectra_lightcurve_average_clv_rm_correction'
    else:
        subroutine_name = 'spectra_lightcurve_average'

    try:
        lightcurve_average = load_from_cpickle(subroutine_name, config_in['output'])
        print("{0:45s}                         {1:s}".format(subroutine_name, 'Plotting'))
    except:
        print("{0:45s}                         {1:s}".format(subroutine_name, 'Plot skipped'))
        return

    C_bands = lightcurve_average['C_bands']

    print()
    for band_key in C_bands:
        print("Average      Band: {1:s}   Delta:{2:8.4f} +- {3:8.4f} [%]".format('   ', band_key,
                                                                                lightcurve_average['average'][band_key][
                                                                                       'delta'][0],
                                                                                lightcurve_average['average'][band_key][
                                                                                       'delta'][1]))

    for band_key in C_bands:
        plt.figure(figsize=(12, 6))
        plt.title('Average spectra lightcurve\n {0:s}'.format(band_key))
        plt.errorbar(lightcurve_average['observations']['phase'],
                    lightcurve_average['observations']['ratio_' + band_key][:,0]*100 -100.,
                     yerr= lightcurve_average['observations']['ratio_' + band_key][:,1]*100 ,
                     fmt='.', c='k', alpha=0.25, label='observations')
        plt.errorbar(lightcurve_average['binned']['observations']['phase'],
                    lightcurve_average['binned']['observations']['ratio_' + band_key][:, 0]*100 -100.,
                     yerr= lightcurve_average['binned']['observations']['ratio_' + band_key][:,1]*100 ,
                     fmt='.', c='k', alpha=1.0, label='observations')

        plt.axvspan(-1, lightcurve_average['bins']['transit_in_bins'][0], alpha=0.25, color='green')
        plt.axvspan(lightcurve_average['bins']['transit_in_bins'][-1], 1., alpha=0.25, color='green')

        plt.axhline(0, c='C1')
        plt.xlim(lightcurve_average['observations']['phase'][0]-0.01,
                 lightcurve_average['observations']['phase'][-1]+0.01)
        plt.xlabel('orbital phase')
        plt.ylabel('$\mathcal{R}$ - 1.')
        plt.legend()
        plt.show()

