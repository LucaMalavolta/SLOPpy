from __future__ import print_function, division
from SLOPpy.subroutines.common import *
from SLOPpy.subroutines.spectral_subroutines import *
from SLOPpy.subroutines.io_subroutines import *
from SLOPpy.subroutines.fit_subroutines import *
from SLOPpy.subroutines.shortcuts import *
from SLOPpy.subroutines.rebin_subroutines import *
from astropy.convolution import convolve, Box1DKernel

__all__ = ['compute_transmission_lightcurve_planetRF',
           'plot_transmission_lightcurve_planetRF',
           'compute_transmission_lightcurve_stellarRF',
           'plot_transmission_lightcurve_stellarRF',
           'compute_transmission_lightcurve_observerRF',
           'plot_transmission_lightcurve_observerRF',
           'compute_transmission_lightcurve',
           'plot_transmission_lightcurve'
           ]


def compute_transmission_lightcurve_planetRF(config_in, lines_label):
    compute_transmission_lightcurve(config_in, lines_label, reference='planetRF')


def plot_transmission_lightcurve_planetRF(config_in, night_input):
    plot_transmission_lightcurve(config_in, night_input, reference='planetRF')


def compute_transmission_lightcurve_stellarRF(config_in, lines_label):
    compute_transmission_lightcurve(config_in, lines_label, reference='stellarRF')


def plot_transmission_lightcurve_stellarRF(config_in, night_input):
    plot_transmission_lightcurve(config_in, night_input, reference='stellarRF')


def compute_transmission_lightcurve_observerRF(config_in, lines_label):
    compute_transmission_lightcurve(config_in, lines_label, reference='observerRF')


def plot_transmission_lightcurve_observerRF(config_in, night_input):
    plot_transmission_lightcurve(config_in, night_input, reference='observerRF')

subroutine_name = 'transmission_lightcurve'
pick_files =  'transmission_spectrum'

def compute_transmission_lightcurve(config_in, lines_label, reference='planetRF'):

    do_average_instead_of_sum = True

    night_dict = from_config_get_nights(config_in)
    #instrument_dict = from_config_get_instrument(config_in)
    #system_dict = from_config_get_system(config_in)
    planet_dict = from_config_get_planet(config_in)

    spectral_lines = from_config_get_spectral_lines(config_in)
    #shared_data = load_from_cpickle('shared', config_in['output'])

    lines_dict = spectral_lines[lines_label] # from_config_get_transmission_lightcurve(config_in)

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

    output_list = ['user',
                    'mcmc_night_MED',
                    'mcmc_night_MAP',
                    'mcmc_global_MED',
                    'mcmc_global_MAP']

    append_list = ['', '_uncorrected', '_clv_model']

    for output_selection in output_list:
        skip_iteration = False

        for night in night_dict:

            print("compute_transmission_lightcurve            Night: ", night)

            try:
                transmission = load_from_cpickle(pick_files+'_'+reference + '_' + output_selection, config_in['output'], night, lines_label)
            except (FileNotFoundError, IOError):
                print('No transmission spectra found for case:{0:s}, be sure to run transmission_spectra before this step'.format(output_selection))
                skip_iteration = True

            if skip_iteration: continue

            try:
                lightcurve = load_from_cpickle(subroutine_name +'_'+reference+ '_' + output_selection, config_in['output'], night, lines_label)
                print()
                continue
            except (FileNotFoundError, IOError):
                print("  No transmission_lightcurve file found, computing now ")
                print()


            """ Retrieving the list of observations"""
            lists = load_from_cpickle('lists', config_in['output'], night)

            #calib_data = load_from_cpickle('calibration_fibA', config_in['output'], night)
            #input_data = retrieve_observations( config_in['output'], night, lists['observations'])
            observational_pams = load_from_cpickle('observational_pams', config_in['output'], night)

            # doublet sodium in the lab reference frame

            """
            C stands for central
            """
            C_bands = {}
            for passband_key, passband_val in lines_dict['passbands'].items():
                C_bands[passband_key] = {}
                for line_key, line_val in lines_dict['lines'].items():
                    C_bands[passband_key][line_key] = (np.abs(transmission['wave'] - line_val) * 2. < passband_val)

            """
            S stands for side
            """
            S_bands = {}
            for band_key, band_val in lines_dict['continuum'].items():
                S_bands[band_key] = (transmission['wave'] >= band_val[0]) & (transmission['wave'] <= band_val[1])


            processed = {
                'subroutine': subroutine_name
            }

            lightcurve = {
                'subroutine': subroutine_name,
                'arrays': {
                    'observations': {
                        'obs_name': np.zeros(len(lists['observations']), dtype=str),
                        'phase': np.zeros(len(lists['observations'])),
                    },
                    'transit_in': {},
                    'transit_full': {},
                    'transit_out': {},
                },
                'C_bands': C_bands,
                'S_bands': S_bands,
                'bins': {
                    'transit_in_bins': transit_in_bins,
                    'transit_in_step': transit_in_step,
                    'transit_full_bins': transit_full_bins,
                    'transit_full_step': transit_full_step
                }
            }
            """ Adding the C-bands arrays to the dictionary"""

            for band_key in C_bands:
                for name_append in append_list:
                    lightcurve['arrays']['observations']['delta_' + band_key + name_append] = np.zeros([len(lists['observations']), 2])

            transit_out_flag = np.zeros(len(lists['observations']), dtype=bool)
            transit_in_flag = np.zeros(len(lists['observations']), dtype=bool)
            transit_full_flag = np.zeros(len(lists['observations']), dtype=bool)

            for n_obs, obs in enumerate( lists['observations']):

                processed[obs] = {}
                lightcurve[obs] = {}


                try:
                    phase_internal = (observational_pams[obs]['BJD'] - night_dict[night]['time_of_transit'][0])/planet_dict['period'][0]
                except:
                    phase_internal = (observational_pams[obs]['BJD'] - night_dict[night]['time_of_transit'])/planet_dict['period'][0]


                processed[obs]['bands'] = {
                    'phase': phase_internal
                }
                processed[obs]['bands_uncorrected'] = {
                    'phase': phase_internal
                }
                processed[obs]['bands_clv_model'] = {
                    'phase': phase_internal
                }

                processed[obs]['s_integrated'] = 0.000
                processed[obs]['s_integrated_uncorrected'] = 0.000
                processed[obs]['s_integrated_clv_model'] = 0.000
                processed[obs]['s_sigmaq_sum'] = 0.000

                n_bands = 0.0
                for band_key, band_val in S_bands.items():
                    if do_average_instead_of_sum:
                        processed[obs]['bands'][band_key] =  \
                            [np.average(transmission[obs]['normalized'][band_val]),
                            np.sum((transmission[obs]['normalized_err'][band_val])**2)
                            /len(transmission[obs]['normalized_err'][band_val])**2]
                        processed[obs]['bands_uncorrected'][band_key] =  \
                            [np.average(transmission[obs]['normalized_uncorrected'][band_val]),
                            np.sum((transmission[obs]['normalized_uncorrected_err'][band_val])**2)
                            /len(transmission[obs]['normalized_uncorrected_err'][band_val])**2]
                        processed[obs]['bands_clv_model'][band_key] =  \
                            [np.average(transmission[obs]['clv_model_rebinned'][band_val]),
                            np.sum((transmission[obs]['normalized_err'][band_val])**2)
                            /len(transmission[obs]['normalized_err'][band_val])**2]
                    else:
                        processed[obs]['bands'][band_key] =  \
                            [np.sum(transmission[obs]['normalized'][band_val]),
                            np.sum((transmission[obs]['normalized_err'][band_val])**2)]
                        processed[obs]['bands_uncorrected'][band_key] =  \
                            [np.sum(transmission[obs]['normalized_uncorrected'][band_val]),
                            np.sum((transmission[obs]['normalized_uncorrected_err'][band_val])**2)]
                        processed[obs]['bands_clv_model'][band_key] =  \
                            [np.sum(transmission[obs]['clv_model_rebinned'][band_val]),
                            np.sum((transmission[obs]['normalized_err'][band_val])**2)]
                    processed[obs]['s_integrated'] += processed[obs]['bands'][band_key][0]
                    processed[obs]['s_integrated_uncorrected'] += processed[obs]['bands_uncorrected'][band_key][0]
                    processed[obs]['s_integrated_clv_model'] += processed[obs]['bands_clv_model'][band_key][0]
                    processed[obs]['s_sigmaq_sum'] += processed[obs]['bands'][band_key][1]
                    n_bands += 1.

                processed[obs]['s_integrated'] /= n_bands
                processed[obs]['s_integrated_uncorrected'] /= n_bands
                processed[obs]['s_integrated_clv_model'] /= n_bands
                processed[obs]['s_sigmaq_sum'] /= n_bands**2

                #processed[obs]['s_factor'] =np.power(s_integrated, -2.0)
                #processed[obs]['s_factor_clv_model'] =np.power(s_integrated, -2.0)
                #processed[obs]['s_factor_uncorrected'] = np.power(s_integrated, -2.0)

                for band_key, band_dict in C_bands.items():
                    processed[obs]['bands'][band_key] = {}
                    processed[obs]['bands_uncorrected'][band_key] = {}
                    processed[obs]['bands_clv_model'][band_key] = {}

                    processed[obs]['c_integrated'] = 0.000
                    processed[obs]['c_integrated_uncorrected'] = 0.000
                    processed[obs]['c_integrated_clv_model'] = 0.000
                    processed[obs]['c_sigmaq_sum'] = 0.000

                    n_bands = 0.0
                    for line_key, line_val in band_dict.items():
                        if do_average_instead_of_sum:
                            processed[obs]['bands'][band_key][line_key] = \
                                [np.average(transmission[obs]['normalized'][line_val]),
                                np.sum((transmission[obs]['normalized_err'][line_val])**2)
                                / len(transmission[obs]['normalized_err'][line_val])**2]
                            processed[obs]['bands_uncorrected'][band_key][line_key] = \
                                [np.average(transmission[obs]['normalized_uncorrected'][line_val]),
                                np.sum((transmission[obs]['normalized_uncorrected_err'][line_val])**2)
                                / len(transmission[obs]['normalized_uncorrected_err'][line_val])**2]
                            processed[obs]['bands_clv_model'][band_key][line_key] = \
                                [np.average(transmission[obs]['clv_model_rebinned'][line_val]),
                                np.sum((transmission[obs]['normalized_err'][line_val])**2)
                                / len(transmission[obs]['normalized_err'][line_val])**2]
                        else:
                            processed[obs]['bands'][band_key][line_key] = \
                                [np.sum(transmission[obs]['normalized'][line_val]),
                                np.sum((transmission[obs]['normalized_err'][line_val]) ** 2)]
                            processed[obs]['bands_uncorrected'][band_key][line_key] = \
                                [np.sum(transmission[obs]['normalized_uncorrected'][line_val]),
                                np.sum((transmission[obs]['normalized_uncorrected_err'][line_val]) ** 2)]
                            processed[obs]['bands_clv_model'][band_key][line_key] = \
                                [np.sum(transmission[obs]['clv_model_rebinned'][line_val]),
                                np.sum((transmission[obs]['normalized_err'][line_val]) ** 2)]

                        processed[obs]['c_integrated'] += processed[obs]['bands'][band_key][line_key][0]
                        processed[obs]['c_integrated_uncorrected'] += processed[obs]['bands_uncorrected'][band_key][line_key][0]
                        processed[obs]['c_integrated_clv_model'] += processed[obs]['bands_clv_model'][band_key][line_key][0]
                        processed[obs]['c_sigmaq_sum'] += processed[obs]['bands'][band_key][line_key][1]
                        n_bands += 1.

                    processed[obs]['c_integrated'] /= n_bands
                    processed[obs]['c_integrated_uncorrected'] /= n_bands
                    processed[obs]['c_integrated_clv_model'] /= n_bands
                    processed[obs]['c_sigmaq_sum'] /= n_bands ** 2

                    for name_append in append_list:

                        lightcurve[obs]['delta_' + band_key + name_append] = [processed[obs]['c_integrated' + name_append] 
                                                                            - processed[obs]['s_integrated' + name_append],
                                                                np.sqrt(processed[obs]['s_sigmaq_sum'] + processed[obs]['c_sigmaq_sum'])]

                        lightcurve['arrays']['observations']['delta_' + band_key + name_append][n_obs, :] = \
                            lightcurve[obs]['delta_' + band_key + name_append][:]

                lightcurve[obs]['phase'] = processed[obs]['bands']['phase']
                lightcurve['arrays']['observations']['obs_name'][n_obs] = obs
                lightcurve['arrays']['observations']['phase'][n_obs] = lightcurve[obs]['phase']

                if obs in lists['transit_out']:
                    transit_out_flag[n_obs] = True
                if obs in lists['transit_in']:
                    transit_in_flag[n_obs] = True
                if obs in lists['transit_full']:
                    transit_full_flag[n_obs] = True

            for band_key in C_bands:
                for name_append in append_list:
                    lightcurve['arrays']['rescaling_' + band_key + name_append] = \
                        np.average(lightcurve['arrays']['observations']['delta_' + band_key + name_append][transit_out_flag, 0], axis=0)

            sorting_index = np.argsort(lightcurve['arrays']['observations']['phase'])

            transit_out_flag = transit_out_flag[sorting_index]
            transit_in_flag = transit_in_flag[sorting_index]
            transit_full_flag = transit_full_flag[sorting_index]

            lightcurve['arrays']['observations']['obs_name'] = lightcurve['arrays']['observations']['obs_name'][sorting_index]
            lightcurve['arrays']['observations']['phase'] = lightcurve['arrays']['observations']['phase'][sorting_index]

            lightcurve['arrays']['transit_in']['obs_name'] = lightcurve['arrays']['observations']['obs_name'][transit_in_flag]
            lightcurve['arrays']['transit_in']['phase'] = lightcurve['arrays']['observations']['phase'][transit_in_flag]

            lightcurve['arrays']['transit_full']['obs_name'] = lightcurve['arrays']['observations']['obs_name'][transit_full_flag]
            lightcurve['arrays']['transit_full']['phase'] = lightcurve['arrays']['observations']['phase'][transit_full_flag]

            lightcurve['arrays']['transit_out']['obs_name'] = lightcurve['arrays']['observations']['obs_name'][transit_out_flag]
            lightcurve['arrays']['transit_out']['phase'] = lightcurve['arrays']['observations']['phase'][transit_out_flag]

            for band_key in C_bands:
                for name_append in append_list:

                    lightcurve['arrays']['observations']['delta_' + band_key + name_append] = \
                        lightcurve['arrays']['observations']['delta_' + band_key + name_append][sorting_index]
                        # / lightcurve['arrays']['rescaling_' + band_key]
                    lightcurve['arrays']['transit_in']['delta_' + band_key + name_append] = \
                        lightcurve['arrays']['observations']['delta_' + band_key + name_append][transit_in_flag]
                    lightcurve['arrays']['transit_full']['delta_' + band_key + name_append] = \
                        lightcurve['arrays']['observations']['delta_' + band_key + name_append][transit_full_flag]
                    lightcurve['arrays']['transit_out']['delta_' + band_key + name_append] = \
                        lightcurve['arrays']['observations']['delta_' + band_key + name_append][transit_out_flag]

            lightcurve['arrays']['observations']['transit_out_flag'] = transit_out_flag
            lightcurve['arrays']['observations']['transit_in_flag'] = transit_in_flag
            lightcurve['arrays']['observations']['transit_full_flag'] = transit_full_flag

            pre_duration = transit_full_bins[0] - lightcurve['arrays']['transit_out']['phase'][0]
            if pre_duration > 0:
                nsteps_pre = int(pre_duration/transit_full_step)
                if pre_duration % transit_full_step > 0.0:
                    nsteps_pre += 1
            else:
                nsteps_pre = 0

            post_duration = lightcurve['arrays']['transit_out']['phase'][-1] - transit_full_bins[-1]
            if post_duration > 0:
                nsteps_post = int(post_duration / transit_full_step)
                if post_duration % transit_full_step > 0.0:
                    nsteps_post += 1
            else:
                nsteps_post = 0

            transit_bins = np.arange(transit_full_bins[0]-nsteps_pre*transit_full_step,
                                    transit_full_bins[-1] + (nsteps_post+1.1) * transit_full_step,
                                    transit_full_step)

            lightcurve['binned'] = {
                'observations': {
                    'phase': np.zeros(len(transit_bins)),
                },
                'transit_in': {},
                'transit_full': {},
                'transit_out': {},
            }
            for band_key in C_bands:
                for name_append in append_list:
                    lightcurve['binned']['observations']['delta_' + band_key + name_append] = np.zeros([len(transit_bins), 2])

            transit_out_flag = np.zeros(len(transit_bins), dtype=bool)
            transit_in_flag = np.zeros(len(transit_bins), dtype=bool)
            transit_full_flag = np.zeros(len(transit_bins), dtype=bool)

            n_a = 0
            for nb in range(0, len(transit_bins)-1):
                sel = (lightcurve['arrays']['observations']['phase'] >= transit_bins[nb]) \
                    & (lightcurve['arrays']['observations']['phase'] < transit_bins[nb+1])

                if np.sum(sel) <= 0: continue
                lightcurve['binned']['observations']['phase'][n_a] = np.average(lightcurve['arrays']['observations']['phase'][sel])

                for band_key in C_bands:
                    for name_append in append_list:
                        lightcurve['binned']['observations']['delta_' + band_key + name_append][n_a, 0], sum_weights = np.average(
                            lightcurve['arrays']['observations']['delta_' + band_key + name_append][sel, 0],
                            weights=1. / lightcurve['arrays']['observations']['delta_' + band_key + name_append][sel, 1]**2,
                            returned=True)

                        lightcurve['binned']['observations']['delta_' + band_key + name_append][n_a, 1] = np.sqrt(1. / sum_weights)

                if np.abs(lightcurve['binned']['observations']['phase'][n_a]) >= \
                    total_transit_duration/2./planet_dict['period'][0]:
                        transit_out_flag[n_a] = True
                elif np.abs(lightcurve['binned']['observations']['phase'][n_a]) >= \
                    full_transit_duration/2./planet_dict['period'][0]:
                        transit_in_flag[n_a] = True
                else:
                    transit_full_flag[n_a] = True

                n_a += 1 # bins actually computed

            lightcurve['binned']['transit_in']['phase'] = lightcurve['binned']['observations']['phase'][transit_in_flag]
            lightcurve['binned']['transit_full']['phase'] = lightcurve['binned']['observations']['phase'][transit_full_flag]
            lightcurve['binned']['transit_out']['phase'] = lightcurve['binned']['observations']['phase'][transit_out_flag]
            lightcurve['binned']['observations']['phase'] = lightcurve['binned']['observations']['phase'][:n_a]

            for band_key in C_bands:
                for name_append in append_list:
                    lightcurve['binned']['transit_in']['delta_' + band_key + name_append] = \
                        lightcurve['binned']['observations']['delta_' + band_key + name_append][transit_in_flag, :]
                    lightcurve['binned']['transit_full']['delta_' + band_key + name_append] = \
                        lightcurve['binned']['observations']['delta_' + band_key + name_append][transit_full_flag, :]
                    lightcurve['binned']['transit_out']['delta_' + band_key + name_append] = \
                        lightcurve['binned']['observations']['delta_' + band_key + name_append][transit_out_flag, :]
                    lightcurve['binned']['observations']['delta_' + band_key + name_append] = \
                        lightcurve['binned']['observations']['delta_' + band_key + name_append][:n_a, :]

            save_to_cpickle(subroutine_name + '_'+reference + '_' + output_selection +'_processed', processed, config_in['output'], night, lines_label)
            save_to_cpickle(subroutine_name + '_'+reference + '_' + output_selection, lightcurve, config_in['output'], night, lines_label)


def plot_transmission_lightcurve(config_in, night_input='', reference='planetRF'):
    import matplotlib.pyplot as plt

    night_dict = from_config_get_nights(config_in)

    if night_input=='':
        night_list = night_dict
    else:
        night_list = np.atleast_1d(night_input)

    for night in night_list:

        print("plot_transmission_lightcurve               Night: ", night)

        """ Retrieving the analysis"""
        try:
            lightcurve = load_from_cpickle(subroutine_name+'_'+reference, config_in['output'], night)
        except:
            print()
            print("No transmission lightcurve dataset, no plots")
            continue

        observational_pams = load_from_cpickle('observational_pams', config_in['output'], night)
        C_bands = lightcurve['C_bands']

        for band_key in C_bands:
            plt.figure(figsize=(12, 6))
            plt.title('Transmission lightcurve - night {0:s} \n {1:s}'.format(night, band_key))
            plt.errorbar(lightcurve['arrays']['observations']['phase'],
                        lightcurve['arrays']['observations']['delta_' + band_key][:,0],
                         yerr= lightcurve['arrays']['observations']['delta_' + band_key][:,1],
                         fmt='.', c='k', alpha=0.25, label='observations')
            plt.errorbar(lightcurve['binned']['observations']['phase'],
                        lightcurve['binned']['observations']['delta_' + band_key][:, 0],
                         yerr= lightcurve['binned']['observations']['delta_' + band_key][:,1],
                         fmt='.', c='k', alpha=1.0, label='observations')

            plt.axvspan(-1, lightcurve['bins']['transit_in_bins'][0], alpha=0.25, color='green')
            plt.axvspan(lightcurve['bins']['transit_in_bins'][-1], 1., alpha=0.25, color='green')

            plt.axhline(0, c='C1')
            plt.xlim(lightcurve['arrays']['observations']['phase'][0]-0.01,
                     lightcurve['arrays']['observations']['phase'][-1]+0.01)
            plt.xlabel('$\lambda$ [$\AA$]')
            plt.ylabel('$\mathcal{R}$ - 1.')
            plt.legend()
            plt.show()




