from __future__ import print_function, division
from SLOPpy.subroutines.common import *
from SLOPpy.subroutines.spectral_subroutines import *
from SLOPpy.subroutines.io_subroutines import *
from SLOPpy.subroutines.fit_subroutines import *
from SLOPpy.subroutines.shortcuts import *
from SLOPpy.subroutines.rebin_subroutines import *
from SLOPpy.subroutines.clv_rm_subroutines import *
from SLOPpy.subroutines.math_functions import *
from astropy.convolution import convolve, Box1DKernel

__all__ = ['compute_spectra_lightcurve',
           'compute_spectra_lightcurve_clv_rm_correction',
           'plot_spectra_lightcurve',
           'plot_spectra_lightcurve_clv_rm_correction']


def compute_spectra_lightcurve_clv_rm_correction(config_in, lines_label):
    compute_spectra_lightcurve(config_in, lines_label)


def plot_spectra_lightcurve_clv_rm_correction(config_in, night_input=''):
    plot_spectra_lightcurve(config_in, night_input)

subroutine_name = 'spectra_lightcurve'
sampler_name = 'emcee'

def compute_spectra_lightcurve(config_in, lines_label):


    results_list_default = ['user',
                    'mcmc_night_MED',
                    'mcmc_night_MAP',
                    'mcmc_global_MED',
                    'mcmc_global_MAP']

    append_list = ['', '_uncorrected', '_clv_model']

    do_average_instead_of_sum = True

    night_dict = from_config_get_nights(config_in)
    #instrument_dict = from_config_get_instrument(config_in)
    #system_dict = from_config_get_system(config_in)
    planet_dict = from_config_get_planet(config_in)

    spectral_lines = from_config_get_spectral_lines(config_in)
    lines_dict = spectral_lines[lines_label]
    clv_rm_correction = lines_dict.get('clv_rm_correction', True)


    # from_config_get_transmission_lightcurve(config_in)
    #lightcurve_dict = from_config_get_transmission_lightcurve(config_in)

    shared_data = load_from_cpickle('shared', config_in['output'])

    """ Using the MCMC fit range to define the transmission spectrum region """

    shared_selection = (shared_data['coadd']['wave'] >= lines_dict['range'][0]) \
        & (shared_data['coadd']['wave'] < lines_dict['range'][1])

    preparation_template = {
        'subroutine': subroutine_name,
        'range': lines_dict['range'],
        'wave': shared_data['coadd']['wave'][shared_selection],
        'step': shared_data['coadd']['step'][shared_selection],
        'size': int(np.sum(shared_selection)),
    }

    if 'full_transit_duration' in planet_dict:
        full_transit_duration = planet_dict['total_transit_duration'][0]
    else:
        full_transit_duration = planet_dict['transit_duration'][0]

    if 'total_transit_duration' in planet_dict:
        total_transit_duration = planet_dict['total_transit_duration'][0]
    else:
        total_transit_duration = planet_dict['transit_duration'][0]

    """
        The transit phase [0-1] is divided in N (=5) bins. Two arrays are computed:
        - transit_in_bins: array with the boundaries of the bins, size=N+1
        - transit_in_step: average size of the bin, size=1
    """

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


    """ Preparation stage  - rebinning of spectra """

    for night in night_dict:

        preparation = None # Free up memory
        try:
            preparation = load_from_cpickle(subroutine_name + '_preparation', config_in['output'], night, lines_label)
            print("{0:45s} Night:{1:15s}   {2:s}".format(subroutine_name+ '_preparation', night, 'Retrieved'))
            continue
        except:
            print("{0:45s} Night:{1:15s}   {2:s}".format(subroutine_name+ '_preparation', night, 'Computing'))
            print()

        preparation = preparation_template.copy()

        """ Retrieving the list of observations"""
        lists = load_from_cpickle('lists', config_in['output'], night)

        calib_data = load_from_cpickle('calibration_fibA', config_in['output'], night)
        input_data = retrieve_observations( config_in['output'], night, lists['observations'])
        observational_pams = load_from_cpickle('observational_pams', config_in['output'], night)


        for n_obs, obs in enumerate( lists['observations']):

            preparation[obs] = {}

            preparation[obs]['rescaling'], \
            preparation[obs]['rescaled'], \
            preparation[obs]['rescaled_err'] = perform_rescaling(
                input_data[obs]['wave'], input_data[obs]['e2ds'], input_data[obs]['e2ds_err'],
                observational_pams['wavelength_rescaling'])

            preserve_flux = input_data[obs].get('absolute_flux', True)

            preparation[obs]['rebinned'] = \
                rebin_2d_to_1d(input_data[obs]['wave'],
                               input_data[obs]['step'],
                               input_data[obs]['e2ds'],
                               calib_data['blaze'],
                               preparation['wave'],
                               preparation['step'],
                               preserve_flux=preserve_flux,
                               rv_shift=observational_pams[obs]['rv_shift_ORF2SRF'])

            preparation[obs]['rebinned_err'] = \
                rebin_2d_to_1d(input_data[obs]['wave'],
                               input_data[obs]['step'],
                               input_data[obs]['e2ds_err'],
                               calib_data['blaze'],
                               preparation['wave'],
                               preparation['step'],
                               preserve_flux=preserve_flux,
                               rv_shift=observational_pams[obs]['rv_shift_ORF2SRF'])


        save_to_cpickle(subroutine_name+'_preparation', preparation, config_in['output'], night, lines_label)

        # Free up memory
        calib_data = None
        input_data = None
        observational_pams = None

    """ Actual computation of spectral lightcurve """
    # doublet sodium in the lab reference frame

    """
    C stands for central
    """
    C_bands = {}
    for passband_key, passband_val in lines_dict['passbands'].items():
        C_bands[passband_key] = {}
        for line_key, line_val in lines_dict['lines'].items():
            C_bands[passband_key][line_key] = (np.abs(preparation['wave'] - line_val) < passband_val / 2.)

    """
    S stands for side
    """
    S_bands = {}
    for band_key, band_val in lines_dict['continuum'].items():
        S_bands[band_key] = (preparation['wave'] >= band_val[0]) & (preparation['wave'] <= band_val[1])



    results_list = results_list_default.copy()

    for results_selection in results_list_default:
        skip_iteration = False

        for night in night_dict:

            print_warning = True
            if skip_iteration: continue

            
            binned_mcmc_night = check_existence_cpickle(
                    'transmission_binned_mcmc_'+sampler_name+'_results', config_in['output'], night, lines_label)
            binned_mcmc_global = check_existence_cpickle(
                    'transmission_binned_mcmc_'+sampler_name+'_results', config_in['output'], lines_label)

            mcmc_night = check_existence_cpickle(
                    'transmission_mcmc_'+sampler_name+'_results', config_in['output'], night, lines_label)
            mcmc_global = check_existence_cpickle(
                    'transmission_mcmc_'+sampler_name+'_results', config_in['output'], lines_label)

            results_list = ['user']
            if (mcmc_night or binned_mcmc_night):
                results_list.append(['mcmc_night_MED', 'mcmc_night_MAP'])
            if (mcmc_global or binned_mcmc_global):
                results_list.append(['mcmc_global_MED', 'mcmc_global_MAP'])
            if results_selection not in results_list:
                print('  {0:s} results not found, skipping iteration'.format(results_selection))
                skip_iteration = True
                continue


            if mcmc_night and results_selection in  ['mcmc_night_MED', 'mcmc_night_MAP']:
                mcmc_results_night = load_from_cpickle(
                    'transmission_mcmc_'+sampler_name+'_results', config_in['output'], night, lines_label)

                print('  Observational parameters from MCMC fit of unbinned data, individual night')

            elif mcmc_global and results_selection in  ['mcmc_global_MED', 'mcmc_global_MAP']:
                mcmc_results_global = load_from_cpickle(
                    'transmission_mcmc_'+sampler_name+'_results', config_in['output'], lines=lines_label)

                print('  Observational parameters from MCMC fit of unbinned data, global fit')

            elif binned_mcmc_night and results_selection in  ['mcmc_night_MED', 'mcmc_night_MAP']:
                mcmc_results_night = load_from_cpickle(
                    'transmission_binned_mcmc_'+sampler_name+'_results', config_in['output'], night, lines_label)

                print('  Observational parameters from MCMC fit of binned data, individual night')

            elif binned_mcmc_global and results_selection in  ['mcmc_global_MED', 'mcmc_global_MAP']:
                mcmc_results_global = load_from_cpickle(
                    'transmission_binned_mcmc_'+sampler_name+'_results', config_in['output'], lines=lines_label)

                print('  Observational parameters from MCMC fit of binned data, global fit')
            else:
                print('  Observational parameters from configuration file')

            try:
                lightcurve = load_from_cpickle(subroutine_name + '_' + results_selection , config_in['output'], night, lines_label)
                print("{0:45s} Night:{1:15s}   {2:s}".format(subroutine_name, night, 'Retrieved'))
                continue
            except:
                print("{0:45s} Night:{1:15s}   {2:s}".format(subroutine_name, night, 'Computing'))
                print()

            preparation = load_from_cpickle(subroutine_name + '_preparation', config_in['output'], night, lines_label)

            """ Retrieving the list of observations"""
            lists = load_from_cpickle('lists', config_in['output'], night)
            observational_pams = load_from_cpickle('observational_pams', config_in['output'], night)

            processed = {
                'subroutine': 'compute_spectra_lightcurve',
                'range': preparation['range'],
                'wave': preparation['wave'],
                'step': preparation['step'],
                'size': preparation['size']
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
                'average': {},
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
                    lightcurve['arrays']['observations']['ratio_' + band_key + name_append] = np.zeros([len(lists['observations']), 2])

            transit_out_flag = np.zeros(len(lists['observations']), dtype=bool)
            transit_in_flag = np.zeros(len(lists['observations']), dtype=bool)
            transit_full_flag = np.zeros(len(lists['observations']), dtype=bool)

            if clv_rm_correction:
                try:
                    clv_rm_models = load_from_cpickle('clv_rm_models', config_in['output'], night, lines_label)
                except (FileNotFoundError, IOError):
                    clv_rm_models = load_from_cpickle('clv_rm_models', config_in['output'], night)



            """ Shift into planetary reference system is the default
            choice"""
            if results_selection == 'user':
                planet_R_factor = observational_pams.get('Rp_factor', 1.00000)
            elif results_selection == 'mcmc_night_MED':
                planet_R_factor = mcmc_results_night['results']['planet_R']
            elif results_selection == 'mcmc_night_MAP':
                planet_R_factor = mcmc_results_night['results_MAP']['planet_R']
            elif results_selection == 'mcmc_global_MED':
                planet_R_factor = mcmc_results_global['results']['planet_R']
            elif results_selection == 'mcmc_global_MAP':
                planet_R_factor = mcmc_results_global['results_MAP']['planet_R']

            for n_obs, obs in enumerate( lists['observations']):

                processed[obs] = {}
                lightcurve[obs] = {}


                processed[obs]['uncorrected'] = preparation[obs]['rebinned']
                processed[obs]['uncorrected_err'] = preparation[obs]['rebinned_err']

                if clv_rm_correction:

                    """" CLV + RM computation in the planetary reference frame """
                    processed[obs]['clv_model_stellarRF'] = interpolate1d_grid_nocheck(planet_R_factor,
                                                                                        clv_rm_models['common']['radius_grid'],
                                                                                        clv_rm_models[obs]['clv_rm_model_convolved_normalized'])

                    processed[obs]['clv_model_rebinned'] = \
                        rebin_1d_to_1d(clv_rm_models['common']['wave'],
                                        clv_rm_models['common']['step'],
                                        processed[obs]['clv_model_stellarRF'],
                                        processed['wave'],
                                        processed['step'],
                                        preserve_flux=False)

                    processed[obs]['rebinned'] = processed[obs]['uncorrected'] / processed[obs]['clv_model_rebinned']
                    processed[obs]['rebinned_err'] = processed[obs]['uncorrected_err'] / processed[obs]['clv_model_rebinned']

                else:
                    processed[obs]['clv_model_rebinned'] = np.ones(processed['size'])
                    processed[obs]['rebinned'] = processed[obs]['uncorrected']
                    processed[obs]['rebinned_err'] = processed[obs]['uncorrected_err']
                    if print_warning:
                        print('   *** No CLV correction')
                        print_warning = False


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

                n_bands = 0.00
                for band_key, band_val in S_bands.items():
                    if do_average_instead_of_sum:
                        processed[obs]['bands'][band_key] =  \
                            [np.average(processed[obs]['rebinned'][band_val]),
                            np.sum((processed[obs]['rebinned_err'][band_val])**2)
                            /len(processed[obs]['rebinned_err'][band_val])**2]
                        processed[obs]['bands_uncorrected'][band_key] =  \
                            [np.average(processed[obs]['uncorrected'][band_val]),
                            np.sum((processed[obs]['uncorrected_err'][band_val])**2)
                            /len(processed[obs]['uncorrected_err'][band_val])**2]
                        processed[obs]['bands_clv_model'][band_key] =  \
                            [np.average(processed[obs]['clv_model_rebinned'][band_val]),
                            np.sum((processed[obs]['rebinned_err'][band_val])**2)
                            /len(processed[obs]['rebinned_err'][band_val])**2]
                    else:
                        processed[obs]['bands'][band_key] =  \
                            [np.sum(processed[obs]['rebinned'][band_val]),
                            np.sum((processed[obs]['rebinned_err'][band_val])**2)]
                        processed[obs]['bands_uncorrected'][band_key] =  \
                            [np.sum(processed[obs]['uncorrected'][band_val]),
                            np.sum((processed[obs]['uncorrected_err'][band_val])**2)]
                        processed[obs]['bands_clv_model'][band_key] =  \
                            [np.sum(processed[obs]['clv_model_rebinned'][band_val]),
                            np.sum((processed[obs]['rebinned_err'][band_val])**2)]

                    processed[obs]['s_integrated'] += processed[obs]['bands'][band_key][0]
                    processed[obs]['s_integrated_uncorrected'] += processed[obs]['bands_uncorrected'][band_key][0]
                    processed[obs]['s_integrated_clv_model'] += processed[obs]['bands_clv_model'][band_key][0]
                    processed[obs]['s_sigmaq_sum'] += processed[obs]['bands'][band_key][1]
                    n_bands += 1.

                #todo:  why a 2 denominator??? 
                processed[obs]['s_integrated'] /= (n_bands / 2.)
                processed[obs]['s_integrated_uncorrected'] /= (n_bands / 2.)
                processed[obs]['s_integrated_clv_model'] /= (n_bands / 2.)
                processed[obs]['s_sigmaq_sum'] /= (n_bands / 2.)**2


                for band_key, band_dict in C_bands.items():
                    processed[obs]['bands'][band_key] = {}
                    processed[obs]['bands_uncorrected'][band_key] = {}
                    processed[obs]['bands_clv_model'][band_key] = {}

                    processed[obs]['c_integrated'] = 0.000
                    processed[obs]['c_integrated_uncorrected'] = 0.000
                    processed[obs]['c_integrated_clv_model'] = 0.000
                    processed[obs]['c_sigmaq_sum'] = 0.000

                    n_bands = 0.00

                    for line_key, line_val in band_dict.items():
                        if do_average_instead_of_sum:
                            processed[obs]['bands'][band_key][line_key] = \
                                [np.average(processed[obs]['rebinned'][line_val]),
                                np.sum((processed[obs]['rebinned_err'][line_val]) ** 2)
                                / len(processed[obs]['rebinned_err'][line_val]) ** 2]
                            processed[obs]['bands_uncorrected'][band_key][line_key] = \
                                [np.average(processed[obs]['uncorrected'][line_val]),
                                np.sum((processed[obs]['uncorrected_err'][line_val]) ** 2)
                                / len(processed[obs]['rebinned_err'][line_val]) ** 2]
                            processed[obs]['bands_clv_model'][band_key][line_key] = \
                                [np.average(processed[obs]['clv_model_rebinned'][line_val]),
                                np.sum((processed[obs]['rebinned_err'][line_val]) ** 2)
                                / len(processed[obs]['rebinned_err'][line_val]) ** 2]
                        else:
                            processed[obs]['bands'][band_key][line_key] = \
                                [np.sum(processed[obs]['rebinned'][line_val]),
                                np.sum((processed[obs]['rebinned_err'][line_val]) ** 2)]

                        processed[obs]['c_integrated'] += processed[obs]['bands'][band_key][line_key][0]
                        processed[obs]['c_integrated_uncorrected'] += processed[obs]['bands_uncorrected'][band_key][line_key][0]
                        processed[obs]['c_integrated_clv_model'] += processed[obs]['bands_clv_model'][band_key][line_key][0]
                        processed[obs]['c_sigmaq_sum'] += processed[obs]['bands'][band_key][line_key][1]
                        n_bands += 1.



                    processed[obs]['c_integrated'] /= (n_bands / 2.)
                    processed[obs]['c_integrated_uncorrected'] /= (n_bands / 2.)
                    processed[obs]['c_integrated_clv_model'] /= (n_bands / 2.)
                    processed[obs]['c_sigmaq_sum'] /= (n_bands / 2.) ** 2


                    for name_append in append_list:

                        ratio = processed[obs]['c_integrated' + name_append] / processed[obs]['s_integrated' + name_append]
                        ratio_err = ratio * np.sqrt(
                            processed[obs]['c_sigmaq_sum'] / processed[obs]['c_integrated' + name_append] ** 2
                            + processed[obs]['s_sigmaq_sum'] /  processed[obs]['s_integrated' + name_append] ** 2)

                        lightcurve[obs]['ratio_' + band_key + name_append] = [ratio, ratio_err]

                        lightcurve['arrays']['observations']['ratio_' + band_key + name_append][n_obs, :] = \
                            lightcurve[obs]['ratio_' + band_key + name_append][:]

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
                        np.average(lightcurve['arrays']['observations']['ratio_' + band_key + name_append][transit_out_flag, 0], axis=0)

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

                    lightcurve['arrays']['observations']['ratio_' + band_key + name_append] = \
                        lightcurve['arrays']['observations']['ratio_' + band_key + name_append][sorting_index] \
                        / lightcurve['arrays']['rescaling_' + band_key + name_append]
                    lightcurve['arrays']['transit_in']['ratio_' + band_key + name_append] = \
                        lightcurve['arrays']['observations']['ratio_' + band_key + name_append][transit_in_flag]
                    lightcurve['arrays']['transit_full']['ratio_' + band_key + name_append] = \
                        lightcurve['arrays']['observations']['ratio_' + band_key + name_append][transit_full_flag]
                    lightcurve['arrays']['transit_out']['ratio_' + band_key + name_append] = \
                        lightcurve['arrays']['observations']['ratio_' + band_key + name_append][transit_out_flag]

                    avg_out, avg_out_sq = \
                        np.average(lightcurve['arrays']['transit_out']['ratio_' + band_key + name_append][:, 0],
                                weights=1./(lightcurve['arrays']['transit_out']['ratio_' + band_key + name_append][:, 1])**2,
                                returned=True)
                    avg_in, avg_in_sq = \
                        np.average(lightcurve['arrays']['transit_in']['ratio_' + band_key + name_append][:, 0],
                                weights=1. / (lightcurve['arrays']['transit_in']['ratio_' + band_key + name_append][:, 1]) ** 2,
                                returned=True)
                    avg_full, avg_full_sq = \
                        np.average(lightcurve['arrays']['transit_full']['ratio_' + band_key + name_append][:, 0],
                                weights=1. / (lightcurve['arrays']['transit_full']['ratio_' + band_key + name_append][:, 1]) ** 2,
                                returned=True)

                    lightcurve['average'][band_key + name_append] = {
                        'average_out': np.asarray([avg_out, 1./np.power(avg_out_sq, 0.5)]),
                        'average_in': np.asarray([avg_in, 1. / np.power(avg_in_sq, 0.5)]),
                        'average_full': np.asarray([avg_full, 1. / np.power(avg_full_sq, 0.5)]),
                    }

                    delta_fac = (lightcurve['average'][band_key + name_append]['average_full'][0]
                                 / lightcurve['average'][band_key + name_append]['average_out'][0])

                    delta_err = delta_fac * np.sqrt(
                        (lightcurve['average'][band_key + name_append]['average_out'][1]
                        / lightcurve['average'][band_key + name_append]['average_out'][0]) ** 2
                        + (lightcurve['average'][band_key + name_append]['average_full'][1]
                        / lightcurve['average'][band_key + name_append]['average_full'][0]) ** 2)

                    lightcurve['average'][band_key + name_append]['delta'] = np.asarray([(1.-delta_fac)*100., delta_err*100.])

            lightcurve['arrays']['observations']['transit_out_flag'] = transit_out_flag
            lightcurve['arrays']['observations']['transit_in_flag'] = transit_in_flag
            lightcurve['arrays']['observations']['transit_full_flag'] = transit_full_flag

            """ Compute the duration of the pre-transit observations, using as scale
                the number of bins, with the same size as those used inside the
                transit.

                The value is given by the difference of the phase of the beginning of the transit minus
                the phase of the first observation, keeping in mind that the centre of the transit has phase = 0

                An additional bin is added if there are observations left out from the actual number of bins
            """
            pre_duration = transit_full_bins[0] - lightcurve['arrays']['transit_out']['phase'][0]
            if pre_duration > 0:
                nsteps_pre = int(pre_duration/transit_full_step)
                if pre_duration % transit_full_step > 0.0:
                    nsteps_pre += 1
            else:
                nsteps_pre = 0

            """ same as pre-transit, but suing the post-transit instead"""
            post_duration = lightcurve['arrays']['transit_out']['phase'][-1] - transit_full_bins[-1]
            if post_duration > 0:
                nsteps_post = int(post_duration / transit_full_step)
                if post_duration % transit_full_step > 0.0:
                    nsteps_post += 1
            else:
                nsteps_post = 0

            """ THe full array with both in-transit and out-transit phase, built in such a way that the 
                - the lower boundary of the first in-transit bin corresponds to the beginning of the transit
                - the upper boundary of the last in-transit bin corresponds to the end of the transit
            """
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
                    lightcurve['binned']['observations']['ratio_' + band_key + name_append] = np.zeros([len(transit_bins), 2])

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

                        lightcurve['binned']['observations']['ratio_' + band_key + name_append][n_a, 0], sum_weights = np.average(
                            lightcurve['arrays']['observations']['ratio_' + band_key + name_append][sel, 0],
                            weights=1. / lightcurve['arrays']['observations']['ratio_' + band_key + name_append][sel, 1]**2,
                            returned=True)

                        lightcurve['binned']['observations']['ratio_' + band_key + name_append][n_a, 1] = np.sqrt(1. / sum_weights)

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
                    lightcurve['binned']['transit_in']['ratio_' + band_key + name_append] = \
                        lightcurve['binned']['observations']['ratio_' + band_key + name_append][transit_in_flag, :]
                    lightcurve['binned']['transit_full']['ratio_' + band_key + name_append] = \
                        lightcurve['binned']['observations']['ratio_' + band_key + name_append][transit_full_flag, :]
                    lightcurve['binned']['transit_out']['ratio_' + band_key + name_append] = \
                        lightcurve['binned']['observations']['ratio_' + band_key + name_append][transit_out_flag, :]
                    lightcurve['binned']['observations']['ratio_' + band_key + name_append] = \
                        lightcurve['binned']['observations']['ratio_' + band_key + name_append][:n_a, :]

            save_to_cpickle(subroutine_name+ '_' + results_selection + '_processed', processed, config_in['output'], night, lines_label)
            save_to_cpickle(subroutine_name+ '_' + results_selection, lightcurve, config_in['output'], night, lines_label)

            # Forcing memory deallocation
            lightcurve = None
            processed = None
            preparation = None
            clv_rm_models = None

def plot_spectra_lightcurve(config_in, night_input='', clv_rm_correction=False):
    import matplotlib.pyplot as plt

    if clv_rm_correction:
        subroutine_name = 'spectra_lightcurve_clv_rm_correction'
    else:
        subroutine_name = 'spectra_lightcurve'

    night_dict = from_config_get_nights(config_in)

    if night_input=='':
        night_list = night_dict
    else:
        night_list = np.atleast_1d(night_input)

    for night in night_list:

        """ Retrieving the analysis"""
        try:
            lightcurve = load_from_cpickle(subroutine_name, config_in['output'], night)
            print("{0:45s} Night:{1:15s}   {2:s}".format(subroutine_name, night, 'Plotting'))
        except:
            print("{0:45s} Night:{1:15s}   {2:s}".format(subroutine_name, night, 'Skipped'))
            continue

        #observational_pams = load_from_cpickle('observational_pams', config_in['output'], night)
        C_bands = lightcurve['C_bands']


        print()
        for band_key in C_bands:
            print("Night: {0:s}   Band: {1:s}   Delta:{2:8.4f} +- {3:8.4f} [%]".format(night, band_key,
                  lightcurve['average'][band_key]['delta'][0], lightcurve['average'][band_key]['delta'][1]))

        for band_key in C_bands:


            plt.figure(figsize=(12, 6))
            plt.title('Spectra lightcurve - night {0:s} \n {1:s}'.format(night, band_key))
            plt.errorbar(lightcurve['arrays']['observations']['phase'],
                        lightcurve['arrays']['observations']['ratio_' + band_key][:,0]*100 -100.,
                         yerr= lightcurve['arrays']['observations']['ratio_' + band_key][:,1]*100 ,
                         fmt='.', c='k', alpha=0.25, label='observations')
            plt.errorbar(lightcurve['binned']['observations']['phase'],
                        lightcurve['binned']['observations']['ratio_' + band_key][:, 0]*100 -100.,
                         yerr= lightcurve['binned']['observations']['ratio_' + band_key][:,1]*100 ,
                         fmt='.', c='k', alpha=1.0, label='observations')

            plt.axvspan(-1, lightcurve['bins']['transit_in_bins'][0], alpha=0.25, color='green')
            plt.axvspan(lightcurve['bins']['transit_in_bins'][-1], 1., alpha=0.25, color='green')

            plt.axhline(0, c='C1')
            plt.xlim(lightcurve['arrays']['observations']['phase'][0]-0.01,
                     lightcurve['arrays']['observations']['phase'][-1]+0.01)
            plt.xlabel('orbital phase')
            plt.ylabel('$\mathcal{R}$ - 1. [%]')
            plt.legend()
            plt.show()

    print()
