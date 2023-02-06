from __future__ import print_function, division
from SLOPpy.subroutines.common import *
from SLOPpy.subroutines.spectral_subroutines import *
from SLOPpy.subroutines.io_subroutines import *
from SLOPpy.subroutines.fit_subroutines import *
from SLOPpy.subroutines.shortcuts import *
from SLOPpy.subroutines.rebin_subroutines import *
from SLOPpy.subroutines.clv_rm_subroutines import *
from astropy.convolution import convolve, Box1DKernel

__all__ = ['compute_spectra_lightcurve',
           'compute_spectra_lightcurve_clv_rm_correction',
           'plot_spectra_lightcurve',
           'plot_spectra_lightcurve_clv_rm_correction']


def compute_spectra_lightcurve_clv_rm_correction(config_in, lines_label):
    compute_spectra_lightcurve(config_in, lines_label)


def plot_spectra_lightcurve_clv_rm_correction(config_in, night_input=''):
    plot_spectra_lightcurve(config_in, night_input)


def compute_spectra_lightcurve(config_in, lines_label):


    subroutine_name = 'spectra_lightcurve'

    sampler_name = 'emcee'

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

    processed_template = {
        'subroutine': subroutine_name,
        'range': lines_dict['range'],
        'wave': shared_data['coadd']['wave'][shared_selection],
        'step': shared_data['coadd']['step'][shared_selection],
        'size': np.int(np.sum(shared_selection)),
    }


    # doublet sodium in the lab reference frame

    """
    C stands for central
    """
    C_bands = {}
    for passband_key, passband_val in spectral_lines['passbands'].items():
        C_bands[passband_key] = {}
        for line_key, line_val in spectral_lines['lines'].items():
            C_bands[passband_key][line_key] = (np.abs(shared_data['coadd']['wave'] - line_val) < passband_val / 2.)

    """
    S stands for side
    """
    S_bands = {}
    for band_key, band_val in spectral_lines['continuum'].items():
        S_bands[band_key] = (shared_data['coadd']['wave'] >= band_val[0]) & (shared_data['coadd']['wave'] <= band_val[1])

    """
        The transit phase [0-1] is divided in N (=5) bins. Two arrays are computed:
        - transit_in_bins: array with the boundaries of the bins, size=N+1
        - transit_in_step: average size of the bin, size=1
    """
    transit_in_bins = np.linspace(
        -planet_dict['transit_duration'][0]/2./planet_dict['period'][0],
        planet_dict['transit_duration'][0]/2./planet_dict['period'][0],
        6
    )
    transit_in_step = np.average(transit_in_bins[1:]-transit_in_bins[:-1])

    for night in night_dict:

        try:
            lightcurve = load_from_cpickle(subroutine_name, config_in['output'], night)
            print("{0:45s} Night:{1:15s}   {2:s}".format(subroutine_name, night, 'Retrieved'))
            continue
        except:
            print("{0:45s} Night:{1:15s}   {2:s}".format(subroutine_name, night, 'Computing'))
            print()


        """ Retrieving the list of observations"""
        lists = load_from_cpickle('lists', config_in['output'], night)

        calib_data = load_from_cpickle('calibration_fibA', config_in['output'], night)
        input_data = retrieve_observations( config_in['output'], night, lists['observations'])
        observational_pams = load_from_cpickle('observational_pams', config_in['output'], night)

        processed = processed_template.copy()

        lightcurve = {
            'subroutine': subroutine_name,
            'arrays': {
                'observations': {
                    'obs_name': np.zeros(len(lists['observations']), dtype=str),
                    'phase': np.zeros(len(lists['observations'])),
                },
                'transit_in': {},
                'transit_out': {},
            },
            'C_bands': C_bands,
            'S_bands': S_bands,
            'average': {},
            'bins': {
                'transit_in_bins': transit_in_bins,
                'transit_in_step': transit_in_step
            }
        }
        """ Adding the C-bands arrays to the dictionary"""

        for band_key in C_bands:
            lightcurve['arrays']['observations']['ratio_' + band_key] = np.zeros([len(lists['observations']), 2])

        transit_out_flag = np.zeros(len(lists['observations']), dtype=bool)
        transit_in_flag = np.zeros(len(lists['observations']), dtype=bool)

        if clv_rm_correction:
            try:
                clv_rm_models = load_from_cpickle('clv_rm_models', config_in['output'], night, lines_label)
            except (FileNotFoundError, IOError):
                clv_rm_models = load_from_cpickle('clv_rm_models', config_in['output'], night)


        for n_obs, obs in enumerate( lists['observations']):

            processed[obs] = {}
            lightcurve[obs] = {}

            processed[obs]['rescaling'], \
            processed[obs]['rescaled'], \
            processed[obs]['rescaled_err'] = perform_rescaling(
                input_data[obs]['wave'], input_data[obs]['e2ds'], input_data[obs]['e2ds_err'],
                observational_pams['wavelength_rescaling'])

            preserve_flux = input_data[obs].get('absolute_flux', True)

            processed[obs]['uncorrected'] = \
                rebin_2d_to_1d(input_data[obs]['wave'],
                               input_data[obs]['step'],
                               input_data[obs]['e2ds'],
                               calib_data['blaze'],
                               processed['wave'],
                               processed['step'],
                               preserve_flux=preserve_flux,
                               rv_shift=observational_pams[obs]['rv_shift_ORF2SRF'])

            processed[obs]['uncorrected_err'] = \
                rebin_2d_to_1d(input_data[obs]['wave'],
                               input_data[obs]['step'],
                               input_data[obs]['e2ds_err'],
                               calib_data['blaze'],
                               processed['wave'],
                               processed['step'],
                               preserve_flux=preserve_flux,
                               rv_shift=observational_pams[obs]['rv_shift_ORF2SRF'])

            if clv_rm_correction:

                rv_shift = 0.0 # we always stay in SRF
                correction, _ = clv_rm_correction_factor_computation(
                    clv_rm_modelling, shared_data['coadd']['wave'], shared_data['coadd']['step'], rv_shift, obs)

                processed[obs]['clv_rm_correction'] = correction
                processed[obs]['corrected'] /= correction
                processed[obs]['corrected_err'] /= correction





            try:
                phase_internal = (observational_pams[obs]['BJD'] - night_dict[night]['time_of_transit'][0])/planet_dict['period'][0]
            except:
                phase_internal = (observational_pams[obs]['BJD'] - night_dict[night]['time_of_transit'])/planet_dict['period'][0]


            processed[obs]['bands'] = {
                'phase': phase_internal
            }

            s_integrated = 0.000
            s_sigmaq_sum = 0.00
            n_bands = 0.00
            for band_key, band_val in S_bands.items():
                if do_average_instead_of_sum:
                    processed[obs]['bands'][band_key] =  \
                        [np.average(processed[obs]['rebinned'][band_val]),
                         np.sum((processed[obs]['rebinned_err'][band_val])**2)
                         / len(processed[obs]['rebinned_err'][band_val])**2]
                else:
                    processed[obs]['bands'][band_key] =  \
                        [np.sum(processed[obs]['rebinned'][band_val]),
                         np.sum((processed[obs]['rebinned_err'][band_val])**2)]

                s_integrated += processed[obs]['bands'][band_key][0]
                s_sigmaq_sum += processed[obs]['bands'][band_key][1]
                n_bands += 1

            s_integrated *= (2. / n_bands)
            s_sigmaq_sum *= (2. / n_bands)**2

            s_factor_term = np.power(s_integrated, -2.0)

            for band_key, band_dict in C_bands.items():
                processed[obs]['bands'][band_key] = {}

                c_integrated = 0.000
                c_sigmaq_sum = 0.000
                n_bands = 0.00

                for line_key, line_val in band_dict.items():
                    if do_average_instead_of_sum:
                        processed[obs]['bands'][band_key][line_key] = \
                            [np.average(processed[obs]['rebinned'][line_val]),
                             np.sum((processed[obs]['rebinned_err'][line_val]) ** 2)
                             / len(processed[obs]['rebinned_err'][line_val]) ** 2]
                    else:
                        processed[obs]['bands'][band_key][line_key] = \
                            [np.sum(processed[obs]['rebinned'][line_val]),
                             np.sum((processed[obs]['rebinned_err'][line_val]) ** 2)]

                    c_integrated += processed[obs]['bands'][band_key][line_key][0]
                    c_sigmaq_sum += processed[obs]['bands'][band_key][line_key][1]
                    n_bands += 1

                c_integrated *= (2. / n_bands)
                c_sigmaq_sum *= (2. / n_bands) ** 2

                ratio = c_integrated / s_integrated
                lightcurve[obs]['ratio_' + band_key] = [ratio, ratio * np.sqrt( c_sigmaq_sum / c_integrated ** 2 +
                                                                                s_sigmaq_sum / s_integrated ** 2)]
                                                        #np.sqrt(s_factor_term
                                                        #      * (c_sigmaq_sum + ratio**2 * s_sigmaq_sum))]

                lightcurve['arrays']['observations']['ratio_' + band_key][n_obs, :] = \
                    lightcurve[obs]['ratio_' + band_key][:]

            lightcurve[obs]['phase'] = processed[obs]['bands']['phase']
            lightcurve['arrays']['observations']['obs_name'][n_obs] = obs
            lightcurve['arrays']['observations']['phase'][n_obs] = lightcurve[obs]['phase']

            if obs in lists['transit_out']:
                transit_out_flag[n_obs] = True
            else:
                transit_in_flag[n_obs] = True


        for band_key in C_bands:
            lightcurve['arrays']['rescaling_' + band_key] = \
                np.average(lightcurve['arrays']['observations']['ratio_' + band_key][transit_out_flag, 0], axis=0)

        sorting_index = np.argsort(lightcurve['arrays']['observations']['phase'])

        transit_out_flag = transit_out_flag[sorting_index]
        transit_in_flag = transit_in_flag[sorting_index]

        lightcurve['arrays']['observations']['obs_name'] = lightcurve['arrays']['observations']['obs_name'][sorting_index]
        lightcurve['arrays']['observations']['phase'] = lightcurve['arrays']['observations']['phase'][sorting_index]

        lightcurve['arrays']['transit_in']['obs_name'] = lightcurve['arrays']['observations']['obs_name'][transit_in_flag]
        lightcurve['arrays']['transit_in']['phase'] = lightcurve['arrays']['observations']['phase'][transit_in_flag]

        lightcurve['arrays']['transit_out']['obs_name'] = lightcurve['arrays']['observations']['obs_name'][transit_out_flag]
        lightcurve['arrays']['transit_out']['phase'] = lightcurve['arrays']['observations']['phase'][transit_out_flag]

        for band_key in C_bands:
            lightcurve['arrays']['observations']['ratio_' + band_key] = \
                lightcurve['arrays']['observations']['ratio_' + band_key][sorting_index] \
                 / lightcurve['arrays']['rescaling_' + band_key]
            lightcurve['arrays']['transit_in']['ratio_' + band_key] = \
                lightcurve['arrays']['observations']['ratio_' + band_key][transit_in_flag]
            lightcurve['arrays']['transit_out']['ratio_' + band_key] = \
                lightcurve['arrays']['observations']['ratio_' + band_key][transit_out_flag]

            avg_out, avg_out_sq = \
                np.average(lightcurve['arrays']['transit_out']['ratio_' + band_key][:, 0],
                           weights=1./(lightcurve['arrays']['transit_out']['ratio_' + band_key][:, 1])**2,
                           returned=True)
            avg_in, avg_in_sq = \
                np.average(lightcurve['arrays']['transit_in']['ratio_' + band_key][:, 0],
                           weights=1. / (lightcurve['arrays']['transit_in']['ratio_' + band_key][:, 1]) ** 2,
                           returned=True)

            lightcurve['average'][band_key] = {
                'average_out': np.asarray([avg_out, 1./np.power(avg_out_sq, 0.5)]),
                'average_in': np.asarray([avg_in, 1. / np.power(avg_in_sq, 0.5)]),
            }

            delta_fac = \
                    lightcurve['average'][band_key]['average_in'][0]/lightcurve['average'][band_key]['average_out'][0]
            delta_err = delta_fac * np.sqrt(
                (lightcurve['average'][band_key]['average_out'][1]
                 / lightcurve['average'][band_key]['average_out'][0]) ** 2
                + (lightcurve['average'][band_key]['average_in'][1]
                   / lightcurve['average'][band_key]['average_in'][0]) ** 2)

            lightcurve['average'][band_key]['delta'] = np.asarray([(1.-delta_fac)*100., delta_err*100.])

        lightcurve['arrays']['observations']['transit_out_flag'] = transit_out_flag
        lightcurve['arrays']['observations']['transit_in_flag'] = transit_in_flag


        """ Compute the duration of the pre-transit observations, using as scale
            the number of bins, with the same size as those used inside the transit.
            
            The value is given by the difference of the phase of the beginning of the transit minus
            the phase of the first observation, keeping in mind that the centre of the transit has phase = 0
            
            An additional bin is added if there are observations left out from the actual number of bins
        """
        pre_duration = transit_in_bins[0] - lightcurve['arrays']['transit_out']['phase'][0]
        if pre_duration > 0:
            nsteps_pre = int(pre_duration/transit_in_step)
            if pre_duration % transit_in_step > 0.0:
                nsteps_pre += 1
        else:
            nsteps_pre = 0

        """ same as pre-transit, but suing the post-transit instead"""
        post_duration = lightcurve['arrays']['transit_out']['phase'][-1] - transit_in_bins[-1]
        if post_duration > 0:
            nsteps_post = int(post_duration / transit_in_step)
            if post_duration % transit_in_step > 0.0:
                nsteps_post += 1
        else:
            nsteps_post = 0

        """ THe full array with both in-transit and out-transit phase, built in such a way that the 
            - the lower boundary of the first in-transit bin corresponds to the beginning of the transit
            - the upper boundary of the last in-transit bin corresponds to the end of the transit
        """
        transit_bins = np.arange(transit_in_bins[0]-nsteps_pre*transit_in_step,
                                 transit_in_bins[-1] + (nsteps_post+1.1) * transit_in_step,
                                 transit_in_step)

        lightcurve['binned'] = {
            'observations': {
                'phase': np.zeros(len(transit_bins)),
            },
            'transit_in': {},
            'transit_out': {},
        }
        for band_key in C_bands:
            lightcurve['binned']['observations']['ratio_' + band_key] = np.zeros([len(transit_bins), 2])

        transit_out_flag = np.zeros(len(transit_bins), dtype=bool)
        transit_in_flag = np.zeros(len(transit_bins), dtype=bool)

        n_a = 0
        for nb in range(0, len(transit_bins)-1):
            sel = (lightcurve['arrays']['observations']['phase'] >= transit_bins[nb]) \
                   & (lightcurve['arrays']['observations']['phase'] < transit_bins[nb+1])

            if np.sum(sel) <= 0: continue
            lightcurve['binned']['observations']['phase'][n_a] = np.average(lightcurve['arrays']['observations']['phase'][sel])

            for band_key in C_bands:

                lightcurve['binned']['observations']['ratio_' + band_key][n_a, 0], sum_weights = np.average(
                    lightcurve['arrays']['observations']['ratio_' + band_key][sel, 0],
                    weights=1. / lightcurve['arrays']['observations']['ratio_' + band_key][sel, 1]**2,
                    returned=True)

                lightcurve['binned']['observations']['ratio_' + band_key][n_a, 1] = np.sqrt(1. / sum_weights)

            if np.abs(lightcurve['binned']['observations']['phase'][n_a]) >= \
                planet_dict['transit_duration'][0]/2./planet_dict['period'][0]:
                transit_out_flag[n_a] = True
            else:
                transit_in_flag[n_a] = True

            n_a += 1 # bins actually computed

        lightcurve['binned']['transit_in']['phase'] = lightcurve['binned']['observations']['phase'][transit_in_flag]
        lightcurve['binned']['transit_out']['phase'] = lightcurve['binned']['observations']['phase'][transit_out_flag]
        lightcurve['binned']['observations']['phase'] = lightcurve['binned']['observations']['phase'][:n_a]

        for band_key in C_bands:
            lightcurve['binned']['transit_in']['ratio_' + band_key] = \
                lightcurve['binned']['observations']['ratio_' + band_key][transit_in_flag, :]
            lightcurve['binned']['transit_out']['ratio_' + band_key] = \
                lightcurve['binned']['observations']['ratio_' + band_key][transit_out_flag, :]
            lightcurve['binned']['observations']['ratio_' + band_key] = \
                lightcurve['binned']['observations']['ratio_' + band_key][:n_a, :]

        save_to_cpickle(subroutine_name+'_processed', processed, config_in['output'], night)
        save_to_cpickle(subroutine_name, lightcurve, config_in['output'], night)


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
