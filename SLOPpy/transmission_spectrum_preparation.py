from __future__ import print_function, division
from SLOPpy.subroutines.common import *
from SLOPpy.subroutines.spectral_subroutines import *
from SLOPpy.subroutines.io_subroutines import *
from SLOPpy.subroutines.fit_subroutines import *
from SLOPpy.subroutines.shortcuts import *
from SLOPpy.subroutines.plot_subroutines import *

from scipy.interpolate import UnivariateSpline


__all__ = ['compute_transmission_spectrum_preparation',
           'plot_transmission_spectrum_preparation']

def compute_transmission_spectrum_preparation(config_in):

    subroutine_name = 'transmission_spectrum_preparation'

    night_dict = from_config_get_nights(config_in)

    for night in night_dict:

        try:
            preparation = load_from_cpickle('transmission_preparation',
                                          config_in['output'],
                                          night)
            print("{0:45s} Night:{1:15s}   {2:s}".format(subroutine_name, night, 'Retrieved'))
            continue
        except:
            print("{0:45s} Night:{1:15s}   {2:s}".format(subroutine_name, night, 'Computing'))
            print()

        """ Retrieving the list of observations"""
        lists = load_from_cpickle('lists', config_in['output'], night)

        calib_data = load_from_cpickle('calibration_fibA', config_in['output'], night)
        input_data = retrieve_observations(config_in['output'], night, lists['observations'])
        observational_pams = load_from_cpickle('observational_pams', config_in['output'], night)

        if config_in['master-out'].get('use_composite', False):
            master_out = load_from_cpickle('master_out_composite', config_in['output'], night)
            print('  Using composite master-out from all nights')
        else:
            master_out = load_from_cpickle('master_out', config_in['output'], night)

        if config_in['master-out'].get('use_smoothed', False):
            master_out['rescaled'] = master_out['smoothed']
            master_out['rescaled_err'] = master_out['smoothed_err']
            print('  Using smoothed master-out')

        preparation = {
            'subroutine': subroutine_name,
        }

        for obs in lists['observations']:

            preparation[obs] = {}

            preparation[obs]['master_out'] = {}
            preparation[obs]['wave'] = input_data[obs]['wave'] #Added for plotting purpose only

            """ Step 1+2): bring back the master-out to the ORF and rebin the 1D master-out to the 2D observation scale"""
            preparation[obs]['master_out']['rebinned'] = \
                rebin_1d_to_2d(master_out['wave'],
                               master_out['step'],
                               master_out['rescaled'],
                               input_data[obs]['wave'],
                               input_data[obs]['step'],
                               rv_shift=-observational_pams[obs]['rv_shift_ORF2SRF_mod'],
                               preserve_flux=False)

            preparation[obs]['master_out']['rebinned_err'] = \
                rebin_1d_to_2d(master_out['wave'],
                               master_out['step'],
                               master_out['rescaled_err'],
                               input_data[obs]['wave'],
                               input_data[obs]['step'],
                               rv_shift=-observational_pams[obs]['rv_shift_ORF2SRF_mod'],
                               preserve_flux=False,
                               is_error=True)

            for order in range(0, observational_pams['n_orders']):
                preparation[obs]['master_out']['rebinned'][order, :], \
                preparation[obs]['master_out']['rebinned_err'][order, :], \
                _ = \
                    replace_values_errors_with_interpolation_1d(preparation[obs]['master_out']['rebinned'][order, :],
                                                                preparation[obs]['master_out']['rebinned_err'][order, :],
                                                                less_than=0.001, greater_than=5.0000)

            #replace_values_errors(preparation[obs]['master_out']['rebinned'],
            #                      preparation[obs]['master_out']['rebinned_err'],
            #                      threshold=0.0001, replacement=1.0000)

            """ Step 3): obtain the unscaled transmission spectrum for this observation """
            preparation[obs]['ratio'] = input_data[obs]['e2ds']/\
                                      preparation[obs]['master_out']['rebinned']
            preparation[obs]['ratio_err'] = preparation[obs]['ratio'] * \
                                          np.sqrt((input_data[obs]['e2ds_err']/
                                                   input_data[obs]['e2ds'])**2 +
                                                  (preparation[obs]['master_out']['rebinned_err']/
                                                   preparation[obs]['master_out']['rebinned'])**2)

            preparation[obs]['ratio_precleaning'] = preparation[obs]['ratio'].copy()
            preparation[obs]['ratio_precleaning_err'] = preparation[obs]['ratio_err'].copy()


        if night_dict[night].get('spline_residuals', True):

            print()
            print('   Cleaning for telluric residuals with Univariate Spline - threshold about 5%')
            # cleaning using spline_univariate
            for order in range(0, observational_pams['n_orders']):
                obs_reference =  lists['observations'][0]

                len_y = len(lists['observations'])
                len_x = len(preparation[obs_reference]['wave'][order, :])

                time_from_transit = np.empty(len_y, dtype=np.double)
                data_array = np.empty([len_y, len_x], dtype=np.double)
                median_array =  np.empty(len_y, dtype=np.double)
                for i_obs, obs in enumerate(lists['observations']):
                    time_from_transit[i_obs] =  input_data[obs]['BJD'] - observational_pams['time_of_transit']
                    median_array[i_obs] = np.median(preparation[obs]['ratio_precleaning'][order ,:])
                    data_array[i_obs, :] = preparation[obs]['ratio_precleaning'][order ,:]/median_array[i_obs]
                    #wave = preparation[obs]['wave'][order, :]

                res = data_array * 1.
                val = np.empty([len_y, len_x], dtype=np.double)

                for ii in range(0, len_x):
                    spl = UnivariateSpline(time_from_transit, data_array[:, ii])
                    val[:,ii] = spl(time_from_transit)
                    res[:,ii] -= val[:,ii]
                    res[:,ii] /= val[:,ii]

                sel = np.abs(res) > 0.05

                for i_obs, obs in enumerate(lists['observations']):
                    if np.sum(sel[i_obs]) > 0:
                        preparation[obs]['ratio'][order, sel[i_obs]] = val[i_obs, sel[i_obs]] * median_array[i_obs]
                        preparation[obs]['ratio_err'][order, sel[i_obs]] *= 10.
        else:
            print()
            print('   Cleaning for telluric residuals NOT performed')


        for obs in lists['observations']:

            preparation[obs]['deblazed'] = preparation[obs]['ratio'] / calib_data['blaze'] / (input_data[obs]['step'] / np.median(input_data[obs]['step']))
            preparation[obs]['deblazed_err'] = preparation[obs]['ratio_err'] / calib_data['blaze'] / (input_data[obs]['step'] / np.median(input_data[obs]['step']))


            if not config_in['settings'].get('full_output', False):
                del preparation[obs]['master_out']
            else:
                # added for plotting purposes only
                preparation[obs]['rescaling'], \
                preparation[obs]['rescaled'], \
                preparation[obs]['rescaled_err'] = perform_rescaling(
                    preparation[obs]['wave'],
                    preparation[obs]['deblazed'],
                    preparation[obs]['deblazed_err'],
                    observational_pams['wavelength_rescaling'])

        save_to_cpickle('transmission_preparation', preparation, config_in['output'], night)

    print()
    """ Keep going from here after preparation, unless the subroutines has been called just
        to preform the data preparation step
    """


def plot_transmission_spectrum_preparation(config_in, night_input=''):

    subroutine_name = 'transmission_spectrum_preparation'

    night_dict = from_config_get_nights(config_in)

    if night_input == '':
        night_list = night_dict
    else:
        night_list = np.atleast_1d(night_input)




    for night in night_list:

        """ Retrieving the list of observations"""
        lists = load_from_cpickle('lists', config_in['output'], night)
        observational_pams = load_from_cpickle('observational_pams', config_in['output'], night)

        # ! To be removed when testing is done
        # ! This plots do not make any sense anymore 
        input_data = retrieve_observations(config_in['output'], night, lists['observations'])


        """ Retrieving the analysis"""
        try:
            preparation = load_from_cpickle('transmission_preparation', config_in['output'], night)
        except:
            print("No transmission spectrum results, no plots")
            print()
            continue


        #from SLOPpy.subroutines.lines_fit_functions import logprob_case12
        from matplotlib.colors import BoundaryNorm
        from matplotlib.ticker import MaxNLocator

        len_y = len(lists['observations'])
        len_x = 4096
        order= 11

        time_from_transit = np.empty(len_y, dtype=np.double)
        plot_data = np.empty([len_y, len_x], dtype=np.double)

        for i_obs, obs in enumerate(lists['observations']):
            time_from_transit[i_obs] =  input_data[obs]['BJD'] - observational_pams['time_of_transit']
            plot_data[i_obs, :] = preparation[obs]['deblazed'][order ,:]/ np.median(preparation[obs]['deblazed'][order ,:])
            wave = preparation[obs]['wave'][order, :]


        wave_meshgrid, time_meshgrid = np.meshgrid(wave, time_from_transit)

        cmap = plt.get_cmap('coolwarm')

        #levels = MaxNLocator(nbins=15).tick_values(
        #    plot_data.min(), plot_data.max())
        levels = MaxNLocator(nbins=21).tick_values(0.90, 1.10)
        norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

        plt.figure(figsize=(15, 10))
        plt.title('Transmission map in observer reference frame\n {0:s}'.format(night))

        PCF = plt.contourf(wave_meshgrid, time_meshgrid,
                            plot_data, levels=levels, cmap=cmap)
        cbar = plt.colorbar(PCF)
        cbar.ax.set_ylabel('Intensity')
        plt.show()


        if night_dict[night].get('spline_residuals', True):
            res = plot_data * 1.
            from scipy.interpolate import UnivariateSpline
            for ii in range(0,4096):
                spl = UnivariateSpline(time_from_transit, plot_data[:, ii])
                val = spl(time_from_transit)
                res[:,ii] -= val
                res[:,ii] /= val



            cmap = plt.get_cmap('coolwarm')

            levels = MaxNLocator(nbins=10).tick_values(-0.05, 0.05)
            norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
            plt.figure(figsize=(15, 10))
            plt.title('Residuals after dividing by UnivariateSpline Spline\n {0:s}'.format(night))

            PCF = plt.contourf(wave_meshgrid, time_meshgrid,
                                res, levels=levels, cmap=cmap)
            cbar = plt.colorbar(PCF)
            cbar.ax.set_ylabel('Intensity')
            plt.show()





        """ Creation of the color array, based on the BJD of the observations
        """
        colors_properties, colors_plot, colors_scatter = make_color_array_matplotlib3(lists, observational_pams)

        fig = plt.figure(figsize=(12, 6))

        gs = GridSpec(1, 2, width_ratios=[50, 1])
        ax1 = plt.subplot(gs[0, 0])

        ax1.set_ylim(0.90, 1.10)
        #ax2 = plt.subplot(gs[1, 0], sharex=ax1, sharey=ax1)
        cbax1 = plt.subplot(gs[:, 1])

        for obs in lists['transit_in']:

            #preparation[obs]['rescaling'], \
            #preparation[obs]['rescaled'], \
            #preparation[obs]['rescaled_err'] = perform_rescaling(
            #    preparation[obs]['wave'],
            #    preparation[obs]['deblazed'] / (input_data[obs]['step'] / np.median(input_data[obs]['step'])),
            #    preparation[obs]['deblazed_err'] / (input_data[obs]['step'] / np.median(input_data[obs]['step'])),
            #    observational_pams['wavelength_rescaling'])

            preparation[obs]['rescaling'], \
            preparation[obs]['rescaled'], \
            preparation[obs]['rescaled_err'] = perform_rescaling(
                preparation[obs]['wave'],
                preparation[obs]['deblazed'],
                preparation[obs]['deblazed_err'],
                observational_pams['wavelength_rescaling'])


            ax1.scatter(preparation[obs]['wave'],
                    preparation[obs]['rescaled'],
                    s=1, alpha=0.25,
                    color=colors_plot['mBJD'][obs])

        sm = plt.cm.ScalarMappable(cmap=colors_properties['cmap'], norm=colors_properties['norm']['mBJD'])
        sm.set_array([])  # You have to set a dummy-array for this to work...
        cbar = plt.colorbar(sm, cax=cbax1)
        cbar.set_label('BJD - 2450000.0')
        fig.subplots_adjust(wspace=0.05, hspace=0.4)
        plt.show()
