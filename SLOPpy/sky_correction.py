from __future__ import print_function, division
from SLOPpy.subroutines.common import *
from SLOPpy.subroutines.spectral_subroutines import *
from SLOPpy.subroutines.fit_subroutines import *
from SLOPpy.subroutines.io_subroutines import *
from SLOPpy.subroutines.plot_subroutines import *

__all__ = ["compute_sky_correction", "plot_sky_correction"]

subroutine_name = 'sky_correction'

def compute_sky_correction(config_in):

    night_dict = from_config_get_nights(config_in)

    for night in night_dict:

        try:
            processed = load_from_cpickle('skycorrected_fibA', config_in['output'], night)
            print("{0:45s} Night:{1:15s}   {2:s}".format(subroutine_name, night, 'Retrieved'))
            continue
        except:

            """ Retrieving the list of observations"""
            lists = load_from_cpickle('lists', config_in['output'], night)

            """ Retrieving the observations and calibration data for fiber B, if they exist"""
            try:
                input_data_B = load_from_cpickle('input_dataset_fibB', config_in['output'], night)
                calib_data_B = load_from_cpickle('calibration_fibB', config_in['output'], night)
                print("{0:45s} Night:{1:15s}   {2:s}".format(subroutine_name, night, 'Computing'))
            except:
                print("{0:45s} Night:{1:15s}   {2:s}".format(subroutine_name, night, 'Skipped'))
                continue

        """ Retrieving the observations and calibration data for fiber A"""
        print()
        print("  Retrieving the data for night ", night)
        input_data_A = load_from_cpickle('input_dataset_fibA', config_in['output'], night)
        calib_data_A = load_from_cpickle('calibration_fibA', config_in['output'], night)

        map_orders_A = calib_data_A['fibAB_orders_match']
        map_orders_B = calib_data_B['fibAB_orders_match']
        """map_orders_A = [0,1,2,3] = map_orders_B"""

        processed = {
            'subroutine': 'sky_correction'
        }

        for obs in lists['observations']:
            processed[obs] = {}

        " computing the ratio between the lamp flux of fiber A and  B"
        print()
        print("  Computing the ratio between the lamp flux of fiber A and  B")

        processed['ratioAB'] = calib_data_A['lamp'][map_orders_A, :]/calib_data_B['lamp'][map_orders_B, :]

        first_obs = lists['observations'][0]

        wave_difference = \
            input_data_A[first_obs]['wave'][map_orders_A, :] - input_data_B[first_obs]['wave'][map_orders_B, :]

        print()
        print("   Wavelength difference between fiber A and B: ", \
            np.average(wave_difference), " +- ", np.std(wave_difference), " \AA")

        if np.abs(np.average(wave_difference)) > 0.006 or np.std(wave_difference) > 0.006:
            raise ValueError("TO BE IMPLEMENTED!!!!!!!   ")
            quit()

        else:
            """ We assume that the relative RV shift between finber A and fiber B in the pixel scale is
            is minimal """

            for obs in lists['observations']:
                processed[obs]['sky_fibA'] = np.zeros([input_data_A['n_orders'], input_data_A['n_pixels']])
                processed[obs]['sky_fibA'][map_orders_A, :] = \
                    processed['ratioAB'] * input_data_B[obs]['e2ds'][map_orders_B, :]
                processed[obs]['e2ds'] = input_data_A[obs]['e2ds'] - processed[obs]['sky_fibA']
                processed[obs]['e2ds_err'] = np.sqrt(
                    input_data_A[obs]['e2ds_err'][map_orders_A, :] ** 2 +
                    (processed['ratioAB'] * input_data_B[obs]['e2ds_err'][map_orders_B, :]) ** 2)

                """ Zero or negative values are identified, flagged and substituted with another value """
                #replacement = 0.1
                #processed[obs]['null'] = (processed[obs]['e2ds'] <= replacement)
                #processed[obs]['e2ds'][processed[obs]['null']] = replacement

        save_to_cpickle('skycorrected_fibA', processed, config_in['output'], night)


def plot_sky_correction(config_in, night_input=''):

    night_dict = from_config_get_nights(config_in)

    if night_input=='':
        night_list = night_dict
    else:
        night_list = np.atleast_1d(night_input)

    for night in night_list:

        print("plot_sky_correction                        Night: ", night)

        """ Retrieving the list of observations"""
        lists = load_from_cpickle('lists', config_in['output'], night)

        """ Retrieving the observations and calibration data for fiber B, if they exist"""
        try:
            input_data_B = load_from_cpickle('input_dataset_fibB', config_in['output'], night)
            calib_data_B = load_from_cpickle('calibration_fibB', config_in['output'], night)
        except:
            print("No fiber_B dataset available, skipping sky correction plot")
            continue

        observational_pams = load_from_cpickle('observational_pams', config_in['output'], night)

        input_data_A = load_from_cpickle('input_dataset_fibA', config_in['output'], night)
        processed = load_from_cpickle('skycorrected_fibA', config_in['output'], night)

        colors_properties, colors_plot, colors_scatter = make_color_array_matplotlib3(lists, observational_pams)
        fig, gs, cbax1, ax1, ax2, ax3 = grid_3plot_small()

        for i, obs in enumerate(lists['observations']):
            for k in range(0, input_data_B[obs]['n_orders']):
                if i == 0 and k == 0:
                    ax2.scatter(input_data_B[obs]['wave'][k, :], input_data_B[obs]['e2ds'][k, :],
                                c=colors_scatter['mBJD'][obs], s=2, alpha=0.5, label='Sky observations (ORF)')
                else:
                    ax2.scatter(input_data_B[obs]['wave'][k, :], input_data_B[obs]['e2ds'][k, :],
                                c=colors_scatter['mBJD'][obs], s=2, alpha=0.5)

            for k in range(0, input_data_A[obs]['n_orders']):

                if i == 0 and k == 0:
                    ax1.scatter(input_data_A[obs]['wave'][k, :], input_data_A[obs]['e2ds'][k, :],
                                c=colors_scatter['mBJD'][obs], s=1, alpha=0.5, label='Target observations (ORF)')

                else:
                    ax1.scatter(input_data_A[obs]['wave'][k, :], input_data_A[obs]['e2ds'][k, :],
                                c=colors_scatter['mBJD'][obs], s=1, alpha=0.5)

                ax3.scatter(input_data_A[obs]['wave'][k, :], processed[obs]['e2ds'][k, :],
                            c=colors_scatter['mBJD'][obs], s=1, alpha=0.5)

        ax1.set_title('Night: {0:s} \n Input spectra'.format(night))
        ax1.legend(loc=1)

        ax2.set_title('Sky spectrum from fiber B')

        ax3.set_title('After Sky correction')
        ax3.set_xlabel('$\lambda$ [$\AA$]')

        sm = plt.cm.ScalarMappable(cmap=colors_properties['cmap'], norm=colors_properties['norm']['mBJD'])
        sm.set_array([])  # You have to set a dummy-array for this to work...
        cbar = plt.colorbar(sm, cax=cbax1)
        cbar.set_label('BJD - 2450000.0')
        plt.show()
