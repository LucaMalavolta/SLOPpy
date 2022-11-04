from __future__ import print_function, division
from SLOPpy.subroutines.common import *
from SLOPpy.subroutines.spectral_subroutines import *
from SLOPpy.subroutines.fit_subroutines import *
from SLOPpy.subroutines.io_subroutines import *
from SLOPpy.subroutines.shortcuts import *
from SLOPpy.subroutines.plot_subroutines import *

__all__ = ["compute_pca_preparation"]

def compute_pca_preparation(config_in, append_name=None):

    if append_name:
        subroutine_name = 'pca_preparation_' + append_name
        filename = 'pca_preparation_' + append_name
    else:
        subroutine_name = 'pca_preparation'
        filename = 'pca_preparation'

    night_dict = from_config_get_nights(config_in)

    preparation_dict = {
        'fit_iters': 5,
        'fit_order': 5,
        'fit_sigma': 4
    }


    for night in night_dict:

        try:
            preparation = load_from_cpickle(filename, config_in['output'], night)
            print("{0:45s} Night:{1:15s}   {2:s}".format(subroutine_name, night, 'Retrieved'))
            continue
        except:
            print("{0:45s} Night:{1:15s}   {2:s}".format(subroutine_name, night, 'Computing'))
            print()

        """ Retrieving the list of observations"""
        lists = load_from_cpickle('lists', config_in['output'], night)

        """ Retrieving input and calibration data """
        calib_data = load_from_cpickle('calibration_fibA', config_in['output'], night)
        input_data = retrieve_observations(config_in['output'], night, lists['observations'],
                                           use_refraction=False, use_telluric=False, use_interstellar=False,
                          use_telluric_spline= False)
        observational_pams = load_from_cpickle('observational_pams', config_in['output'], night)


        obs_ref =lists['observations'][0]
        n_obs = len(lists['observations'])

        n_orders = input_data[obs_ref]['n_orders']
        n_pixels = input_data[obs_ref]['n_pixels']

        stack_wave = np.zeros([n_obs, n_orders, n_pixels], dtype=np.double)
        stack_e2ds = np.zeros([n_obs, n_orders, n_pixels], dtype=np.double)
        stack_e2ds_err = np.zeros([n_obs, n_orders, n_pixels], dtype=np.double)
        stack_bjd = np.zeros(n_obs, dtype=np.double)
        stack_airmass = np.zeros(n_obs, dtype=np.double)

        import matplotlib.pyplot as pl


        for i_obs, obs in enumerate(lists['observations']):

            stack_wave[i_obs, :, :] = input_data[obs]['wave']
            stack_e2ds[i_obs, :, :] = input_data[obs]['e2ds']
            stack_e2ds_err[i_obs, :, :] = input_data[obs]['e2ds_err']

            stack_bjd[i_obs] = input_data[obs]['BJD']
            stack_airmass[i_obs] = input_data[obs]['AIRMASS']

            median = np.nanmedian(stack_e2ds[i_obs, :, :], axis=0)
            stack_e2ds[i_obs, :, :] /= median
            stack_e2ds_err[i_obs, :, :] /= median

            print(i_obs, obs)
            plt.scatter(stack_wave[i_obs, 22, :], stack_e2ds[i_obs, 22, :], s=2)
            plt.scatter(stack_wave[i_obs, 22, :], stack_e2ds_err[i_obs, 22, :], s=2)

            plt.show()

        poly_flag = (stack_e2ds > 0.01)
        poly_flag[:, :, :50]= False
        poly_flag[:, :, :50]= False

        stack_polyfit = np.zeros_like(stack_e2ds)

        for i_orders in range(0,n_orders):

            order_wave = stack_wave[:, i_orders, :]
            order_e2ds = stack_e2ds[:, i_orders, :]
            order_flag = poly_flag[:, i_orders, :]

            for n_iter in range(0, preparation_dict['fit_iters']):

                coeff_order = np.polynomial.chebyshev.chebfit(
                    order_wave[order_flag],
                    order_e2ds[order_flag],
                    preparation_dict['fit_order'])

                fit_order = \
                    np.polynomial.chebyshev.chebval(order_wave, coeff_order)
                fit_shaped = np.reshape(fit_order, np.shape(order_wave))
                residuals = order_e2ds - fit_shaped

                if n_iter < preparation_dict['fit_iters'] - 1:
                    std = np.std(residuals[order_flag])
                    order_flag = (order_flag) & (residuals > -preparation_dict['fit_sigma'] * std)

            if i_orders==22:
                plt.scatter(stack_wave[:, i_orders,:], stack_e2ds[:, i_orders, :]-0.5, s=2)
                plt.scatter(stack_wave[:, i_orders,:], stack_e2ds[:, i_orders, :]/fit_shaped, s=2)

            stack_e2ds[:, i_orders, :]/=fit_shaped
            stack_e2ds_err[:, i_orders, :]/=fit_shaped
            stack_polyfit[:, i_orders, :] =fit_shaped

        preparation = {
            'stack_e2ds': stack_e2ds,
            'stack_e2ds_err': stack_e2ds_err,
            'stack_bjd': stack_e2ds,
            'stack_airmass': stack_e2ds,
            'stack_polyfit': stack_polyfit,
            'frame': {
                     'n_obs': n_obs,
                     'n_orders': n_orders,
                     'n_pixels': n_pixels,
                },
            'fit_pams': preparation_dict
        }

        plt.show()

        quit()