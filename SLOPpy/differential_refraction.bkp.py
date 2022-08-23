from __future__ import print_function, division
from SLOPpy.subroutines.common import *
from SLOPpy.subroutines.spectral_subroutines import *
from SLOPpy.subroutines.fit_subroutines import *
from SLOPpy.subroutines.io_subroutines import *
from SLOPpy.subroutines.shortcuts import *

__all__ = ["compute_differential_refraction", "plot_differential_refraction"]

subroutine_name = 'differential_refraction'


def compute_differential_refraction(config_in):
    night_dict = from_config_get_nights(config_in)

    print()

    for night in night_dict:

        try:
            refraction = load_from_cpickle('refraction', config_in['output'], night)
            print("{0:45s} Night:{1:15s}   {2:s}".format(subroutine_name, night, 'Retrieved'))
            continue
        except:
            print("{0:45s} Night:{1:15s}   {2:s}".format(subroutine_name, night, 'Computing'))
            print()

        """ Retrieving the list of observations"""
        lists = load_from_cpickle('lists', config_in['output'], night)

        calib_data = load_from_cpickle('calibration_fibA', config_in['output'], night)
        input_data = retrieve_observations(config_in['output'], night, lists['observations'],
                                           use_refraction=False, use_telluric=False)
        observational_pams = load_from_cpickle('observational_pams', config_in['output'], night)

        try:
            processed = load_from_cpickle('refraction_processed_halfway', config_in['output'], night)
            refraction = load_from_cpickle('refraction_halfway', config_in['output'], night)
            print("  Starting from intermediate step ")
        except:

            processed = {
                'subroutine': subroutine_name,
                'coadd': {
                    'wave': input_data['coadd']['wave'],
                    'size': input_data['coadd']['size'],
                    'step': input_data['coadd']['step'],
                    #'flux': np.zeros(input_data['coadd']['size'], dtype=np.double),
                    #'flux_err': np.zeros(input_data['coadd']['size'], dtype=np.double)
                }
            }

            refraction = {
                'subroutine': 'differential_refraction',
                'wave': processed['coadd']['wave']
            }

            total_flux = np.empty([len(lists['observations']), input_data['coadd']['size']], dtype=np.double)
            total_wght = np.zeros([len(lists['observations']), input_data['coadd']['size']], dtype=np.double)
            total_mask = np.ones([len(lists['observations']), input_data['coadd']['size']], dtype=bool)

            print("  Chebyshev polynomial order for differential refraction fit: ",
                  observational_pams['refraction_poly_order'])

            print("  Number of iterations: ",
                  observational_pams['refraction_poly_iters'])
            print()

            """ Rebinning of all the spectra """
            for n_obs, obs in enumerate(lists['observations']):

                print("  Spectral rebinning - Processing: ", obs)
                processed[obs] = {}

                """ Rebinning of the spectra in the SRF, except for a fixed constant in order to minimize
                    the difference between  """

                processed[obs]['flux_rebinned_stellarRF'] = \
                    rebin_2d_to_1d(input_data[obs]['wave'],
                                   input_data[obs]['step'],
                                   input_data[obs]['e2ds'],
                                   calib_data['blaze'],
                                   processed['coadd']['wave'],
                                   processed['coadd']['step'],
                                   rv_shift=observational_pams[obs]['rv_shift_ORF2SRF_mod'])

                processed[obs]['err_flux_rebinned_SRF'] = \
                    rebin_2d_to_1d(input_data[obs]['wave'],
                                   input_data[obs]['step'],
                                   input_data[obs]['e2ds_err'],
                                   calib_data['blaze'],
                                   processed['coadd']['wave'],
                                   processed['coadd']['step'],
                                   rv_shift=observational_pams[obs]['rv_shift_ORF2SRF_mod'],
                                   is_error=True)

                """ Zero or negative values are identified, flagged and substituted with another value """
                processed[obs]['flux_rebinned_stellarRF'], \
                processed[obs]['err_flux_rebinned_SRF'], \
                processed[obs]['flux_rebinned_SRF_null'] = \
                    replace_values_errors_with_interpolation_1d(processed[obs]['flux_rebinned_stellarRF'],
                                                                processed[obs]['err_flux_rebinned_SRF'],
                                                                force_positive=True)

                processed[obs]['rescaling'], processed[obs]['rescaled'], processed[obs]['rescaled_err'] = \
                    perform_rescaling(processed['coadd']['wave'],
                                      processed[obs]['flux_rebinned_stellarRF'],
                                      processed[obs]['err_flux_rebinned_SRF'],
                                      observational_pams['wavelength_rescaling'])

                processed[obs]['rescaled_blazed'] = input_data[obs]['e2ds'] \
                                                    / processed[obs]['rescaling'] \
                                                    / calib_data['blaze']

                if obs in lists['telluric']:
                    total_flux[n_obs, :] = processed[obs]['rescaled']
                    total_mask[n_obs, :] = processed[obs]['flux_rebinned_SRF_null']
                    total_wght[n_obs, :] = 1. / (processed[obs]['rescaled_err'] ** 2)

                    # processed['coadd']['flux'] += processed[obs]['flux_rebinned_stellarRF']
                    # """ SNR (assumed to be the square root of the flux) is added in quadrature """
                    # processed['coadd']['flux_err'] += processed[obs]['err_flux_rebinned_SRF'] ** 2
                    print("      Observation added to reference spectrum")

            #masked_array = np.ma.array(total_flux, mask=total_mask)
            #processed['coadd']['rescaled'], sum_weights = np.ma.average(masked_array,
            #                                                            weights=total_wght,
            #                                                            axis=0,
            #                                                            returned=True)
            # processed['coadd']['rescaled'][sum_weights <= 0.0001] = 1.000
            # sum_weights[sum_weights <= 0.0001] = 0.0001
            # processed['coadd']['rescaled_err'] = 1. / np.sqrt(sum_weights)

            masked_array = np.ma.array(total_flux, mask=total_mask)
            rescaled_mask, sum_weights = np.ma.average(masked_array,
                                                       weights=total_wght,
                                                       axis=0,
                                                       returned=True)

            processed['coadd']['rescaled'] = rescaled_mask.filled(0.00)
            sum_weights[sum_weights <= 0.0] = 1.0
            processed['coadd']['rescaled_err'] = 1. / np.sqrt(sum_weights)

            processed['coadd']['rescaled'], processed['coadd']['rescaled_err'], processed['coadd']['null'] = \
                replace_values_errors_with_interpolation_1d(processed['coadd']['rescaled'],
                                                            processed['coadd']['rescaled_err'],
                                                            force_positive=True)

            save_to_cpickle('refraction_processed_halfway', processed, config_in['output'], night)
            save_to_cpickle('refraction_halfway', refraction, config_in['output'], night)

        """ Now each observation is divided by the reference spectrum, after being redshifted in the observer RF
            The result is then used to model the flux variation
        """

        for obs in lists['observations']:

            print("  Division by reference spectrum and fit of the flux variation: ", obs)

            """ Going back to the observer RF and rebinning the spectrum into the observed orders """
            processed[obs]['master_flux'] = \
                rebin_1d_to_2d(processed['coadd']['wave'],
                               processed['coadd']['step'],
                               processed['coadd']['rescaled'],
                               input_data[obs]['wave'],
                               input_data[obs]['step'],
                               rv_shift=-observational_pams[obs]['rv_shift_ORF2SRF_mod'])

            processed[obs]['master_ferr'] = \
                rebin_1d_to_2d(processed['coadd']['wave'],
                               processed['coadd']['step'],
                               processed['coadd']['rescaled_err'],
                               input_data[obs]['wave'],
                               input_data[obs]['step'],
                               rv_shift=-observational_pams[obs]['rv_shift_ORF2SRF_mod'],
                               is_error=True)

            """ Zero or negative values are identified, flagged and substituted with another value """
            processed[obs]['master_flux'], processed[obs]['master_ferr'], processed[obs]['master_null'] = \
                replace_values_errors_with_interpolation_2d(processed[obs]['master_flux'],
                                                            processed[obs]['master_ferr'],
                                                            less_than=0.001)

            processed[obs]['ratio'] = processed[obs]['rescaled_blazed'] / processed[obs]['master_flux']
            """
            processed[obs]['ratio'] = input_data[obs]['e2ds']\
                                      /processed[obs]['rescaling']\
                                      / (processed[obs]['master_flux'] * calib_data['blaze'])
            """

            refraction[obs] = {}
            refraction[obs]['polyfit_e2ds'] = np.zeros([input_data[obs]['n_orders'], input_data[obs]['n_pixels']])
            processed[obs]['residuals'] = np.zeros([input_data[obs]['n_orders'], input_data[obs]['n_pixels']])
            refraction[obs]['poly_flag'] = np.zeros([input_data[obs]['n_orders'], input_data[obs]['n_pixels']],
                                                    dtype=bool)

            for order in range(0, input_data[obs]['n_orders']):

                order_coeff_name = 'order_' + repr(order)
                refraction[obs]['poly_flag'][order, :] = (processed[obs]['ratio'][order, :] > 0.1)
                refraction[obs]['poly_flag'][order, :50] = False
                refraction[obs]['poly_flag'][order, -50:] = False

                for n_iter in range(0, observational_pams['refraction_poly_iters']):

                    refraction[obs][order_coeff_name] = np.polynomial.chebyshev.chebfit(
                        input_data[obs]['wave'][order, refraction[obs]['poly_flag'][order, :]],
                        processed[obs]['ratio'][order, refraction[obs]['poly_flag'][order, :]],
                        observational_pams['refraction_poly_order'])

                    refraction[obs]['polyfit_e2ds'][order, :] = \
                        np.polynomial.chebyshev.chebval(input_data[obs]['wave'][order, :],
                                                        refraction[obs][order_coeff_name])

                    processed[obs]['residuals'][order, :] = refraction[obs]['polyfit_e2ds'][order, :]\
                                                            - processed[obs]['ratio'][order, :]

                    if n_iter < observational_pams['refraction_poly_iters'] - 1:
                        std = np.std(processed[obs]['residuals'][order, :])
                        refraction[obs]['poly_flag'][order, :] = (refraction[obs]['poly_flag'][order, :]) \
                                                                 & (np.abs(processed[obs]['residuals'][order, :]) <
                                                                    observational_pams['refraction_poly_sigma'] * std)

            processed[obs]['e2ds_corrected'] = input_data[obs]['e2ds'] / refraction[obs]['polyfit_e2ds']
            processed[obs]['e2ds_corrected_err'] = input_data[obs]['e2ds_err'] / refraction[obs]['polyfit_e2ds']

            processed[obs]['flux_rebinned_stellarRF_corrected'] = \
                rebin_2d_to_1d(input_data[obs]['wave'],
                               input_data[obs]['step'],
                               processed[obs]['e2ds_corrected'],
                               calib_data['blaze'],
                               processed['coadd']['wave'],
                               processed['coadd']['step'],
                               rv_shift=observational_pams[obs]['rv_shift_ORF2SRF_mod'])

            processed[obs]['err_flux_rebinned_SRF_corrected'] = \
                rebin_2d_to_1d(input_data[obs]['wave'],
                               input_data[obs]['step'],
                               processed[obs]['e2ds_corrected_err'],
                               calib_data['blaze'],
                               processed['coadd']['wave'],
                               processed['coadd']['step'],
                               rv_shift=observational_pams[obs]['rv_shift_ORF2SRF_mod'],
                               is_error=True)

            processed[obs]['flux_rebinned_stellarRF_corrected'], \
            processed[obs]['err_flux_rebinned_SRF_corrected'], _ = \
                replace_values_errors_with_interpolation_1d(processed[obs]['flux_rebinned_stellarRF_corrected'],
                                                            processed[obs]['err_flux_rebinned_SRF_corrected'],
                                                            less_than=0.001)

        save_to_cpickle('refraction_processed', processed, config_in['output'], night)
        save_to_cpickle('refraction', refraction, config_in['output'], night)


def plot_differential_refraction(config_in, night_input=''):
    night_dict = from_config_get_nights(config_in)

    if night_input == '':
        night_list = night_dict
    else:
        night_list = np.atleast_1d(night_input)

    for night in night_list:

        """ Retrieving the list of observations"""
        lists = load_from_cpickle('lists', config_in['output'], night)

        """ Retrieving the observations"""
        input_data = retrieve_observations(config_in['output'], night, lists['observations'],
                                           use_refraction=False, use_telluric=False)
        #observational_pams = load_from_cpickle('observational_pams', config_in['output'], night)


        try:
            """ Retrieving the analysis"""
            processed = load_from_cpickle('refraction_processed', config_in['output'], night)
            refraction = load_from_cpickle('refraction', config_in['output'], night)
        except:
            print("  Failed in retrieving processed data")
            return

        """ Creation of the color array, based on the BJD of the observations
        """
        bjd = []
        am = []

        for obs in lists['observations']:
            bjd.append(input_data[obs]['BJD'] - 2450000.0)
            am.append(input_data[obs]['AIRMASS'])

        color_cmap = plt.cm.viridis
        color_norm = plt.Normalize(vmin=bjd[0], vmax=bjd[-1])
        colors = color_cmap(color_norm(np.asarray(bjd)))

        offset = 0.10

        y_limits = [0.8, 1.2]


        """
        fig = plt.figure(figsize=(12, 6))
        gs = GridSpec(2, 2, width_ratios=[50, 1])

        ax1 = plt.subplot(gs[0, 0])
        ax2 = plt.subplot(gs[1, 0], sharex=ax1, sharey=ax1)
        cbax1 = plt.subplot(gs[:, 1])

        for i, obs in enumerate(lists['observations']):
            shift = i/10.0

            for order in range(0, input_data[obs]['n_orders']):
                ax1.scatter(input_data[obs]['wave'][order, :],
                         processed[obs]['rescaled_blazed'][order, :] - shift,
                         c=line_colors[i], s=1, alpha=0.5)

                ax1.plot(input_data[obs]['wave'][order, :],
                         processed[obs]['master_flux'][order, :], - shift,
                         c='k', lw=1)

                ax2.scatter(input_data[obs]['wave'][order, :],
                         processed[obs]['rescaled_blazed'][order, :]/refraction[obs]['polyfit_e2ds'][order, :] - shift,
                         c=line_colors[i], s=1, alpha=0.5)

                ax2.plot(input_data[obs]['wave'][order, :],
                         processed[obs]['master_flux'][order, :], - shift,
                         c='k', lw=1)


        ax1.set_xlim(processed['coadd']['wave'][0], processed['coadd']['wave'][-1])
        ax1.legend(loc=3)
        ax1.set_title('Night: ' + night)

        ax2.set_xlabel('$\lambda$ [$\AA$]')

        sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=colors[0], vmax=colors[-1]))
        sm.set_array([])  # You have to set a dummy-array for this to work...
        cbar = plt.colorbar(sm, cax=cbax1)
        cbar.set_label('BJD - 2450000.0')
        fig.subplots_adjust(wspace=0.05, hspace=0.4)
        plt.show()
        """

        fig = plt.figure(figsize=(12, 6))
        gs = GridSpec(2, 2, width_ratios=[50, 1])

        ax1 = plt.subplot(gs[0, 0])
        ax2 = plt.subplot(gs[1, 0], sharex=ax1, sharey=ax1)
        cbax1 = plt.subplot(gs[:, 1])

        for i, obs in enumerate(lists['observations']):
            color = [colors[i][:-1]]

            if i==0:
                ax1.scatter(processed['coadd']['wave'],
                         processed[obs]['flux_rebinned_stellarRF'] / processed[obs]['rescaling'],
                         c=color, s=2, alpha=0.2, label='observation')
            else:
                ax1.scatter(processed['coadd']['wave'],
                         processed[obs]['flux_rebinned_stellarRF'] / processed[obs]['rescaling'],
                         c=color, s=2, alpha=0.2)

            ax2.scatter(processed['coadd']['wave'],
                     processed[obs]['flux_rebinned_stellarRF_corrected'] / processed[obs]['rescaling'],
                     c=color, s=3, alpha=0.2)

        ax1.plot(processed['coadd']['wave'], processed['coadd']['rescaled'], c='k', lw=1, label='reference spectrum')
        ax2.plot(processed['coadd']['wave'], processed['coadd']['rescaled'], c='k', lw=1)


        ax1.set_xlim(processed['coadd']['wave'][0], processed['coadd']['wave'][-1])
        ax1.set_ylim(y_limits)
        ax2.set_ylim(y_limits)
        ax1.legend(loc=1)
        ax1.set_title('Night: {0:s} \n Input spectra'.format(night))
        ax2.set_title('Corrected spectra')
        ax2.set_xlabel('$\lambda$ [$\AA$]')

        sm = plt.cm.ScalarMappable(cmap=color_cmap, norm=color_norm)
        sm.set_array([])  # You have to set a dummy-array for this to work...
        cbar = plt.colorbar(sm, cax=cbax1)
        cbar.set_label('BJD - 2450000.0')
        fig.subplots_adjust(wspace=0.05, hspace=0.4)
        plt.show()


        """
        PLOT
        """

        fig = plt.figure(figsize=(12, 6))
        gs = GridSpec(1, 2, width_ratios=[50, 1])
        ax = plt.subplot(gs[0, 0])
        cbax1 = plt.subplot(gs[:, 1])

        for i, obs in enumerate(lists['observations']):
            if i == 0:
                offset = np.std(processed[obs]['ratio'][refraction[obs]['poly_flag']].flatten()) * 6
                average = np.average(processed[obs]['ratio'][refraction[obs]['poly_flag']].flatten())
                y_limits = [average-offset, average+offset]

            color = [colors[i][:-1]]

            for order in range(0, input_data[obs]['n_orders']):

                ax.scatter(input_data[obs]['wave'][refraction[obs]['poly_flag']],
                           processed[obs]['ratio'][refraction[obs]['poly_flag']] + offset*i,
                           s=1, c=color, alpha=0.50, zorder=2)

                ax.scatter(input_data[obs]['wave'][~refraction[obs]['poly_flag']],
                           processed[obs]['ratio'][~refraction[obs]['poly_flag']] + offset*i,
                           s=2, c='k', alpha=0.05, zorder=1)
                ax.plot(input_data[obs]['wave'][order, :],
                        refraction[obs]['polyfit_e2ds'][order, :] + offset*i,
                        c='k', lw=1, zorder=5)

        y_limits_offset = [min(y_limits[0] + offset * i, y_limits[0]),
                           max(y_limits[1] + offset * i, y_limits[1])]

        ax.set_ylim(y_limits_offset)
        ax.set_xlabel('$\lambda$ [$\AA$]')
        ax.legend(loc=3)
        ax.set_title('Night: {0:s} \n Fit of the ratio obs/master'.format(night))

        sm = plt.cm.ScalarMappable(cmap=color_cmap, norm=color_norm)
        sm.set_array([])  # You have to set a dummy-array for this to work...
        cbar = plt.colorbar(sm, cax=cbax1)
        cbar.set_label('BJD - 2450000.0')
        fig.subplots_adjust(wspace=0.05, hspace=0.4)
        plt.show()


        """ 
        PLOT: residuals of the fit
        """

        fig = plt.figure(figsize=(12, 6))
        gs = GridSpec(1, 2, width_ratios=[50, 1])
        ax = plt.subplot(gs[0, 0])
        cbax1 = plt.subplot(gs[:, 1])

        for i, obs in enumerate(lists['observations']):
            if i == 0:
                median = np.median(processed[obs]['residuals'][refraction[obs]['poly_flag']].flatten())
                offset = np.std(processed[obs]['residuals'][refraction[obs]['poly_flag']].flatten()) * 6
                y_limits = [median-offset, median+offset]

            color = colors[i][:-1]

            for order in range(0, input_data[obs]['n_orders']):

                # Workaround to damn stupid matplotlib error I didn't manage to solve
                ax.scatter(input_data[obs]['wave'][refraction[obs]['poly_flag']],
                           processed[obs]['residuals'][refraction[obs]['poly_flag']] + offset*i,
                           s=1, c=[color], alpha=0.50, zorder=2)

                ax.scatter(input_data[obs]['wave'][~refraction[obs]['poly_flag']],
                           processed[obs]['residuals'][~refraction[obs]['poly_flag']] + offset*i,
                           s=2, c='k', alpha=0.05, zorder=1)

            ax.axhline(offset*i, c='k', zorder=3)

        y_limits_offset = [min(y_limits[0] + offset * i, y_limits[0]),
                           max(y_limits[1] + offset * i, y_limits[1])]

        ax.set_ylim(y_limits_offset)
        ax.set_xlabel('$\lambda$ [$\AA$]')
        ax.legend(loc=3)
        ax.set_title('Night: {0:s} \n Residuals of the fit on ratio obs/master'.format(night))

        sm = plt.cm.ScalarMappable(cmap=color_cmap, norm=color_norm)
        sm.set_array([])  # You have to set a dummy-array for this to work...
        cbar = plt.colorbar(sm, cax=cbax1)
        cbar.set_label('BJD - 2450000.0')
        fig.subplots_adjust(wspace=0.05, hspace=0.4)
        plt.show()

        """ 
        PLOT: corrected e2ds
        """

        fig = plt.figure(figsize=(12, 6))
        gs = GridSpec(1, 2, width_ratios=[50, 1])
        ax = plt.subplot(gs[0, 0])
        cbax1 = plt.subplot(gs[:, 1])

        for i, obs in enumerate(lists['observations']):
            color = [colors[i][:-1]]

            ax.scatter(input_data[obs]['wave'][refraction[obs]['poly_flag']],
                       processed[obs]['e2ds_corrected'][refraction[obs]['poly_flag']]/processed[obs]['rescaling'],
                       s=2, c=color, alpha=0.10)

            ax.scatter(input_data[obs]['wave'][~refraction[obs]['poly_flag']],
                       processed[obs]['e2ds_corrected'][~refraction[obs]['poly_flag']]/processed[obs]['rescaling'],
                       s=2, c='k', alpha=0.05)

            #for order in range(0, np.size(input_data[obs]['wave'][:, 0])):
            #
            #    ax.plot(input_data[obs]['wave'][order, :],
            #            refraction[obs]['polyfit_e2ds'][order, :],
            #            c=color_array, lw=1)

        ax.set_xlabel('$\lambda$ [$\AA$]')
        ax.set_title('Night: {0:s} \n Corrected and rescaled e2ds spectra'.format(night))

        sm = plt.cm.ScalarMappable(cmap=color_cmap, norm=color_norm)
        sm.set_array([])  # You have to set a dummy-array for this to work...
        cbar = plt.colorbar(sm, cax=cbax1)
        cbar.set_label('BJD - 2450000.0')
        fig.subplots_adjust(wspace=0.05, hspace=0.4)
        plt.show()