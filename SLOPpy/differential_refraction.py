from __future__ import print_function, division
from SLOPpy.subroutines.common import *
from SLOPpy.subroutines.spectral_subroutines import *
from SLOPpy.subroutines.fit_subroutines import *
from SLOPpy.subroutines.io_subroutines import *
from SLOPpy.subroutines.shortcuts import *
from SLOPpy.subroutines.plot_subroutines import *
from SLOPpy.differential_refraction_preparation import compute_differential_refraction_preparation

__all__ = ["compute_differential_refraction",
           "plot_differential_refraction",
           "compute_differential_refraction_update",
           "plot_differential_refraction_update"]


def compute_differential_refraction_update(config_in):
    compute_differential_refraction(config_in, append_name='update')


def plot_differential_refraction_update(config_in, night_input=''):
    plot_differential_refraction(config_in, night_input, append_name='update')


def compute_differential_refraction(config_in, append_name=None):

    if append_name:
        subroutine_name = 'differential_refraction_' + append_name
        filename = 'refraction_' + append_name
    else:
        subroutine_name = 'differential_refraction'
        filename = 'refraction'

    compute_differential_refraction_preparation(config_in, append_name)

    night_dict = from_config_get_nights(config_in)

    for night in night_dict:

        try:
            refraction = load_from_cpickle(filename, config_in['output'], night)
            print("{0:45s} Night:{1:15s}   {2:s}".format(subroutine_name, night, 'Retrieved'))
            continue
        except:
            print("{0:45s} Night:{1:15s}   {2:s}".format(subroutine_name, night, 'Computing'))
            print()

        refraction_dict = from_config_refraction(config_in, night)

        """ Retrieving the list of observations"""
        lists = load_from_cpickle('lists', config_in['output'], night)

        """ Retrieving input and calibration data """
        calib_data = load_from_cpickle('calibration_fibA', config_in['output'], night)
        input_data = retrieve_observations(config_in['output'], night, lists['observations'],
                                           use_refraction=False)
        observational_pams = load_from_cpickle('observational_pams', config_in['output'], night)

        preparation = load_from_cpickle(filename + '_preparation', config_in['output'], night)

        defined_reference = night_dict[night]['refraction'].get('reference', False)

        if defined_reference:
            reference = load_from_cpickle(filename + '_reference', config_in['output'])
            preparation['coadd']['wave'] = reference['wave']
            preparation['coadd']['step'] = reference['step']
            preparation['coadd']['rescaled'] = reference['rescaled']

            if not preparation.get('absolute_SRF', False):
                print("  Observations and reference spectra are in different reference system ")
                quit()

        processed = {
            'subroutine': subroutine_name,
            'coadd': {
                'wave': preparation['coadd']['wave'],
                'step': preparation['coadd']['step'],
                'size': preparation['coadd']['size']
            }
        }

        refraction = {
            'subroutine': subroutine_name,
            'wave': preparation['coadd']['wave'],
            'binned': {}
        }

        refraction['binned']['wave'] = np.arange(preparation['coadd']['wave'][0],
                                                 preparation['coadd']['wave'][-11],
                                                 20.*preparation['coadd']['step'][0], dtype=np.double)

        refraction['binned']['size'] = np.size(refraction['binned']['wave'])
        refraction['binned']['step'] = np.ones(refraction['binned']['size'], dtype=np.double) \
            * 20. * preparation['coadd']['step'][0]

        if refraction_dict['approach'] == 'full_spectrum':
            print("  Differential refraction performed over the full spectrum")
        elif refraction_dict['approach'] == 'individual_order':
            print("   Differential refraction performed order-by-order ")
        else:
            raise ValueError("ERROR: fitting approach for differential refraction not implemented")

        if refraction_dict['method'] == 'spline':
            print("  Modelling performed with spline")
            print("  Spline order for differential refraction fit: ",
                  refraction_dict['fit_order'])
            print("  Knots spacing (in Angstrom) for refraction fit: ",
                  refraction_dict['knots_spacing'])
        elif refraction_dict['method'] == 'polynomial':
            print("  Modelling performed  with polynomials")
            print("  Chebyshev polynomial order for differential refraction fit: ",
                  refraction_dict['fit_order'])
        else:
            raise ValueError("ERROR: fitting method for differential refraction not implemented")

        print("  Number of iterations: ",
              refraction_dict['fit_iters'])
        print()

        """ Now each observation is divided by the reference spectrum, after being doppler-shifted to the observer RF
            The result is then used to model the flux variation
        """
        approach = refraction_dict.get('approach', 'full_spectrum')

        if approach == 'full_spectrum':

            for obs in lists['observations']:
                print("  Division by reference spectrum and fit of the flux variation: ", obs)

                refraction[obs] = {}
                processed[obs] = {}

                if preparation['absolute_SRF']:
                    rv_shift = observational_pams[obs]['rv_shift_ORF2SRF']
                else:
                    rv_shift = observational_pams[obs]['rv_shift_ORF2SRF_mod']

                processed[obs]['ratio'] = preparation[obs]['rescaled'] / preparation['coadd']['rescaled']

                processed[obs]['ratio_err'] = processed[obs]['ratio'] * \
                    np.sqrt((preparation[obs]['rescaled_err']/preparation[obs]['rescaled'])**2
                            + (preparation['coadd']['rescaled_err']/preparation['coadd']['rescaled'])**2)

                refraction[obs]['fit_flag'] = (processed[obs]['ratio'] > 0.01)

                if refraction_dict['method'] == 'spline':
                    for n_iter in range(0, refraction_dict['fit_iters']):
                        wave = processed['coadd']['wave'][refraction[obs]['fit_flag']]

                        """ picking the number of knots """
                        nknots = ((np.amax(wave) - np.amin(wave)) / refraction_dict['knots_spacing'])
                        """ picking the indices of the knots"""
                        idx_knots = (np.arange(1, len(wave) - 1, (len(wave) - 2.) / nknots)).astype('int')
                        """ passing from indices to knots values """
                        processed[obs]['knots'] = wave[idx_knots]

                        refraction[obs]['coeff'] = \
                            sci_int.splrep(
                                processed['coadd']['wave'][refraction[obs]['fit_flag']],
                                processed[obs]['ratio'][refraction[obs]['fit_flag']],
                                task=-1,
                                k=refraction_dict['fit_order'],
                                t=processed[obs]['knots'])

                        refraction[obs]['fit_s1d'] = sci_int.splev(processed['coadd']['wave'], refraction[obs]['coeff'])
                        processed[obs]['residuals'] = processed[obs]['ratio'] - refraction[obs]['fit_s1d']

                        if n_iter < refraction_dict['fit_iters'] - 1:
                            std = np.std(processed[obs]['residuals'])
                            refraction[obs]['fit_flag'] = (refraction[obs]['fit_flag']) \
                                & (np.abs(processed[obs]['residuals']) <
                                   refraction_dict['fit_sigma'] * std)

                elif refraction_dict['method'] == 'polynomial':

                    refraction[obs]['fit_flag'][:50] = False
                    refraction[obs]['fit_flag'][-50:] = False

                    for n_iter in range(0, refraction_dict['fit_iters']):

                        refraction[obs]['coeff'] = np.polynomial.chebyshev.chebfit(
                            processed['coadd']['wave'][refraction[obs]['fit_flag']],
                            processed[obs]['ratio'][refraction[obs]['fit_flag']],
                            refraction_dict['fit_order'])

                        refraction[obs]['fit_s1d'] = \
                            np.polynomial.chebyshev.chebval(processed['coadd']['wave'], refraction[obs]['coeff'])

                        processed[obs]['residuals'] = processed[obs]['ratio'] - refraction[obs]['fit_s1d']

                        if n_iter < refraction_dict['fit_iters'] - 1:
                            std = np.std(processed[obs]['residuals'])
                            refraction[obs]['fit_flag'] = (refraction[obs]['fit_flag']) \
                                & (np.abs(processed[obs]['residuals']) <
                                   refraction_dict['fit_sigma'] * std)

                """ Going back to the observer RF and rebinning the polynomial fit into the observed orders """
                refraction[obs]['fit_e2ds'] = \
                    rebin_1d_to_2d(processed['coadd']['wave'],
                                   input_data['coadd']['step'],
                                   refraction[obs]['fit_s1d'],
                                   input_data[obs]['wave'],
                                   input_data[obs]['step'],
                                   rv_shift=-rv_shift,
                                   preserve_flux=False,
                                   )

                """ Zero or negative values are identified, flagged and substituted with another value """
                refraction[obs]['null'] = np.zeros([input_data[obs]['n_orders'], input_data[obs]['n_pixels']],
                                                   dtype=bool)
                for order in range(0, observational_pams['n_orders']):
                    refraction[obs]['fit_e2ds'][order, :], _, refraction[obs]['null'][order, :] = \
                        replace_values_errors_with_interpolation_1d(refraction[obs]['fit_e2ds'][order, :],
                                                                    refraction[obs]['fit_e2ds'][order, :],
                                                                    less_than=0.001)

                processed[obs]['flux_rebinned_stellarRF_corrected'] = preparation[obs]['flux_rebinned_stellarRF'] \
                    / refraction[obs]['fit_s1d']
                refraction[obs]['binned_residuals'] = \
                    rebin_1d_to_1d(processed['coadd']['wave'],
                                   processed['coadd']['step'],
                                   processed[obs]['residuals'],
                                   refraction['binned']['wave'],
                                   refraction['binned']['step'],
                                   preserve_flux=False)

        elif approach == 'individual_order':

            for obs in lists['observations']:

                print("  Division by reference spectrum and fit of the flux variation: ", obs)

                refraction[obs] = {}
                processed[obs] = {}

                if preparation['absolute_SRF']:
                    rv_shift = observational_pams[obs]['rv_shift_ORF2SRF']
                else:
                    rv_shift = observational_pams[obs]['rv_shift_ORF2SRF_mod']

                """ Going back to the observer RF and rebinning the spectrum into the observed orders """
                preserve_flux = input_data[obs].get('absolute_flux', True)
                # processed[obs]['master_flux'] = \
                processed[obs]['master_flux'] = \
                    rebin_1d_to_2d(preparation['coadd']['wave'],
                                   preparation['coadd']['step'],
                                   preparation['coadd']['rescaled'],
                                   input_data[obs]['wave'],
                                   input_data[obs]['step'],
                                   preserve_flux=preserve_flux,
                                   rv_shift=-rv_shift)

                # processed[obs]['master_ferr'] = \
                processed[obs]['master_ferr'] = \
                    rebin_1d_to_2d(preparation['coadd']['wave'],
                                   preparation['coadd']['step'],
                                   preparation['coadd']['rescaled_err'],
                                   input_data[obs]['wave'],
                                   input_data[obs]['step'],
                                   preserve_flux=preserve_flux,
                                   rv_shift=-rv_shift,
                                   is_error=True)

                """ Zero or negative values are identified, flagged and substituted with another value """
                processed[obs]['master_flux'], processed[obs]['master_ferr'], processed[obs]['master_null'] = \
                    replace_values_errors_with_interpolation_2d(processed[obs]['master_flux'],
                                                                processed[obs]['master_ferr'],
                                                                less_than=0.001)

                # processed[obs]['ratio'] = preparation[obs]['rescaled_blazed'] / master_flux

                processed[obs]['ratio'] = input_data[obs]['e2ds'] \
                    / preparation[obs]['rescaling'] \
                    / (processed[obs]['master_flux'] * calib_data['blaze'])

                processed[obs]['ratio_err'] = processed[obs]['ratio'] * \
                    np.sqrt(preparation[obs]['rescaling']/input_data[obs]['e2ds']
                            + (processed[obs]['master_ferr']/processed[obs]['master_flux'])**2)

                refraction[obs]['fit_e2ds'] = np.zeros([input_data[obs]['n_orders'], input_data[obs]['n_pixels']])
                processed[obs]['residuals'] = np.zeros([input_data[obs]['n_orders'], input_data[obs]['n_pixels']])
                refraction[obs]['fit_flag'] = np.zeros([input_data[obs]['n_orders'], input_data[obs]['n_pixels']],
                                                       dtype=bool)

                for order in range(0, input_data[obs]['n_orders']):

                    order_coeff_name = 'order_' + repr(order)
                    refraction[obs]['fit_flag'][order, :] = (processed[obs]['ratio'][order, :] > 0.1)

                    if refraction_dict['method'] == 'spline':

                        for n_iter in range(0, refraction_dict['fit_iters']):

                            wave = input_data[obs]['wave'][order, refraction[obs]['fit_flag'][order, :]]

                            """ picking the number of knots """
                            nknots = ((np.amax(wave) - np.amin(wave)) / refraction_dict['knots_spacing'])
                            """ picking the indices of the knots"""
                            idx_knots = (np.arange(1, len(wave) - 1, (len(wave) - 2.) / nknots)).astype('int')
                            """ passing from indices to knots values """
                            refraction[obs]['knots'] = wave[idx_knots]

                            refraction[obs][order_coeff_name] = \
                                sci_int.splrep(
                                    input_data[obs]['wave'][order, refraction[obs]['fit_flag'][order, :]],
                                    processed[obs]['ratio'][order, refraction[obs]['fit_flag'][order, :]],
                                    task=-1,
                                    k=refraction_dict['fit_order'],
                                    t=refraction[obs]['knots'])

                            refraction[obs]['fit_e2ds'][order, :] = \
                                sci_int.splev(input_data[obs]['wave'][order, :],
                                              refraction[obs][order_coeff_name])

                            processed[obs]['residuals'][order, :] = processed[obs]['ratio'][order, :] \
                                - refraction[obs]['fit_e2ds'][order, :]

                            if n_iter < refraction_dict['fit_iters'] - 1:
                                std = np.std(processed[obs]['residuals'][order, :])
                                refraction[obs]['fit_flag'][order, :] = (refraction[obs]['fit_flag'][order, :]) \
                                    & (np.abs(processed[obs]['residuals'][order, :]) <
                                       refraction_dict['fit_sigma'] * std)

                    elif refraction_dict['method'] == 'polynomial':
                        refraction[obs]['fit_flag'][order, :50] = False
                        refraction[obs]['fit_flag'][order, -50:] = False

                        for n_iter in range(0, refraction_dict['fit_iters']):

                            refraction[obs][order_coeff_name] = np.polynomial.chebyshev.chebfit(
                                input_data[obs]['wave'][order, refraction[obs]['fit_flag'][order, :]],
                                processed[obs]['ratio'][order, refraction[obs]['fit_flag'][order, :]],
                                refraction_dict['fit_order'])

                            refraction[obs]['fit_e2ds'][order, :] = \
                                np.polynomial.chebyshev.chebval(input_data[obs]['wave'][order, :],
                                                                refraction[obs][order_coeff_name])

                            processed[obs]['residuals'][order, :] = processed[obs]['ratio'][order, :] \
                                - refraction[obs]['fit_e2ds'][order, :]

                            if n_iter < refraction_dict['fit_iters'] - 1:
                                std = np.std(processed[obs]['residuals'][order, :])
                                refraction[obs]['fit_flag'][order, :] = (refraction[obs]['fit_flag'][order, :]) \
                                    & (np.abs(processed[obs]['residuals'][order, :]) <
                                       refraction_dict['fit_sigma'] * std)

                e2ds_corrected = input_data[obs]['e2ds'] / refraction[obs]['fit_e2ds']
                e2ds_corrected_err = input_data[obs]['e2ds_err'] / refraction[obs]['fit_e2ds']

                preserve_flux = input_data[obs].get('absolute_flux', True)
                processed[obs]['flux_rebinned_stellarRF_corrected'] = \
                    rebin_2d_to_1d(input_data[obs]['wave'],
                                   input_data[obs]['step'],
                                   e2ds_corrected,
                                   calib_data['blaze'],
                                   processed['coadd']['wave'],
                                   input_data['coadd']['step'],
                                   preserve_flux=preserve_flux,
                                   rv_shift=rv_shift)

                processed[obs]['err_flux_rebinned_stellarRF_corrected'] = \
                    rebin_2d_to_1d(input_data[obs]['wave'],
                                   input_data[obs]['step'],
                                   e2ds_corrected_err,
                                   calib_data['blaze'],
                                   processed['coadd']['wave'],
                                   input_data['coadd']['step'],
                                   preserve_flux=preserve_flux,
                                   rv_shift=rv_shift,
                                   is_error=True)

                processed[obs]['flux_rebinned_stellarRF_corrected'], \
                    processed[obs]['err_flux_rebinned_stellarRF_corrected'], _ = \
                    replace_values_errors_with_interpolation_1d(processed[obs]['flux_rebinned_stellarRF_corrected'],
                                                                processed[obs]['err_flux_rebinned_stellarRF_corrected'],
                                                                less_than=0.001)

                refraction[obs]['binned_residuals'] = \
                    rebin_2d_to_1d(input_data[obs]['wave'],
                                   input_data[obs]['step'],
                                   processed[obs]['residuals'],
                                   np.ones_like(processed[obs]['residuals']),
                                   refraction['binned']['wave'],
                                   refraction['binned']['step'],
                                   rv_shift=0.0000,
                                   preserve_flux=False)

        else:
            print("   Please choose either full_spectrum or individual_order as preferred approach")
            quit()

        if not config_in['settings'].get('full_output', False):
            for obs in lists['observations']:

                del processed[obs]['ratio_err']

                try:
                    del processed[obs]['err_flux_rebinned_stellarRF_corrected']
                    del processed[obs]['master_flux']
                    del processed[obs]['master_ferr']
                    del processed[obs]['master_null']
                except:
                    pass

        if append_name:
            save_to_cpickle('refraction_' + append_name + '_processed', processed, config_in['output'], night)
            save_to_cpickle('refraction_' + append_name, refraction, config_in['output'], night)
        else:
            save_to_cpickle('refraction_processed', processed, config_in['output'], night)
            save_to_cpickle('refraction', refraction, config_in['output'], night)


def plot_differential_refraction(config_in, night_input='', append_name=None):
    night_dict = from_config_get_nights(config_in)

    if night_input == '':
        night_list = night_dict
    else:
        night_list = np.atleast_1d(night_input)

    for night in night_list:

        refraction_dict = from_config_refraction(config_in, night)

        """ Retrieving the list of observations"""
        lists = load_from_cpickle('lists', config_in['output'], night)

        """ Retrieving the observations"""
        input_data = retrieve_observations(config_in['output'], night, lists['observations'],
                                           use_refraction=False, use_telluric=False)
        observational_pams = load_from_cpickle('observational_pams', config_in['output'], night)

        try:
            """ Retrieving the analysis"""
            if append_name:
                processed = load_from_cpickle('refraction_' + append_name + '_processed', config_in['output'], night)
                preparation = load_from_cpickle('refraction_' + append_name + '_preparation', config_in['output'],
                                                night)
                refraction = load_from_cpickle('refraction_' + append_name, config_in['output'], night)
            else:
                processed = load_from_cpickle('refraction_processed', config_in['output'], night)
                preparation = load_from_cpickle('refraction_preparation', config_in['output'], night)
                refraction = load_from_cpickle('refraction', config_in['output'], night)
        except:
            print("  Failed in retrieving processed data")
            return

        approach = refraction_dict.get('approach', 'full_spectrum')

        colors_properties, colors_plot, colors_scatter = make_color_array_matplotlib3(lists, observational_pams)

        offset = 0.10

        y_limits = [0.8, 1.2]

        flag_e2ds = {}
        flag_coadd = {}
        for i, obs in enumerate(lists['observations']):

            shrink_factor = 4
            if input_data[obs]['n_orders'] > shrink_factor:
                factor = (input_data[obs]['n_orders'] * input_data[obs]['n_pixels']) \
                    // (input_data[obs]['n_pixels'] * shrink_factor)
                flag_e2ds[obs] = (np.random.choice(a=([False] * (factor-1)) + [True],
                                                   size=(input_data[obs]['n_orders'], input_data[obs]['n_pixels'])))

                flag_coadd[obs] = \
                    np.random.choice(a=([False] * factor) + [True],  size=input_data['coadd']['size'])

            else:
                flag_e2ds[obs] = np.ones([input_data[obs]['n_orders'], input_data[obs]['n_pixels']], dtype=bool)
                flag_coadd[obs] = np.ones(input_data['coadd']['size'], dtype=bool)

        fig = plt.figure(figsize=(12, 6))
        gs = GridSpec(2, 2, width_ratios=[50, 1])

        ax1 = plt.subplot(gs[0, 0])
        ax2 = plt.subplot(gs[1, 0], sharex=ax1, sharey=ax1)
        cbax1 = plt.subplot(gs[:, 1])

        for i, obs in enumerate(lists['observations']):

            """  We slim down the plot  """

            if i == 0:
                ax1.scatter(preparation['coadd']['wave'][flag_coadd[obs]],
                            preparation[obs]['flux_rebinned_stellarRF'][flag_coadd[obs]] /
                            preparation[obs]['rescaling'],
                            c=colors_scatter['mBJD'][obs], s=2, alpha=0.2, label='observation (SRF)')
            else:
                ax1.scatter(preparation['coadd']['wave'][flag_coadd[obs]],
                            preparation[obs]['flux_rebinned_stellarRF'][flag_coadd[obs]] /
                            preparation[obs]['rescaling'],
                            c=colors_scatter['mBJD'][obs], s=2, alpha=0.2)

            ax2.scatter(processed['coadd']['wave'][flag_coadd[obs]],
                        processed[obs]['flux_rebinned_stellarRF_corrected'][flag_coadd[obs]] /
                        preparation[obs]['rescaling'],
                        c=colors_scatter['mBJD'][obs], s=3, alpha=0.2)

        ax1.plot(preparation['coadd']['wave'], preparation['coadd']['rescaled'], c='k', lw=1, alpha=0.5,
                 label='reference spectrum')
        ax2.plot(preparation['coadd']['wave'], preparation['coadd']['rescaled'], c='k', lw=1, alpha=0.5)

        ax1.set_xlim(processed['coadd']['wave'][0], processed['coadd']['wave'][-1])
        ax1.set_ylim(y_limits)
        ax2.set_ylim(y_limits)
        ax1.legend(loc=1)
        ax1.set_title('Night: {0:s} \n Input spectra'.format(night))
        ax2.set_title('After differential refraction correction')
        ax2.set_xlabel('$\lambda$ [$\AA$]')

        sm = plt.cm.ScalarMappable(cmap=colors_properties['cmap'], norm=colors_properties['norm']['mBJD'])
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
                offset = np.std(processed[obs]['ratio'][refraction[obs]['fit_flag']].flatten()) * 6
                average = np.average(processed[obs]['ratio'][refraction[obs]['fit_flag']].flatten())
                y_limits = [average - offset, average + offset]

            if approach == 'full_spectrum':
                flag = flag_coadd[obs] & refraction[obs]['fit_flag']
                wave = refraction['wave']

            elif approach == 'individual_order':
                flag = flag_e2ds[obs] & refraction[obs]['fit_flag']
                wave = input_data[obs]['wave']

            ax.scatter(wave[flag],
                       processed[obs]['ratio'][flag] + offset * i,
                       c=colors_scatter['mBJD'][obs], s=1, alpha=0.50, zorder=2)

            ax.scatter(wave[~refraction[obs]['fit_flag']],
                       processed[obs]['ratio'][~refraction[obs]['fit_flag']] + offset * i,
                       c='k', s=2, alpha=0.1, zorder=1)

            for order in range(0, input_data[obs]['n_orders']):

                ax.plot(input_data[obs]['wave'][order, :],
                        refraction[obs]['fit_e2ds'][order, :] + offset * i,
                        c='k', lw=1, alpha=0.5, zorder=5)

        y_limits_offset = [min(y_limits[0] + offset * i, y_limits[0]),
                           max(y_limits[1] + offset * i, y_limits[1])]

        ax.set_ylim(y_limits_offset)
        ax.set_xlabel('$\lambda$ [$\AA$]')
        # ax.legend(loc=3)
        ax.set_title('Night: {0:s} \n Differential refraction correction - Fit of the ratio obs/master'.format(night))

        sm = plt.cm.ScalarMappable(cmap=colors_properties['cmap'], norm=colors_properties['norm']['mBJD'])
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

        approach = refraction_dict.get('approach', 'full_spectrum')

        for i, obs in enumerate(lists['observations']):

            if i == 0:
                median = np.median(processed[obs]['residuals'][refraction[obs]['fit_flag']].flatten())
                offset = np.std(processed[obs]['residuals'][refraction[obs]['fit_flag']].flatten()) * 6
                y_limits = [median - offset, median + offset]

            if approach == 'full_spectrum':
                flag = flag_coadd[obs] & refraction[obs]['fit_flag']
                wave = refraction['wave']

            elif approach == 'individual_order':
                flag = flag_e2ds[obs] & refraction[obs]['fit_flag']
                wave = input_data[obs]['wave']

            ax.scatter(wave[flag],
                       processed[obs]['residuals'][flag] + offset * i,
                       c=colors_scatter['mBJD'][obs], s=1, alpha=0.50, zorder=2)

            ax.scatter(wave[~refraction[obs]['fit_flag']],
                       processed[obs]['residuals'][~refraction[obs]['fit_flag']] + offset * i,
                       c='k', s=2, alpha=0.1, zorder=1)

            ax.axhline(offset * i, c='k', zorder=3)

        y_limits_offset = [min(y_limits[0] + offset * i, y_limits[0]),
                           max(y_limits[1] + offset * i, y_limits[1])]

        ax.set_ylim(y_limits_offset)
        ax.set_xlabel('$\lambda$ [$\AA$]')
        # ax.legend(loc=3)
        ax.set_title(
            'Night: {0:s} \n Differential refraction correction - Residuals of the fit on ratio obs/master'.format(
                night))

        sm = plt.cm.ScalarMappable(cmap=colors_properties['cmap'], norm=colors_properties['norm']['mBJD'])
        sm.set_array([])  # You have to set a dummy-array for this to work...
        cbar = plt.colorbar(sm, cax=cbax1)
        cbar.set_label('BJD - 2450000.0')
        fig.subplots_adjust(wspace=0.05, hspace=0.4)
        plt.show()

        continue

        """
        PLOT: corrected e2ds
        """

        fig = plt.figure(figsize=(12, 6))
        gs = GridSpec(1, 2, width_ratios=[50, 1])
        ax = plt.subplot(gs[0, 0])
        cbax1 = plt.subplot(gs[:, 1])

        for i, obs in enumerate(lists['observations']):

            e2ds_corrected = input_data[obs]['e2ds'] / refraction[obs]['fit_e2ds']

            ax.scatter(input_data[obs]['wave'],
                       e2ds_corrected / preparation[obs]['rescaling'],
                       c=colors_scatter['mBJD'][obs], s=2, alpha=0.20)

            # for order in range(0, np.size(input_data[obs]['wave'][:, 0])):
            #
            #    ax.plot(input_data[obs]['wave'][order, :],
            #            refraction[obs]['fit_e2ds'][order, :],
            #            c=color_array, lw=1)

        ax.set_xlabel('$\lambda$ [$\AA$]')
        ax.set_title('Night: {0:s} \n Rescaled e2ds spectra after differential refraction correction'.format(night))

        sm = plt.cm.ScalarMappable(cmap=colors_properties['cmap'], norm=colors_properties['norm']['mBJD'])
        sm.set_array([])  # You have to set a dummy-array for this to work...
        cbar = plt.colorbar(sm, cax=cbax1)
        cbar.set_label('BJD - 2450000.0')
        fig.subplots_adjust(wspace=0.05, hspace=0.4)
        plt.show()
