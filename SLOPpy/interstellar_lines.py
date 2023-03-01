from __future__ import print_function, division
from SLOPpy.subroutines.common import *
from SLOPpy.subroutines.spectral_subroutines import *
from SLOPpy.subroutines.fit_subroutines import *

from SLOPpy.subroutines.plot_subroutines import *
from SLOPpy.subroutines.shortcuts import *

__all__ = ["compute_interstellar_lines", "plot_interstellar_lines"]

subroutine_name = 'interstellar_lines'


# def plot_identify_stellar_lines(config_in)

def compute_interstellar_lines(config_in):

    night_dict = from_config_get_nights(config_in)

    instrument_dict = from_config_get_instrument(config_in)

    interstellar_lines = from_config_get_interstellar_lines(config_in)

    shared_data = load_from_cpickle('shared', config_in['output'])

    if not interstellar_lines:
        return

    for night in night_dict:
        print()
        print("compute_interstellar_lines                 Night: ", night)

        try:
            interstellar = load_from_cpickle('interstellar_lines_processed', config_in['output'], night)
            skip_lineselection = True
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

        instrument = night_dict[night]['instrument']

        processed = {
            'subroutine': subroutine_name,
            'line_rebin': {},
            'line_shift': {}
        }

        interstellar = {
            'subroutine': subroutine_name,
        }

        for obs in lists['observations']:
            processed[obs] = {
                'n_orders': input_data[obs]['n_orders'],
                'n_pixels': input_data[obs]['n_pixels']
            }

            interstellar[obs] = {
                'wave_BRF': shift_wavelength_array(input_data[obs]['wave'], observational_pams[obs]['rv_shift_ORF2BRF']),
                'correction':  np.ones(np.shape(input_data[obs]['wave']))
            }

            """ for plotting purpose only"""
            processed[obs]['flux'] = input_data[obs]['e2ds'] / calib_data['blaze'] / input_data[obs]['step']
            processed[obs]['flux_err'] = np.sqrt(input_data[obs]['e2ds']) / calib_data['blaze'] / input_data[obs][
                'step']

        for line_name, line in interstellar_lines.items():

            processed[line_name] = {
                'line_rebin': {},
                'poly_coeff': {},
                'normalized': {},
                'line_shift': {
                    'selected_points': []
                }
            }

            processed[line_name]['min_wave'] = max(shared_data['coadd']['wavelength_range'][0], line[0] - line[2]*2)
            processed[line_name]['max_wave'] = min(shared_data['coadd']['wavelength_range'][1], line[0] + line[2]*2)

            processed[line_name]['wave'] = np.arange(processed[line_name]['min_wave'],
                                                     processed[line_name]['max_wave'],
                                                     instrument_dict[instrument]['wavelength_step'])

            processed[line_name]['size'] = np.size(processed[line_name]['wave'], axis=0)
            processed[line_name]['step'] = np.ones(processed[line_name]['size'])\
                * instrument_dict[instrument]['wavelength_step']

            processed[line_name]['correction'] = np.ones(processed[line_name]['size'])

            for obs in lists['observations']:

                preserve_flux = input_data[obs].get('absolute_flux', True)

                """ shifting a chunk of the spectra to the Solar System Barycenter reference """
                processed[line_name]['line_rebin'][obs] = \
                    rebin_2d_to_1d(input_data[obs]['wave'],
                                   input_data[obs]['step'],
                                   input_data[obs]['e2ds'],
                                   calib_data['blaze'],
                                   processed[line_name]['wave'],
                                   processed[line_name]['step'],
                                   preserve_flux=preserve_flux,
                                   rv_shift=observational_pams[obs]['rv_shift_ORF2BRF'])

                argmin_sel = np.argmin(processed[line_name]['line_rebin'][obs][5:-5]) + 5
                wave_sel = processed[line_name]['wave'][argmin_sel]

                processed[line_name]['line_shift']['selected_points'].append(wave_sel)

            processed[line_name]['line_shift']['position'] = np.median(
                processed[line_name]['line_shift']['selected_points'])
            processed[line_name]['line_shift']['delta_lambda'] = processed[line_name]['line_shift']['position'] - line[0]

            try:
                if interstellar_lines['automatic'] == True:
                    interstellar[line_name] = processed[line_name]['line_shift']['position']
                else:
                    interstellar[line_name] = line[0]
            except:
                interstellar[line_name] = line[0]

            """ selection of spectral range for continuum normalization and interstellar line modelling"""

            processed[line_name]['wavelength_selection'] = \
                (np.abs(processed[line_name]['wave'] - interstellar[line_name]) < line[1])

            processed[line_name]['continuum_selection'] = \
                (~processed[line_name]['wavelength_selection']) \
                & (np.abs(processed[line_name]['wave'] - interstellar[line_name]) < line[2])

            processed[line_name]['interstellar_selection'] = \
                (processed[line_name]['wavelength_selection'] | processed[line_name]['continuum_selection'])

            spline_wave_points = []
            spline_norm_points = []

            # TODO
            #! 1) Rescale by median each observation and collect all the value
            #! 2) Perform a continuum normalization on the collect values, with
            #! iterative sigma-clipping
            #! 3) Perform a spline / gaussian fit of the spectral line

            for obs in lists['telluric']:

                # sel1 = (np.abs(processed[line_name]['wave'] - interstellar[line_name]) < line[1])
                # sel2 = (~sel1) & (np.abs(processed[line_name]['wave'] - interstellar[line_name]) < line[2])
                # sel3 = (sel1 | sel2)

                """ normalization around the interstellar line """
                processed[line_name]['poly_coeff'][obs] = \
                    np.polyfit(processed[line_name]['wave'][processed[line_name]['continuum_selection']],
                               processed[line_name]['line_rebin'][obs][processed[line_name]['continuum_selection']],
                               2)

                processed[line_name]['normalized'][obs] = \
                    processed[line_name]['line_rebin'][obs][processed[line_name]['interstellar_selection']] \
                    / np.polyval(processed[line_name]['poly_coeff'][obs],
                                 processed[line_name]['wave'][processed[line_name]['interstellar_selection']])

                spline_wave_points.extend(processed[line_name]['wave'][processed[line_name]['interstellar_selection']])
                spline_norm_points.extend(processed[line_name]['normalized'][obs])

            """ sorting the array to avoid problems with the spline function"""
            spline_sorting_index = np.argsort(spline_wave_points)
            spline_wave_points = np.asarray(spline_wave_points)[spline_sorting_index]
            spline_norm_points = np.asarray(spline_norm_points)[spline_sorting_index]

            processed[line_name]['spline_eval'], \
                processed[line_name]['spline_coeff'], \
                processed[line_name]['spline_knots'] = \
                compute_spline(spline_wave_points, spline_norm_points, 0.08, knot_order=3)

            processed[line_name]['correction'][processed[line_name]['wavelength_selection']] = \
                sci_int.splev(processed[line_name]['wave'][processed[line_name]['wavelength_selection']],
                              processed[line_name]['spline_coeff'])

            for obs in lists['observations']:

                interstellar[obs]['wavelength_selection'] = \
                    (np.abs(interstellar[obs]['wave_BRF']-interstellar[line_name]) < line[1])
                interstellar[obs]['continuum_selection'] = \
                    (~interstellar[obs]['wavelength_selection']) \
                    & (np.abs(interstellar[obs]['wave_BRF']-interstellar[line_name]) < line[2])
                interstellar[obs]['interstellar_selection'] = \
                    (interstellar[obs]['wavelength_selection'] | interstellar[obs]['continuum_selection'])

                interstellar[obs]['correction'][interstellar[obs]['wavelength_selection']] = \
                    sci_int.splev(interstellar[obs]['wave_BRF'][interstellar[obs]['wavelength_selection']],
                                  processed[line_name]['spline_coeff'])

        save_to_cpickle('interstellar_lines_processed', processed, config_in['output'], night)
        save_to_cpickle('interstellar_lines', interstellar, config_in['output'], night)


def plot_interstellar_lines(config_in, night_input=''):
    import matplotlib.pyplot as plt

    night_dict = from_config_get_nights(config_in)
    instrument_dict = from_config_get_instrument(config_in)
    system_dict = from_config_get_system(config_in)

    interstellar_lines = from_config_get_interstellar_lines(config_in)

    if not interstellar_lines:
        return

    if night_input == '':
        night_list = night_dict
    else:
        night_list = np.atleast_1d(night_input)

    for night in night_list:

        print("plot_interstellar_lines                    Night: ", night)

        """ Retrieving the list of observations"""
        lists = load_from_cpickle('lists', config_in['output'], night)
        observational_pams = load_from_cpickle('observational_pams', config_in['output'], night)

        """ Retrieving the analysis"""
        try:
            processed = load_from_cpickle('interstellar_lines_processed', config_in['output'], night)
            interstellar = load_from_cpickle('interstellar_lines', config_in['output'], night)
        except:
            print()
            print('No interstellar correction, no plots')
            continue

        colors_properties, colors_plot, colors_scatter = make_color_array_matplotlib3(lists, observational_pams)

        # fig = plt.figure(figsize=(12, 6))
        # gs = GridSpec(2, 2, width_ratios=[50, 1])
        #
        # ax1 = plt.subplot(gs[0, 0])
        # ax2 = plt.subplot(gs[1, 0], sharex=ax1)
        # cbax1 = plt.subplot(gs[:, 1])
        fig, gs, cbax1, ax1, ax2 = grid_2plot()

        for i, obs in enumerate(lists['observations']):

            """rescaling"""
            processed[obs]['flux_rescaling'], processed[obs]['flux_rescaled'], processed[obs]['flux_rescaled_err'] = \
                perform_rescaling(interstellar[obs]['wave_BRF'],
                                  processed[obs]['flux'],
                                  processed[obs]['flux_err'],
                                  observational_pams['wavelength_rescaling'])

            if i == 0:
                ax1.scatter(interstellar[obs]['wave_BRF'], processed[obs]['flux_rescaled'],
                            c=colors_scatter['mBJD'][obs], s=1, alpha=0.5, label='observations (BRF)')
            else:
                ax1.scatter(interstellar[obs]['wave_BRF'], processed[obs]['flux_rescaled'],
                            c=colors_scatter['mBJD'][obs], s=1, alpha=0.5)

            # ax1.plot(interstellar['wave'], interstellar['correction'], c='black')

            ax2.scatter(interstellar[obs]['wave_BRF'], processed[obs]['flux_rescaled']/interstellar[obs]['correction'],
                        c=colors_scatter['mBJD'][obs], s=1, alpha=0.5)

        for line_name, line in interstellar_lines.items():

            ax1.axvline(interstellar[line_name]-line[1], c='b')
            ax1.axvline(interstellar[line_name]+line[1], c='b')
            ax1.axvline(interstellar[line_name]-line[2], c='g')
            ax1.axvline(interstellar[line_name]+line[2], c='g')

            ax2.axvline(interstellar[line_name]-line[1], c='b')
            ax2.axvline(interstellar[line_name]+line[1], c='b')
            ax2.axvline(interstellar[line_name]-line[2], c='g')
            ax2.axvline(interstellar[line_name]+line[2], c='g')

            # ax1.plot(processed[line_name]['wave'], processed[line_name]['flux_rescaled'], c='b')

            try:
                wave_min = min(wave_min, interstellar[line_name])
                wave_max = max(wave_max, interstellar[line_name])
                range_max = max(range_max, line[2])
            except:
                wave_min = interstellar[line_name]
                wave_max = interstellar[line_name]
                range_max = line[2]

        ax1.set_title('Night: {0:s} \n Input spectra'.format(night))
        ax1.set_xlim(wave_min-2*range_max, wave_max+2*range_max)
        ax1.legend(loc=1)

        ax2.set_title('After interstellar line correction')
        ax2.set_xlabel('$\lambda$ [$\AA$]')

        sm = plt.cm.ScalarMappable(cmap=colors_properties['cmap'], norm=colors_properties['norm']['mBJD'])
        sm.set_array([])  # You have to set a dummy-array for this to work...
        cbar = plt.colorbar(sm, cax=cbax1)
        cbar.set_label('BJD - 2450000.0')
        # fig.subplots_adjust(wspace=0.05, hspace=0.4)
        plt.show()
