from __future__ import print_function, division
from SLOPpy.subroutines.common import *
from SLOPpy.subroutines.spectral_subroutines import *
from SLOPpy.subroutines.fit_subroutines import *

from SLOPpy.subroutines.plot_subroutines import *
from SLOPpy.subroutines.shortcuts import *

__all__ = ["compute_interstellar_lines", "plot_interstellar_lines"]

subroutine_name = 'interstellar_lines'


#def plot_identify_stellar_lines(config_in)

def compute_interstellar_lines(config_in):

    night_dict = from_config_get_nights(config_in)

    interstellar_lines = from_config_get_interstellar_lines(config_in)

    if not interstellar_lines:
        return

    for night in night_dict:

        print()
        print("compute_interstellar_lines                 Night: ", night)

        try:
            interstellar = load_from_cpickle('interstellar_lines', config_in['output'], night)
            continue
        except:
            print()
            print("  No interstellar correction file found, computing now ")

        """ Retrieving the list of observations"""
        lists = load_from_cpickle('lists', config_in['output'], night)

        calib_data = load_from_cpickle('calibration_fibA', config_in['output'], night)
        input_data = retrieve_observations(config_in['output'], night, lists['observations'])
        observational_pams = load_from_cpickle('observational_pams', config_in['output'], night)

        processed = {
            'subroutine': subroutine_name
        }

        interstellar = {
            'subroutine': subroutine_name,
        }

        import matplotlib.pyplot as plt

        for obs in lists['observations']:

            processed[obs] = {
                'n_orders': input_data[obs]['n_orders'],
                'n_pixels': input_data[obs]['n_pixels']
            }
            interstellar[obs] = {}

            """ for plotting purpose only"""
            processed[obs]['wave'] = input_data[obs]['wave']

            processed[obs]['flux'] = input_data[obs]['e2ds']/calib_data['blaze']/input_data[obs]['step']
            processed[obs]['flux_err'] = np.sqrt(input_data[obs]['e2ds'])/calib_data['blaze']/input_data[obs]['step']

            if obs in lists['telluric']:

                try:
                    interstellar['flux_total'] += processed[obs]['flux']
                    interstellar['flux_total_err'] += processed[obs]['flux_err']**2
                except:
                    interstellar['wave'] = input_data[obs]['wave']
                    interstellar['flux_total'] = processed[obs]['flux'][:,:]
                    interstellar['flux_total_err'] = processed[obs]['flux_err']**2

        interstellar['flux_total_err'] = np.sqrt(interstellar['flux_total_err'])

        """ Zero or negative values are identified, flagged and substituted with another value """
        interstellar['flux_total'], interstellar['flux_total_err'], interstellar['null'] = \
                replace_values_errors(interstellar['flux_total'], interstellar['flux_total_err'], 0.0001)

        """rescaling"""
        interstellar['flux_rescaling'], interstellar['flux_rescaled'],interstellar['flux_rescaled_err']  = \
                perform_rescaling(interstellar['wave'],
                                  interstellar['flux_total'],
                                  interstellar['flux_total_err'],
                                  observational_pams['wavelength_rescaling'])

        interstellar['correction'] = np.ones(np.shape(interstellar['wave']))

        for line_name, line in interstellar_lines.items():

            interstellar[line_name] = {}

            sel1 = (np.abs(interstellar['wave']-line[0])<line[1])
            sel2 = (~sel1) & (np.abs(interstellar['wave']-line[0])<line[2])
            sel3 = (sel1 | sel2)

            poly_coeff = np.polyfit(interstellar['wave'][sel2], interstellar['flux_rescaled'][sel2], 2)

            normalized = interstellar['flux_rescaled'][sel3]/np.polyval(poly_coeff, interstellar['wave'][sel3])

            interstellar[line_name]['spline_eval'], \
            interstellar[line_name]['spline_coeff'], \
            interstellar[line_name]['spline_knots'] = \
                compute_spline(interstellar['wave'][sel3], normalized, 0.04)

            interstellar['correction'][sel1] = sci_int.splev(interstellar['wave'][sel1], interstellar[line_name]['spline_coeff'])

            interstellar[line_name]['sel1'] = sel1
            interstellar[line_name]['sel2'] = sel2
            interstellar[line_name]['sel3'] = sel3

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

    if night_input=='':
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

        colors, cmap, line_colors = make_color_array(lists, observational_pams)

        fig = plt.figure(figsize=(12, 6))
        gs = GridSpec(2, 2, width_ratios=[50, 1])

        ax1 = plt.subplot(gs[0, 0])
        ax2 = plt.subplot(gs[1, 0], sharex=ax1)
        cbax1 = plt.subplot(gs[:, 1])

        for i, obs in enumerate(lists['observations']):

            """rescaling"""
            processed[obs]['flux_rescaling'], processed[obs]['flux_rescaled'], processed[obs]['flux_rescaled_err']  = \
                perform_rescaling(interstellar['wave'],
                                  processed[obs]['flux'],
                                  processed[obs]['flux_err'],
                                  observational_pams['wavelength_rescaling'])
            ax1.scatter(interstellar['wave'], processed[obs]['flux_rescaled'],
                            s=1, c=line_colors[i])

            #ax1.plot(interstellar['wave'], interstellar['correction'], c='black')

            ax2.scatter(interstellar['wave'], processed[obs]['flux_rescaled']/interstellar['correction'],
                        s=1, c=line_colors[i])

        for line_name, line in interstellar_lines.items():

            #ax1.axvline(line[0], c='k')
            ax1.axvline(line[0]-line[1], c='b')
            ax1.axvline(line[0]+line[1], c='b')
            ax1.axvline(line[0]-line[2], c='g')
            ax1.axvline(line[0]+line[2], c='g')

            #ax2.axvline(line[0], c='k')
            ax2.axvline(line[0]-line[1], c='b')
            ax2.axvline(line[0]+line[1], c='b')
            ax2.axvline(line[0]-line[2], c='g')
            ax2.axvline(line[0]+line[2], c='g')

            try:
                wave_min = min(wave_min, line[0])
                wave_max = max(wave_max, line[0])
                range_max = max(range_max, line[2])
            except:
                wave_min = line[0]
                wave_max = line[0]
                range_max = line[2]

        #ax1.legend(loc=3)
        ax1.set_title('Night: ' + night)
        ax1.set_xlim(wave_min-2*range_max, wave_max+2*range_max)

        ax2.set_xlabel('$\lambda$ [$\AA$]')
        ax1.set_xlabel('$\lambda$ [$\AA$]')

        sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=colors[0], vmax=colors[-1]))
        sm.set_array([])  # You have to set a dummy-array for this to work...
        cbar = plt.colorbar(sm, cax=cbax1)
        cbar.set_label('BJD - 2450000.0')
        fig.subplots_adjust(wspace=0.05, hspace=0.4)
        plt.show()