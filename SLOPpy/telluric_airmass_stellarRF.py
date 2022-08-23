from __future__ import print_function, division
from SLOPpy.subroutines.common import *
from SLOPpy.subroutines.spectral_subroutines import *
from SLOPpy.subroutines.io_subroutines import *
from SLOPpy.subroutines.fit_subroutines import *

from SLOPpy.subroutines.plot_subroutines import *
from SLOPpy.subroutines.shortcuts import *

__all__ = ["compute_telluric_airmass_stellarRF",
           "plot_telluric_airmass_stellarRF",
           "compute_telluric_airmass_reference_stellarRF",
           "plot_telluric_airmass_reference_stellarRF"]

subroutine_name = 'telluric_airmass_stellarRF'


def compute_telluric_airmass_stellarRF(config_in):

    compute_telluric_stellarRF(config_in,
                                use_reference_airmass=False,
                                subroutine_name='telluric_airmass_stellarRF')


def compute_telluric_airmass_reference_stellarRF(config_in):
    compute_telluric_stellarRF(config_in,
                                use_reference_airmass=True,
                                subroutine_name='telluric_airmass_reference_stellarRF')


def plot_telluric_airmass_reference_stellarRF(config_in, night_input):
    """ Alias to simplify the configuration file"""
    plot_telluric_airmass_stellarRF(config_in, night_input)


def compute_telluric_stellarRF(config_in, **kwargs):

    night_dict = from_config_get_nights(config_in)

    for night in night_dict:

        print()
        print("compute_telluric_airmass_stellarRF         Night: ", night)

        try:
            telluric = load_from_cpickle('telluric', config_in['output'], night)
            continue
        except:
            print()
            print("  No telluric correction file found, computing now ")

        """ Retrieving the list of observations"""
        lists = load_from_cpickle('lists', config_in['output'], night)

        calib_data = load_from_cpickle('calibration_fibA', config_in['output'], night)
        input_data = retrieve_observations(config_in['output'], night, lists['observations'], use_telluric=False)
        observational_pams = load_from_cpickle('observational_pams', config_in['output'], night)

        processed = {
            'subroutine': subroutine_name
        }

        telluric = {
            'subroutine': kwargs['subroutine_name'],
            'reference_frame': 'stellar'
        }

        """ Reference airmass for iterative correction of airmass"""
        if kwargs['use_reference_airmass']:
            airmass_temp = np.zeros(lists['n_transit_in'])
            for n_obs, obs in enumerate(lists['transit_in']):
                # This is to ensure that airmass, berv and rvc are associated to the correct spectra
                airmass_temp[n_obs] = input_data[obs]['AIRMASS']
            processed['airmass_ref'] = np.average(airmass_temp)
        else:
            processed['airmass_ref'] = 0.000


        # There must be a more elegant way to do this, but I'm, not aware of it
        for obs in lists['observations']:

            processed[obs] = {
                'n_orders': input_data[obs]['n_orders'],
                'n_pixels': input_data[obs]['n_pixels']
            }
            telluric[obs] = {}

            """ for plotting purpose only"""
            processed[obs]['wave'] = input_data[obs]['wave']
            processed[obs]['e2ds'] = input_data[obs]['e2ds']
            processed[obs]['e2ds_err'] = input_data[obs]['e2ds_err']

            processed[obs]['flux'] = input_data[obs]['e2ds']/calib_data['blaze']/input_data[obs]['step']
            processed[obs]['flux_err'] = np.sqrt(input_data[obs]['e2ds'])/calib_data['blaze']/input_data[obs]['step']

            processed[obs]['flux_SRF'] = \
                rebin_2d_to_1d(input_data[obs]['wave'],
                               input_data[obs]['step'],
                               input_data[obs]['e2ds'],
                               calib_data['blaze'],
                               input_data['coadd']['wave'],
                               input_data['coadd']['step'],
                               rv_shift=observational_pams[obs]['rv_shift_ORF2SRF_mod'])

            processed[obs]['flux_SRF_err'] = \
                rebin_2d_to_1d(input_data[obs]['wave'],
                               input_data[obs]['step'],
                               input_data[obs]['e2ds_err'],
                               calib_data['blaze'],
                               input_data['coadd']['wave'],
                               input_data['coadd']['step'],
                               rv_shift=observational_pams[obs]['rv_shift_ORF2SRF_mod'],
                               is_error=True)

            """ Zero or negative values are identified, flagged and substituted with another value """
            processed[obs]['flux_SRF'], processed[obs]['flux_SRF_err'], processed[obs]['null'] = \
                replace_values_errors(processed[obs]['flux_SRF'], processed[obs]['flux_SRF_err'], 0.0001)

            """rescaling"""
            processed[obs]['flux_SRF_rescaling'], processed[obs]['flux_SRF_rescaled'], processed[obs]['flux_SRF_rescaled_err']  = \
                perform_rescaling(input_data['coadd']['wave'],
                                  processed[obs]['flux_SRF'],
                                  processed[obs]['flux_SRF_err'],
                                  observational_pams['wavelength_rescaling'])

            processed[obs]['logI'] = np.log(processed[obs]['flux_SRF_rescaled'])
            processed[obs]['logI_err'] = processed[obs]['flux_SRF_rescaled_err'] / processed[obs]['flux_SRF_rescaled']

        processed['telluric'] = {}

        n_coadd = np.size(input_data['coadd']['wave'])

        abs_slope = np.ones(n_coadd, dtype=np.double)
        line_shift = np.ones(n_coadd, dtype=np.double)
        zero_point = np.ones(n_coadd, dtype=np.double)
        pearson_r = np.zeros(n_coadd, dtype=np.double)
        pearson_p = np.zeros(n_coadd, dtype=np.double)

        airmass = np.zeros(lists['n_tellurics'], dtype=np.double)
        berv = np.zeros(lists['n_tellurics'], dtype=np.double)
        rvc = np.zeros(lists['n_tellurics'], dtype=np.double)

        for n_obs, obs in enumerate(lists['telluric']):
            # This is to ensure that airmass, berv and rvc are associated to the correct spectra
            processed['telluric'][obs] = {'n_obs': n_obs}
            airmass[n_obs] = observational_pams[obs]['AIRMASS']
            berv[n_obs] = observational_pams[obs]['BERV']
            rvc[n_obs] = observational_pams[obs]['RVC']

        logi_array = np.empty([lists['n_tellurics'], n_coadd], dtype=np.double)
        sigi_array = np.empty([lists['n_tellurics'], n_coadd], dtype=np.double)

        for obs in lists['telluric']:

            n_obs = processed['telluric'][obs]['n_obs']
            logi_array[n_obs, :] = processed[obs]['logI'][:]
            sigi_array[n_obs, :] = processed[obs]['logI_err'][:]

            """ The user has the option to select between different approaches to
                extract the telluric absorption spectrum
                To-Do: move this section to a subroutine for cythonization"""

        if observational_pams['linear_fit_method'] == 'linear_curve_fit':
            abs_slope, zero_point = \
                airmass_linear_curve_fit(airmass, logi_array, sigi_array, n_coadd)

        else:

            abs_slope, zero_point = \
                airmass_linear_lstsq(airmass, logi_array)

        telluric['stellarRF'] = {
            'wave': input_data['coadd']['wave'],
            'step': input_data['coadd']['step']
        }

        telluric['stellarRF']['spectrum'] = np.exp(abs_slope)

        telluric['stellarRF']['emission'] = (telluric['stellarRF']['spectrum'] > 1.00000)
        telluric['stellarRF']['spectrum_fixed'] = telluric['stellarRF']['spectrum'][:]
        telluric['stellarRF']['spectrum_fixed'][telluric['stellarRF']['emission']]= 1.000

        telluric['stellarRF']['spline_eval'], \
        telluric['stellarRF']['spline_coeff'], \
        telluric['stellarRF']['spline_knots'] = \
        compute_spline(input_data['coadd']['wave'],
                       telluric['stellarRF']['spectrum_fixed'],
                       0.05)

        telluric['airmass_ref'] = processed['airmass_ref']

        """ Moving the spline to the observerRF in the e2ds"""
        for obs in lists['observations']:
            """ 1) shifting the telluric correction spline to the observer RV"""

            wave_ORF = shift_wavelength_array(np.asarray(telluric['stellarRF']['spline_coeff'][0]),
                                              -observational_pams[obs]['rv_shift_ORF2SRF_mod'])

            """ 2) new spline coefficients """
            tck1 = [shift_wavelength_array(np.asarray(telluric['stellarRF']['spline_coeff'][0]),
                                              - observational_pams[obs]['rv_shift_ORF2SRF_mod']),
                    telluric['stellarRF']['spline_coeff'][1],
                    telluric['stellarRF']['spline_coeff'][2]]

            """ 3) computation of the spline at the location of the spectra, taking care of the regions 
            out of the coadded spectrum """

            inside_spline = (input_data[obs]['wave'] > wave_ORF[0]) & (input_data[obs]['wave'] < wave_ORF[-1])

            telluric[obs]['spline_noairmass'] = np.ones([input_data[obs]['n_orders'],
                                                         input_data[obs]['n_pixels']],
                                                        dtype=np.double)

            for order in range(0, input_data[obs]['n_orders']):

                if np.sum(inside_spline[order, :])>0 :
                    telluric[obs]['spline_noairmass'][order, inside_spline[order, :]] = \
                        sci_int.splev(input_data[obs]['wave'][order, inside_spline[order, :]], tck1)

            telluric[obs]['spline'] = np.power(telluric[obs]['spline_noairmass'],
                                               observational_pams[obs]['AIRMASS'] - processed['airmass_ref'])

            """ Now a similar approach is followed for the telluric spectrum before spline fit 
            """

            telluric[obs]['spectrum_noairmass'] = \
                rebin_1d_to_2d(input_data['coadd']['wave'],
                               input_data['coadd']['step'],
                               telluric['stellarRF']['spectrum'],
                               input_data[obs]['wave'],
                               input_data[obs]['step'],
                               rv_shift=-observational_pams[obs]['rv_shift_ORF2SRF_mod'],
                               preserve_flux=False)

            telluric[obs]['null'] = telluric[obs]['spectrum_noairmass'] < 0.001
            telluric[obs]['spectrum_noairmass'][ telluric[obs]['null']] = 1.0

            telluric[obs]['spectrum'] = np.power(telluric[obs]['spectrum_noairmass'],
                                               input_data[obs]['AIRMASS'] - processed['airmass_ref'])

            telluric[obs]['airmass'] = input_data[obs]['AIRMASS']
            telluric[obs]['airmass_ref'] = processed['airmass_ref']
            telluric[obs]['rv_shift_ORF2SRF_mod'] = observational_pams[obs]['rv_shift_ORF2SRF_mod']

        save_to_cpickle('telluric_processed', processed, config_in['output'], night)
        save_to_cpickle('telluric', telluric, config_in['output'], night)

        print()
        print("Night ", night, " completed")


def plot_telluric_airmass_stellarRF(config_in, night_input=''):
    import matplotlib.pyplot as plt

    night_dict = from_config_get_nights(config_in)
    instrument_dict = from_config_get_instrument(config_in)
    system_dict = from_config_get_system(config_in)

    if night_input == '':
        night_list = night_dict
    else:
        night_list = np.atleast_1d(night_input)

    for night in night_list:

        print("plot_telluric_airmass_stellarRF            Night: ", night)

        """ Retrieving the list of observations"""
        lists = load_from_cpickle('lists', config_in['output'], night)
        observational_pams = load_from_cpickle('observational_pams', config_in['output'], night)

        """ Retrieving the analysis"""
        try:
            processed = load_from_cpickle('telluric_processed', config_in['output'], night)
            telluric = load_from_cpickle('telluric', config_in['output'], night)
        except:
            print()
            print("No telluric correction, no plots")
            continue

        colors, cmap, line_colors = make_color_array(lists, observational_pams)

        fig = plt.figure(figsize=(12, 6))
        gs = GridSpec(2, 2, width_ratios=[50, 1])

        ax1 = plt.subplot(gs[0, 0])
        ax2 = plt.subplot(gs[1, 0], sharex=ax1)
        cbax1 = plt.subplot(gs[:, 1])

        lift_spectrum = 0.25

        for i, obs in enumerate(lists['observations']):
            color_array = cmap(i / len(lists['observations']))

            _, e2ds_rescaled , _ = \
                perform_rescaling(processed[obs]['wave'],
                                  processed[obs]['e2ds'],
                                  processed[obs]['e2ds_err'],
                                  observational_pams['wavelength_rescaling'])

            e2ds_rescaled_corrected_spectrum = e2ds_rescaled / telluric[obs]['spectrum']
            e2ds_rescaled_corrected_spline = e2ds_rescaled / telluric[obs]['spline']

            for order in range(0, processed[obs]['n_orders']):

                if order == 0 and i==0:
                    ax1.plot(processed[obs]['wave'][order, :],
                             e2ds_rescaled[order, :],
                             c=color_array, lw=1, alpha=0.5, label='uncorrected')
                    ax1.scatter(processed[obs]['wave'][order, :],
                                e2ds_rescaled_corrected_spectrum[order, :],
                                s=1, c=np.atleast_2d(color_array), label='corrected')
                else:
                    ax1.plot(processed[obs]['wave'][order, :],
                                e2ds_rescaled[order, :],
                                c=color_array, lw=1, alpha=0.5)
                    ax1.scatter(processed[obs]['wave'][order, :],
                                e2ds_rescaled_corrected_spectrum[order, :],
                                s=1, c=np.atleast_2d(color_array))

                #ax1.plot(processed[obs]['wave'][order, :],
                #            e2ds_rescaled[order, :]+lift_spectrum,
                #            c=color_array, lw=1, alpha=0.5)
                #ax1.scatter(processed[obs]['wave'][order, :],
                #            e2ds_rescaled_corrected_spline[order, :]+lift_spectrum,
                #            s=1, c=np.atleast_2d(color_array))

                ax2.plot(processed[obs]['wave'][order, :],
                         telluric[obs]['spectrum'][order, :],
                         c=color_array)
                ax2.axhline(1.00, c='k')

                #ax2.plot(processed[obs]['wave'][order, :],
                #         telluric[obs]['spline'][order, :]+lift_spectrum,
                #         c=color_array)
                #ax2.axhline(1.00+lift_spectrum, c='k')

        #ax2.plot(input_data['coadd']['wave'],telluric['stellarRF']['spline_eval']+0.1,c='k')
        #ax2.scatter(input_data['coadd']['wave'],telluric['stellarRF']['spectrum']+0.1,c='r', s=2)

        ax1.legend(loc=3)
        ax1.set_title('Night: ' + night)

        ax2.set_xlabel('$\lambda$ [$\AA$]')

        sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=colors[0], vmax=colors[-1]))
        sm.set_array([])  # You have to set a dummy-array for this to work...
        cbar = plt.colorbar(sm, cax=cbax1)
        cbar.set_label('BJD - 2450000.0')
        fig.subplots_adjust(wspace=0.05, hspace=0.4)
        plt.show()

