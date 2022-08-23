from __future__ import print_function, division
from SLOPpy.subroutines.common import *
from SLOPpy.subroutines.io_subroutines import get_filename


def get_eso_sckycalc_harps(obs_ref, wave_range, ra, dec, night, output):

    # https://www.eso.org/observing/etc/doc/skycalc/helpskycalccli.html

    time_tag = obs_ref[6:25]
    wave_range_nm = wave_range/10.
    wdelta = 0.00075

    input_filename = get_filename('skycalc_input_' + repr(wave_range_nm[0]) + '_' + repr(wave_range_nm[1]),
                                  output, night, extension=".JSON")

    alman_filename = get_filename('skycalc_alman_' + time_tag + '_' + repr(wave_range_nm[1]),
                                  output, night, extension=".JSON")

    output_filename = get_filename('skycalc_output_' + repr(wave_range_nm[0]) + '_' + repr(wave_range_nm[1]) + time_tag,
                                  output, night, extension=".fits")

    if os.path.isfile(output_filename):
        skycalc_hdu = fits.open(night + '.fits')
        data = skycalc_hdu[1].data

        skycalc_hdu.close()
        return data.field(0) * 10., \
               np.ones(len(data.field(0))) * wdelta, \
               data.field(14), \
               np.ones(len(data.field(0)))

    if not os.path.isfile(input_filename):
        input_pams = {
            "pwv_mode": "pwv",
            "incl_moon": "N", # No moon contamination
            "incl_starlight": "N", # No starlight
            "incl_zodiacal": "N", # No zodiacal light
            "incl_loweratm": "Y",
            "incl_upperatm": "Y",
            "incl_airglow": "N", # No airglow
            "incl_therm": "N",
            "vacair": "air", # compute in the air
            "wmin": wave_range[0],
            "wmax": wave_range[1],
            "wgrid_mode": "fixed_wavelength_step",
            "wdelta": wdelta,
            "wres": 20000,
            "lsf_type": "Gaussian", # Gaussian lsf
            "lsf_gauss_fwhm": 5.5,
        }

        with open(input_filename, 'w') as outfile:
            json.dump(input_pams, outfile)

    if not os.path.isfile(alman_filename):
        almanac_pams = {
            "ra": ra,
            "dec": dec,
            "date": time_tag,
            "observatory": "lasilla"
        }

        with open(alman_filename, 'w') as outfile:
            json.dump(almanac_pams, outfile)

    os.system('skycalc_cli'
               + ' --in ' + input_filename
               + ' --alm ' + alman_filename
               + ' --out ' + output_filename)

    skycalc_hdu = fits.open(output_filename)
    data = skycalc_hdu[1].data

    skycalc_hdu.close()

    return data.field(0) * 10., \
           np.ones(len(data.field(0))) * wdelta, \
           data.field(14), \
           np.ones(len(data.field(0)))
