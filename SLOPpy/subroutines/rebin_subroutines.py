from __future__ import print_function, division
import numpy as np
from scipy.interpolate import interp1d
from SLOPpy.subroutines.constants import *


def shift_wavelength(wave, step, rv_shift):
    wave_shift = rv_shift / speed_of_light_km + 1.00000
    return wave * wave_shift, step * wave_shift


def shift_wavelength_array(wave, rv_shift):
    wave_shift = rv_shift / speed_of_light_km + 1.00000
    return wave * wave_shift


def shift_wavelength_to_rest(wave, step, rv_shift):
    inverse_wave_shift = (-rv_shift) / speed_of_light_km + 1.00000
    return wave / inverse_wave_shift, step / inverse_wave_shift


def rebin_exact_flux(wave_in, step_in, flux_in, wave_out, step_out,
                     quadrature=False,
                     preserve_flux=True):
    """
    Previously named rebin_order
    :param wave_in:
    :param step_in:
    :param flux_in:
    :param wave_out:
    :param step_out:
    :param quadrature:
    :param preserve_flux:
    :return:

    Spectral rebinning with flux conservation
    """
    if quadrature:
        flux_in = flux_in**2.

    flux_out = np.zeros(np.shape(wave_out), dtype=np.double)
    n1 = np.size(wave_in)
    n2 = np.size(wave_out)
    ns_prv = 0

    for i in range(0, n2):
        # print i, ' of ', n2
        # Starting and ending point of the bin
        wlb = wave_out[i] - step_out[i] / 2.000
        wle = wave_out[i] + step_out[i] / 2.000

        # Normalized flux value within the bin
        fl_nm = 0.00

        # b->blue and r->red side of the original spectrum which include the bin
        # ib and ie are initialized with values close to the ones of the last iteration to save time
        ib = ns_prv
        ir = ns_prv

        for ns in range(ns_prv, n1 - 1):
            # simple algorithm to search the closest indexes near the bin boundaries
            if wave_in[ib] + step_in[ib] / 2.00 < wlb: ib += 1
            if wave_in[ir] + step_in[ir] / 2.00 < wle: ir += 1

            # when we are close to the boundary of the spectra, we stop
            if ir < ns - 3: break

        # Fail-safe checks
        if ib > ns_prv: ns_prv = ib - 3
        if ib < 0 or ir > n1: continue
        if ib > ir: continue
        if ns_prv < 0: ns_prv = 0

        # Now the true rebinning section
        if ib == ir:
            pix_s = (wle - wlb) / step_in[ib]  # fraction
            pix_e = 0.
            flux_out[i] += pix_s * flux_in[ib]
            fl_nm += pix_s
        elif ib + 1 == ir:
            pix_s = (wave_in[ib] + step_in[ib] * 0.5 - wlb) / step_in[ib]
            pix_e = (wle - (wave_in[ir] - step_in[ir] * 0.5)) / step_in[ir]
            flux_out[i] += (pix_s * flux_in[ib] + pix_e * flux_in[ir])
            fl_nm += (pix_s + pix_e)
        else:
            pix_s = (wave_in[ib] + step_in[ib] * 0.5 - wlb) / step_in[ib]
            pix_e = (wle - (wave_in[ir] - step_in[ir] * 0.5)) / step_in[ir]
            flux_out[i] += (pix_s * flux_in[ib] + pix_e * flux_in[ir])
            fl_nm += (pix_s + pix_e)
            for j in range(ib + 1, ir):
                flux_out[i] += flux_in[j]
                fl_nm += 1.00
        if (not preserve_flux) and fl_nm > 0.0:
            if quadrature:
                fl_nm *= fl_nm
            flux_out[i] /= fl_nm

    if quadrature:
        return np.sqrt(flux_out)
    else:
        return flux_out

def rebin_with_interpolation(wave_in, step_in, flux_in, wave_out, step_out,
                             quadrature=False,
                             preserve_flux=True,
                             interp_kind='cubic'):

    ndata = len(wave_in)

    normalization_factor = 1.0
    if preserve_flux:
        step_in_internal = np.ones(ndata)
        step_out_internal = np.ones(len(step_out))
    else:
        step_in_internal = step_in
        step_out_internal = step_out
        if quadrature:
            normalization_factor = (np.median(step_out) / np.median(step_in))
    if quadrature:
        flux_in = np.power(flux_in, 2.)

    wave_in_cumul = np.zeros(ndata+1)
    flux_in_cumul = np.zeros(ndata+1)

    flux_in_cumul[0] = 0.0

    wave_in_cumul[0] = wave_in[0] - step_in[0] / 2.0

    for i in range(1, ndata):
        flux_in_cumul[i] = flux_in_cumul[i - 1] + flux_in[i - 1] * step_in_internal[i - 1]
        # wave_in_cumul[i] = wave_in[i]-step_in[i]/2.
        wave_in_cumul[i] = wave_in[i] - (wave_in[i] - wave_in[i - 1]) / 2.

    flux_in_cumul[ndata] = flux_in_cumul[ndata - 1] + flux_in[ndata - 1] * step_in_internal[ndata - 1]
    wave_in_cumul[ndata] = wave_in[ndata - 1] + step_in[ndata - 1] / 2.

    flux_cumul_interp1d = interp1d(wave_in_cumul, flux_in_cumul, kind=interp_kind, bounds_error=False, fill_value=0.000)

    flux_out = (flux_cumul_interp1d(wave_out + step_out / 2.) - flux_cumul_interp1d(
        wave_out - step_out / 2.)) / step_out_internal

    if quadrature:
        return np.sqrt(flux_out) / normalization_factor
    else:
        return flux_out


def rebin_1d_to_1d(wave_in, step_in, flux_in, wave_out, step_out,
                   rv_shift=None,
                   is_error=False,
                   quadrature=False,
                   preserve_flux=True,
                   method='cubic_interpolation',
                   reference_value=None):

    if is_error:
        quadrature = True
        method='exact_flux'

    if rv_shift:
        wave_in, step_in = shift_wavelength(wave_in, step_in, rv_shift)

    if method == 'exact_flux':
        flux_out = rebin_exact_flux(wave_in, step_in, flux_in, wave_out, step_out,
                                    quadrature=quadrature, preserve_flux=preserve_flux)
    elif method == 'cubic_interpolation':
        flux_out = rebin_with_interpolation(wave_in, step_in, flux_in, wave_out, step_out,
                                            quadrature=quadrature, preserve_flux=preserve_flux, interp_kind='cubic')


    else:
        raise ValueError("method ", method, 'not supported by rebinning subroutine')


    if reference_value :
        wave_sel = (wave_out<=wave_in[0] +0.005 ) | (wave_out>=wave_in[-1] -0.005)
        flux_out[wave_sel] = reference_value

    return flux_out


def rebin_2d_to_1d(wave_in, step_in, flux_in, blaze_in, wave_out, step_out,
                   rv_shift=None,
                   is_error=False,
                   quadrature=False,
                   preserve_flux=True,
                   skip_blaze_correction=False,
                   method='cubic_interpolation',
                   reference_value=None):

    """
    :param wave_in:
    :param step_in:
    :param flux_in:
    :param blaze_in:
    :param wave_out:
    :param step_out:
    :param rv_shift:
    :param is_error:
    :param quadrature:
    :param preserve_flux:
    :param skip_blaze_correction:
    :param method:
    :return: flux_out
    """

    if rv_shift:
        wave_in, step_in = shift_wavelength(wave_in, step_in, rv_shift)

    o_axis, f_axis = np.shape(wave_in)
    n_rebin = np.size(wave_out)

    if skip_blaze_correction or blaze_in is None:
        flux_deblazed_in = flux_in
    else:
        flux_deblazed_in = flux_in / blaze_in

    if is_error:
        quadrature = True
        method = 'exact_flux'

    # Rebinning of the individual orders. We keep track of starting
    # and ending points of each order in the rebinned solution
    flux_rebin_pix = np.zeros([o_axis, n_rebin], dtype=np.double)

    counter_is = np.zeros(o_axis, dtype=np.int64) - 1
    counter_ie = np.zeros(o_axis, dtype=np.int64) - 1
    skip_order = np.ones(o_axis, dtype=bool)
    for ii in range(0, o_axis):

        counter_is[ii] = np.argmin(abs(wave_in[ii, 0] - wave_out))
        counter_ie[ii] = np.argmin(abs(wave_in[ii, -1] - wave_out))

        if wave_in[ii, 0] > np.amax(wave_out) or wave_in[ii, -1] < np.amin(wave_out):
            continue

        skip_order[ii] = False

        i = counter_is[ii]
        j = counter_ie[ii]+1

        flux_rebin_pix[ii, i:j] = rebin_1d_to_1d(wave_in[ii, :],
                                                 step_in[ii, :],
                                                 flux_deblazed_in[ii, :],
                                                 wave_out[i:j],
                                                 step_out[i:j],
                                                 quadrature=quadrature,
                                                 is_error=is_error,
                                                 preserve_flux=preserve_flux,
                                                 method=method,
                                                 reference_value=reference_value)

    flux_out = np.zeros(n_rebin, dtype=np.double)
    if reference_value:
        flux_out += reference_value

    if quadrature or is_error:
        flux_rebin_pix = np.power(flux_rebin_pix, 2)
        flux_out[counter_is[0]:counter_ie[0]] = flux_rebin_pix[0, counter_is[0]:counter_ie[0]]

        for ii in range(1, o_axis):
            if skip_order[ii]: continue

            p_ie = counter_ie[ii - 1]
            j_is = counter_is[ii]
            j_ie = counter_ie[ii] + 1 # adding one because it is used in interval definition - Python quirks

            if p_ie > j_is:
                nr_joint = float(p_ie - j_is)
                ij = np.arange(j_is, p_ie, 1, dtype=np.int64)
                ij_fraction = np.power((ij-j_is) / nr_joint, 4)
                flux_out[ij] = flux_rebin_pix[ii,ij] * ij_fraction + flux_rebin_pix[ii-1,ij] * (1. - ij_fraction)
                flux_out[p_ie:j_ie] = flux_rebin_pix[ii, p_ie:j_ie]
            else:
                flux_out[j_is:j_ie] = flux_rebin_pix[ii, j_is:j_ie]

        return np.sqrt(flux_out)

    else:
        flux_out[counter_is[0]:counter_ie[0]] = flux_rebin_pix[0, counter_is[0]:counter_ie[0]]
        for ii in range(1, o_axis):
            if skip_order[ii]: continue

            p_ie = counter_ie[ii - 1]
            j_is = counter_is[ii]
            j_ie = counter_ie[ii] + 1

            if p_ie > j_is:
                nr_joint = float(p_ie - j_is)
                ij = np.arange(j_is, p_ie, 1, dtype=np.int64)
                ij_fraction = np.power((ij-j_is) / nr_joint, 2.)
                flux_out[ij] = flux_rebin_pix[ii,ij] * ij_fraction + flux_rebin_pix[ii-1,ij] * (1. - ij_fraction)
                flux_out[p_ie:j_ie] = flux_rebin_pix[ii, p_ie:j_ie]
            else:
                flux_out[j_is:j_ie] = flux_rebin_pix[ii, j_is:j_ie]
        return flux_out


def rebin_1d_to_2d(wave_in, step_in, flux_in, wave_out, step_out,
                   rv_shift=None,
                   is_error=False,
                   quadrature=False,
                   preserve_flux=True,
                   method='cubic_interpolation',
                   reference_value=None):
    """

    :param wave_in:
    :param step_in:
    :param flux_in:
    :param wave_out:
    :param step_out:
    :param rv_shift:
    :param is_error:
    :param quadrature:
    :param preserve_flux:
    :param method:
    :return:
    """

    if rv_shift:
        wave_in, step_in = shift_wavelength(wave_in, step_in, rv_shift)

    if is_error:
        quadrature = True
        method='exact_flux'

    o_axis, f_axis = np.shape(wave_out)
    flux_out = np.zeros([o_axis, f_axis], dtype=np.double)
    if reference_value:
        flux_out += reference_value

    for ii in range(0, o_axis):
        flux_out[ii, :]= rebin_1d_to_1d(wave_in,
                                        step_in,
                                        flux_in,
                                        wave_out[ii, :],
                                        step_out[ii, :],
                                        quadrature=quadrature,
                                        is_error=is_error,
                                        preserve_flux=preserve_flux,
                                        method=method,
                                        reference_value=reference_value)


    return flux_out


def rebin_2d_to_2d(wave_in, step_in, flux_in, wave_out, step_out,
                   rv_shift=None,
                   is_error=False,
                   quadrature=False,
                   preserve_flux=True,
                   method='cubic_interpolation',
                   reference_value=None):

    """
    :param wave_in: 1D wavelength array of the input spectrum
    :param step_in: 1D step-size array of the input spectrum
    :param flux_in: 1D flux array of the input spectrum
    :param wave_out: 2D (order-by-order) wavelength array of the rebinned spectrum
    :param step_out: 2D (order-by-order) step-size array of the rebinned spectrum
    :param rv_shift:
    :param is_error:
    :param quadrature:
    :param preserve_flux:
    :param method:
    :return: flux_out, 2D (order-by-order) flux array, with same size as wave_out
    """

    if rv_shift:
        wave_in, step_in = shift_wavelength(wave_in, step_in, rv_shift)

    if is_error:
        quadrature = True
        method='exact_flux'

    o_axis_input, f_axis_input = np.shape(wave_in)
    o_axis, f_axis = np.shape(wave_out)

    if o_axis_input != o_axis:
        raise ValueError("Mismatch between input and output number of orders in rebin_2d_to_2d")

    flux_out = np.zeros([o_axis, f_axis], dtype=np.double)
    if reference_value:
        flux_out += reference_value

    for ii in range(0, o_axis):
        flux_out[ii, :] = rebin_1d_to_1d(wave_in[ii, :],
                                         step_in[ii, :],
                                         flux_in[ii, :],
                                         wave_out[ii, :],
                                         step_out[ii, :],
                                         quadrature=quadrature,
                                         is_error=is_error,
                                         preserve_flux=preserve_flux,
                                         method=method,
                                         reference_value=reference_value)

    return flux_out