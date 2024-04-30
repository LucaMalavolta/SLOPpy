import numpy as np
from SLOPpy.subroutines.math_functions import interpolate2d_grid_nocheck
from SLOPpy.subroutines.constants import *

def compute_single_line(wave_meshgrid, planet_RVsinusoid, planet_K, line_array):

    """ computing the spectral shift in RV """
    rv_shift = planet_K * planet_RVsinusoid

    line_model = np.ones(np.shape(wave_meshgrid), dtype= np.double)

    """ line_array is always structured in this way:
        line_array = np.empty(n_pams, n_lines)
        line_array[0, 0] = wavelength
        line_array[1, 0] = contrast
        line_array[2, 0] = FWHM
        line_array[3, 0] = winds
    """

    sigma =  line_array[2] / sigma2fwhm * line_array[0] / speed_of_light_km
    line_shifted = line_array[0] +  (rv_shift + line_array[3])  * line_array[0] / speed_of_light_km
    line_model -= line_array[1]  * np.exp(-(1./(2*sigma**2))*(wave_meshgrid-line_shifted)**2)

    return line_model


def compute_multiple_lines(wave_meshgrid, planet_RVsinusoid, planet_K, lines_array, n_lines):

    """ computing the spectral shift in RV """
    rv_shift = planet_K * planet_RVsinusoid

    lines_model = np.ones(np.shape(wave_meshgrid), dtype= np.double)

    for ii in range(0, n_lines):
        sigma =  lines_array[2,ii] / sigma2fwhm * lines_array[0,ii] / speed_of_light_km
        line_shifted = lines_array[0,ii] +  (rv_shift + lines_array[3,ii])  * lines_array[0,ii] / speed_of_light_km
        lines_model -= lines_array[1,ii]  * np.exp(-(1./(2*sigma**2))*(wave_meshgrid-line_shifted)**2)

    return lines_model


""" case  0: only one spectral line, default line parameters are contrast, FWHM, rv_shift """
def logprob_case00(theta,
            boundaries,
            wave_meshgrid,
            transmission_spec,
            transmission_spec_err,
            clv_rm_radius,
            clv_rm_grid,
            planet_RVsinusoid,
            lines_center,
            jitter_index,
            priors_dict,
            return_models=False
            ):

    """ check the boundaries """
    if ((theta < boundaries[:,0]) | (theta > boundaries[:,1])).any():
        return -np.inf

    """ unfolding jitter parameters """
    if jitter_index is None:
        i_j = -1
        jitter_array = 0.
        jitter_pams = 0.
    elif len(jitter_index) > 0:
        jitter_array = jitter_index * 0.
        n_jitter = int(np.amax(jitter_index) + 1)
        jitter_pams = np.empty(n_jitter)
        for i_j in range(0, n_jitter):
            sel = (jitter_index == i_j)
            jitter_array[sel] = theta[-1-i_j]
            jitter_pams[i_j] = theta[-1-i_j]
    else:
        i_j = 0
        jitter_array = wave_meshgrid*0. + theta[-1]
        jitter_pams = theta[-1]

    """ the second-last value after jitter is always the semi-amplitude of the planet
        the third-last value after jitter is always the planetary radius, if included in the model
    """
    planet_K = theta[-2-i_j]
    planet_R = theta[-3-i_j]

    """ computing interpolated model spectrum
    """
    clv_model = interpolate2d_grid_nocheck(planet_R, clv_rm_radius, clv_rm_grid)

    """ line_array is always structured in this way:
        line_array = np.empty(n_pams, n_lines)
        line_array[0, 0] = wavelength
        line_array[1, 0] = contrast
        line_array[2, 0] = FWHM
        line_array[3, 0] = winds
    """

    line_array = [lines_center[0], theta[0], theta[1], theta[2]]
    line_model = compute_single_line(wave_meshgrid, planet_RVsinusoid, planet_K, line_array)

    flux_res = transmission_spec / clv_model / line_model - 1.
    ferr_res = transmission_spec_err / clv_model / line_model

    if return_models:
        lines_array = np.empty([4, 1])
        lines_array[:, 0] = line_array
        return line_model, clv_model, lines_array, planet_K, planet_R, jitter_pams

    var_dict = {
        'planet_K': planet_K,
        'planet_R': planet_R,
        'FWHM': line_array[2],
        'winds': line_array[3],
    }
    log_prior = 0.
    for key_name, key_vals in priors_dict.items():
        log_prior +=  (-(var_dict[key_name] - key_vals[0]) ** 2 / (2 * key_vals[1] ** 2) - 0.5 * np.log(2*np.pi) - np.log(key_vals[1]))

    """ adding the jitter to error esitamets """
    env = 1.0 / (jitter_array ** 2.0 + ferr_res ** 2.0)
    return log_prior -0.5 * (np.size(flux_res) * np.log(2 * np.pi) +
                    np.sum((flux_res) ** 2 * env - np.log(env)))


""" case  1: only one spectral line, no winds """
def logprob_case01(theta,
            boundaries,
            wave_meshgrid,
            transmission_spec,
            transmission_spec_err,
            clv_rm_radius,
            clv_rm_grid,
            planet_RVsinusoid,
            lines_center,
            jitter_index,
            priors_dict,
            return_models=False
            ):

    """ check the boundaries """
    if ((theta < boundaries[:,0]) | (theta > boundaries[:,1])).any():
        return -np.inf

    """ unfolding jitter parameters """
    if jitter_index is None:
        i_j = -1
        jitter_array = 0.
        jitter_pams = 0.
    elif len(jitter_index) > 0:
        jitter_array = jitter_index * 0.
        n_jitter = int(np.amax(jitter_index) + 1)
        jitter_pams = np.empty(n_jitter)
        for i_j in range(0, n_jitter):
            sel = (jitter_index == i_j)
            jitter_array[sel] = theta[-1-i_j]
            jitter_pams[i_j] = theta[-1-i_j]
    else:
        i_j = 0
        jitter_array = wave_meshgrid*0. + theta[-1]
        jitter_pams = theta[-1]

    """ the second-last value after jitter is always the semi-amplitude of the planet
        the third-last value after jitter is always the planetary radius, if included in the model
    """
    planet_K = theta[-2-i_j]
    planet_R = theta[-3-i_j]

    """ computing interpolated model spectrum
    """
    clv_model = interpolate2d_grid_nocheck(planet_R, clv_rm_radius, clv_rm_grid)

    """ line_array is always structured in this way:
        line_array = np.empty(n_pams, n_lines)
        line_array[0, 0] = wavelength
        line_array[1, 0] = contrast
        line_array[2, 0] = FWHM
        line_array[3, 0] = winds
    """

    line_array = [lines_center[0], theta[0], theta[1], 0.000]
    line_model = compute_single_line(wave_meshgrid, planet_RVsinusoid, planet_K, line_array)

    flux_res = transmission_spec / clv_model / line_model - 1.
    ferr_res = transmission_spec_err / clv_model / line_model

    if return_models:
        lines_array = np.empty([4, 1])
        lines_array[:,0] = line_array
        return line_model, clv_model, lines_array, planet_K, planet_R, jitter_pams

    var_dict = {
        'planet_K': planet_K,
        'planet_R': planet_R,
        'FWHM': line_array[2],
        'winds': line_array[3],
    }
    log_prior = 0.
    for key_name, key_vals in priors_dict.items():
        log_prior +=  (-(var_dict[key_name] - key_vals[0]) ** 2 / (2 * key_vals[1] ** 2) - 0.5 * np.log(2*np.pi) - np.log(key_vals[1]))

    """ the last value is always the jitter """
    env = 1.0 / (jitter_array ** 2.0 + ferr_res ** 2.0)
    return log_prior -0.5 * (np.size(flux_res) * np.log(2 * np.pi) +
                    np.sum((flux_res) ** 2 * env - np.log(env)))


""" case  2: only one spectral line, no planetary radius dependance """
def logprob_case02(theta,
            boundaries,
            wave_meshgrid,
            transmission_spec,
            transmission_spec_err,
            clv_model,
            planet_RVsinusoid,
            lines_center,
            jitter_index,
            priors_dict,
            return_models=False
            ):

    """ check the boundaries """
    if ((theta < boundaries[:,0]) | (theta > boundaries[:,1])).any():
        return -np.inf

    """ unfolding jitter parameters """
    if jitter_index is None:
        i_j = -1
        jitter_array = 0.
        jitter_pams = 0.
    elif len(jitter_index) > 0:
        jitter_array = jitter_index * 0.
        n_jitter = int(np.amax(jitter_index) + 1)
        jitter_pams = np.empty(n_jitter)
        for i_j in range(0, n_jitter):
            sel = (jitter_index == i_j)
            jitter_array[sel] = theta[-1-i_j]
            jitter_pams[i_j] = theta[-1-i_j]
    else:
        i_j = 0
        jitter_array = wave_meshgrid*0. + theta[-1]
        jitter_pams = theta[-1]

    """ the second-last value after jitter is always the semi-amplitude of the planet
    """
    planet_K = theta[-2-i_j]

    line_array = [lines_center[0], theta[0], theta[1], theta[2]]
    line_model = compute_single_line(wave_meshgrid, planet_RVsinusoid, planet_K, line_array)

    flux_res = transmission_spec / clv_model / line_model - 1.
    ferr_res = transmission_spec_err / clv_model / line_model

    if return_models:
        lines_array = np.empty([4, 1])
        lines_array[:,0] = line_array
        return line_model, clv_model, lines_array, planet_K, 1.00000, jitter_pams

    var_dict = {
        'planet_K': planet_K,
        'planet_R': 1.000000,
        'FWHM': line_array[2],
        'winds': line_array[3],
    }
    log_prior = 0.
    for key_name, key_vals in priors_dict.items():
        log_prior +=  (-(var_dict[key_name] - key_vals[0]) ** 2 / (2 * key_vals[1] ** 2) - 0.5 * np.log(2*np.pi) - np.log(key_vals[1]))

    """ the last value is always the jitter """
    env = 1.0 / (jitter_array ** 2.0 + ferr_res ** 2.0)
    return log_prior -0.5 * (np.size(flux_res) * np.log(2 * np.pi) +
                    np.sum((flux_res) ** 2 * env - np.log(env)))


""" case  3: only one spectral line, no winds and no planetary radius dependance """
def logprob_case03(theta,
            boundaries,
            wave_meshgrid,
            transmission_spec,
            transmission_spec_err,
            clv_model,
            planet_RVsinusoid,
            lines_center,
            jitter_index,
            priors_dict,
            return_models=False
            ):

    """ check the boundaries """
    if ((theta < boundaries[:,0]) | (theta > boundaries[:,1])).any():
        return -np.inf

    """ unfolding jitter parameters """
    if jitter_index is None:
        i_j = -1
        jitter_array = 0.
        jitter_pams = 0.
    elif len(jitter_index) > 0:
        jitter_array = jitter_index * 0.
        n_jitter = int(np.amax(jitter_index) + 1)
        jitter_pams = np.empty(n_jitter)
        for i_j in range(0, n_jitter):
            sel = (jitter_index == i_j)
            jitter_array[sel] = theta[-1-i_j]
            jitter_pams[i_j] = theta[-1-i_j]
    else:
        i_j = 0
        jitter_array = wave_meshgrid*0. + theta[-1]
        jitter_pams = theta[-1]

    """ the second-last value after jitter is always the semi-amplitude of the planet
    """
    planet_K = theta[-2-i_j]

    line_array = [lines_center[0], theta[0], theta[1], 0.000]
    line_model = compute_single_line(wave_meshgrid, planet_RVsinusoid, planet_K, line_array)

    flux_res = transmission_spec / clv_model / line_model - 1.
    ferr_res = transmission_spec_err / clv_model / line_model

    if return_models:
        lines_array = np.empty([4, 1])
        lines_array[:,0] = line_array
        return line_model, clv_model, lines_array, planet_K, 1.00000, jitter_pams

    var_dict = {
        'planet_K': planet_K,
        'planet_R': 1.000000,
        'FWHM': line_array[2],
        'winds': line_array[3],
    }
    log_prior = 0.
    for key_name, key_vals in priors_dict.items():
        log_prior +=  (-(var_dict[key_name] - key_vals[0]) ** 2 / (2 * key_vals[1] ** 2) - 0.5 * np.log(2*np.pi) - np.log(key_vals[1]))

    """ the last value is always the jitter """
    env = 1.0 / (jitter_array ** 2.0 + ferr_res ** 2.0)
    return log_prior -0.5 * (np.size(flux_res) * np.log(2 * np.pi) +
                    np.sum((flux_res) ** 2 * env - np.log(env)))


""" case 10: more than one spectral lines, all line parameters are free and independent """
def logprob_case10(theta,
            boundaries,
            wave_meshgrid,
            transmission_spec,
            transmission_spec_err,
            clv_rm_radius,
            clv_rm_grid,
            planet_RVsinusoid,
            lines_center,
            jitter_index,
            priors_dict,
            return_models=False
            ):

    """ check the boundaries """
    if ((theta < boundaries[:,0]) | (theta > boundaries[:,1])).any():
        return -np.inf

    """ unfolding jitter parameters """
    if jitter_index is None:
        i_j = -1
        jitter_array = 0.
        jitter_pams = 0.
    elif len(jitter_index) > 0:
        jitter_array = jitter_index * 0.
        n_jitter = int(np.amax(jitter_index) + 1)
        jitter_pams = np.empty(n_jitter)
        for i_j in range(0, n_jitter):
            sel = (jitter_index == i_j)
            jitter_array[sel] = theta[-1-i_j]
            jitter_pams[i_j] = theta[-1-i_j]
    else:
        i_j = 0
        jitter_array = wave_meshgrid*0. + theta[-1]
        jitter_pams = theta[-1]

    """ the second-last value after jitter is always the semi-amplitude of the planet
        the third-last value after jitter is always the planetary radius, if included in the model
    """
    planet_K = theta[-2-i_j]
    planet_R = theta[-3-i_j]

    """ computing interpolated model spectrum
    """
    clv_model = interpolate2d_grid_nocheck(planet_R, clv_rm_radius, clv_rm_grid)

    """ line_array is always structured in this way:
        line_array = np.empty(n_pams, n_lines)
        line_array[0, 0] = wavelength
        line_array[1, 0] = contrast
        line_array[2, 0] = FWHM
        line_array[3, 0] = winds
    """

    n_lines = len(lines_center)
    lines_array = np.empty([4, n_lines])
    i_pams = 0
    for ii in range(0, n_lines):
        lines_array[0, ii] = lines_center[ii]
        lines_array[1, ii] = theta[i_pams]
        i_pams += 1
        lines_array[2, ii] = theta[i_pams]
        i_pams += 1
        lines_array[3, ii] = theta[i_pams]
        i_pams += 1

    lines_model = compute_multiple_lines(wave_meshgrid, planet_RVsinusoid, planet_K, lines_array, n_lines)

    flux_res = transmission_spec / clv_model / lines_model - 1.
    ferr_res = transmission_spec_err / clv_model / lines_model

    if return_models:
        return lines_model, clv_model, lines_array, planet_K, planet_R, jitter_pams

    var_dict = {
        'planet_K': planet_K,
        'planet_R': planet_R,
        'FWHM': lines_array[2,:],
        'winds': lines_array[3,:],
    }
    log_prior = 0.
    for key_name, key_vals in priors_dict.items():
        for var_value in np.atleast_1d(var_dict[key_name]):
            log_prior +=  (-(var_value - key_vals[0]) ** 2 / (2 * key_vals[1] ** 2) - 0.5 * np.log(2*np.pi) - np.log(key_vals[1]))

    """ the last value is always the jitter """
    env = 1.0 / (jitter_array ** 2.0 + ferr_res ** 2.0)
    return log_prior -0.5 * (np.size(flux_res) * np.log(2 * np.pi) +
                    np.sum((flux_res) ** 2 * env - np.log(env)))


""" case 11: more than one spectral lines, all lines are affected by the same wind """
def logprob_case11(theta,
            boundaries,
            wave_meshgrid,
            transmission_spec,
            transmission_spec_err,
            clv_rm_radius,
            clv_rm_grid,
            planet_RVsinusoid,
            lines_center,
            jitter_index,
            priors_dict,
            return_models=False
            ):

    """ check the boundaries """
    if ((theta < boundaries[:,0]) | (theta > boundaries[:,1])).any():
        return -np.inf

    """ unfolding jitter parameters """
    if jitter_index is None:
        i_j = -1
        jitter_array = 0.
        jitter_pams = 0.
    elif len(jitter_index) > 0:
        jitter_array = jitter_index * 0.
        n_jitter = int(np.amax(jitter_index) + 1)
        jitter_pams = np.empty(n_jitter)
        for i_j in range(0, n_jitter):
            sel = (jitter_index == i_j)
            jitter_array[sel] = theta[-1-i_j]
            jitter_pams[i_j] = theta[-1-i_j]
    else:
        i_j = 0
        jitter_array = wave_meshgrid*0. + theta[-1]
        jitter_pams = theta[-1]

    """ the second-last value after jitter is always the semi-amplitude of the planet
        the third-last value after jitter is always the planetary radius, if included in the model
    """
    planet_K = theta[-2-i_j]
    planet_R = theta[-3-i_j]

    """ computing interpolated model spectrum
    """
    clv_model = interpolate2d_grid_nocheck(planet_R, clv_rm_radius, clv_rm_grid)

    """ line_array is always structured in this way:
        line_array = np.empty(n_pams, n_lines)
        line_array[0, 0] = wavelength
        line_array[1, 0] = contrast
        line_array[2, 0] = FWHM
        line_array[3, 0] = winds
    """

    n_lines = len(lines_center)
    lines_array = np.empty([4, n_lines])
    i_pams = 0
    for ii in range(0, n_lines):
        lines_array[0, ii] = lines_center[ii]
        lines_array[1, ii] = theta[i_pams]
        i_pams += 1
        lines_array[2, ii] = theta[i_pams]
        i_pams += 1
        lines_array[3, ii] = theta[-4-i_j]

    lines_model = compute_multiple_lines(wave_meshgrid, planet_RVsinusoid, planet_K, lines_array, n_lines)

    flux_res = transmission_spec / clv_model / lines_model - 1.
    ferr_res = transmission_spec_err / clv_model / lines_model

    if return_models:
        return lines_model, clv_model, lines_array, planet_K, planet_R, jitter_pams

    var_dict = {
        'planet_K': planet_K,
        'planet_R': planet_R,
        'FWHM': lines_array[2,:],
        'winds': lines_array[3,:],
    }
    log_prior = 0.
    for key_name, key_vals in priors_dict.items():
        for var_value in np.atleast_1d(var_dict[key_name]):
            log_prior +=  (-(var_value - key_vals[0]) ** 2 / (2 * key_vals[1] ** 2) - 0.5 * np.log(2*np.pi) - np.log(key_vals[1]))

    """ the last value is always the jitter """
    env = 1.0 / (jitter_array ** 2.0 + ferr_res ** 2.0)
    return log_prior -0.5 * (np.size(flux_res) * np.log(2 * np.pi) +
                    np.sum((flux_res) ** 2 * env - np.log(env)))


""" case 12: more than one spectral lines, all lines have same FWHM """
def logprob_case12(theta,
            boundaries,
            wave_meshgrid,
            transmission_spec,
            transmission_spec_err,
            clv_rm_radius,
            clv_rm_grid,
            planet_RVsinusoid,
            lines_center,
            jitter_index,
            priors_dict,
            return_models=False
            ):

    """ check the boundaries """
    if ((theta < boundaries[:,0]) | (theta > boundaries[:,1])).any():
        return -np.inf

    """ unfolding jitter parameters """
    if jitter_index is None:
        i_j = -1
        jitter_array = 0.
        jitter_pams = 0.
    elif len(jitter_index) > 0:
        jitter_array = jitter_index * 0.
        n_jitter = int(np.amax(jitter_index) + 1)

        jitter_pams = np.empty(n_jitter)
        for i_j in range(0, n_jitter):
            sel = (jitter_index == i_j)

            jitter_array[sel] = theta[-1-i_j]
            jitter_pams[i_j] = theta[-1-i_j]
    else:
        i_j = 0
        jitter_array = wave_meshgrid*0. + theta[-1]
        jitter_pams = theta[-1]

    """ the second-last value after jitter is always the semi-amplitude of the planet
        the third-last value after jitter is always the planetary radius, if included in the model
    """
    planet_K = theta[-2-i_j]
    planet_R = theta[-3-i_j]

    """ computing interpolated model spectrum
    """
    clv_model = interpolate2d_grid_nocheck(planet_R, clv_rm_radius, clv_rm_grid)

    """ line_array is always structured in this way:
        line_array = np.empty(n_pams, n_lines)
        line_array[0, 0] = wavelength
        line_array[1, 0] = contrast
        line_array[2, 0] = FWHM
        line_array[3, 0] = winds
    """

    n_lines = len(lines_center)
    lines_array = np.empty([4, n_lines])
    i_pams = 0
    for ii in range(0, n_lines):
        lines_array[0, ii] = lines_center[ii]
        lines_array[1, ii] = theta[i_pams]
        i_pams += 1
        lines_array[2, ii] = theta[-4-i_j]
        lines_array[3, ii] = theta[i_pams]
        i_pams += 1

    lines_model = compute_multiple_lines(wave_meshgrid, planet_RVsinusoid, planet_K, lines_array, n_lines)

    flux_res = transmission_spec / clv_model / lines_model - 1.
    ferr_res = transmission_spec_err / clv_model / lines_model

    if return_models:
        return lines_model, clv_model, lines_array, planet_K, planet_R, jitter_pams

    var_dict = {
        'planet_K': planet_K,
        'planet_R': planet_R,
        'FWHM': lines_array[2,:],
        'winds': lines_array[3,:],
    }
    log_prior = 0.
    for key_name, key_vals in priors_dict.items():
        for var_value in np.atleast_1d(var_dict[key_name]):
            log_prior +=  (-(var_value - key_vals[0]) ** 2 / (2 * key_vals[1] ** 2) - 0.5 * np.log(2*np.pi) - np.log(key_vals[1]))

    """ the last value is always the jitter """
    env = 1.0 / (jitter_array ** 2.0 + ferr_res ** 2.0)
    return log_prior -0.5 * (np.size(flux_res) * np.log(2 * np.pi) +
                    np.sum((flux_res) ** 2 * env - np.log(env)))


""" case 13: more than one spectral lines, all lines are affected by the same wind and have same FWHM """
def logprob_case13(theta,
            boundaries,
            wave_meshgrid,
            transmission_spec,
            transmission_spec_err,
            clv_rm_radius,
            clv_rm_grid,
            planet_RVsinusoid,
            lines_center,
            jitter_index,
            priors_dict,
            return_models=False
            ):

    """ check the boundaries """
    if ((theta < boundaries[:,0]) | (theta > boundaries[:,1])).any():
        return -np.inf

    """ unfolding jitter parameters """
    if jitter_index is None:
        i_j = -1
        jitter_array = 0.
        jitter_pams = 0.
    elif len(jitter_index) > 0:
        jitter_array = jitter_index * 0.
        n_jitter = int(np.amax(jitter_index) + 1)
        jitter_pams = np.empty(n_jitter)
        for i_j in range(0, n_jitter):
            sel = (jitter_index == i_j)
            jitter_array[sel] = theta[-1-i_j]
            jitter_pams[i_j] = theta[-1-i_j]
    else:
        i_j = 0
        jitter_array = wave_meshgrid*0. + theta[-1]
        jitter_pams = theta[-1]

    """ the second-last value after jitter is always the semi-amplitude of the planet
        the third-last value after jitter is always the planetary radius, if included in the model
    """
    planet_K = theta[-2-i_j]
    planet_R = theta[-3-i_j]

    """ computing interpolated model spectrum
    """
    clv_model = interpolate2d_grid_nocheck(planet_R, clv_rm_radius, clv_rm_grid)

    """ line_array is always structured in this way:
        line_array = np.empty(n_pams, n_lines)
        line_array[0, 0] = wavelength
        line_array[1, 0] = contrast
        line_array[2, 0] = FWHM
        line_array[3, 0] = winds
    """

    n_lines = len(lines_center)
    lines_array = np.empty([4, n_lines])
    i_pams = 0
    for ii in range(0, n_lines):
        lines_array[0, ii] = lines_center[ii]
        lines_array[1, ii] = theta[i_pams]
        i_pams += 1
        lines_array[2, ii] = theta[-5-i_j]
        lines_array[3, ii] = theta[-4-i_j]

    lines_model = compute_multiple_lines(wave_meshgrid, planet_RVsinusoid, planet_K, lines_array, n_lines)

    flux_res = transmission_spec / clv_model / lines_model - 1.
    ferr_res = transmission_spec_err / clv_model / lines_model

    if return_models:
        return lines_model, clv_model, lines_array, planet_K, planet_R, jitter_pams

    var_dict = {
        'planet_K': planet_K,
        'planet_R': planet_R,
        'FWHM': lines_array[2,:],
        'winds': lines_array[3,:],
    }
    log_prior = 0.
    for key_name, key_vals in priors_dict.items():
        for var_value in np.atleast_1d(var_dict[key_name]):
            log_prior +=  (-(var_value - key_vals[0]) ** 2 / (2 * key_vals[1] ** 2) - 0.5 * np.log(2*np.pi) - np.log(key_vals[1]))

    """ the last value is always the jitter """
    env = 1.0 / (jitter_array ** 2.0 + ferr_res ** 2.0)
    return log_prior -0.5 * (np.size(flux_res) * np.log(2 * np.pi) +
                    np.sum((flux_res) ** 2 * env - np.log(env)))


""" case 14: more than one spectral lines, no winds """
def logprob_case14(theta,
            boundaries,
            wave_meshgrid,
            transmission_spec,
            transmission_spec_err,
            clv_rm_radius,
            clv_rm_grid,
            planet_RVsinusoid,
            lines_center,
            jitter_index,
            priors_dict,
            return_models=False
            ):

    """ check the boundaries """
    if ((theta < boundaries[:,0]) | (theta > boundaries[:,1])).any():
        return -np.inf

    """ unfolding jitter parameters """
    if jitter_index is None:
        i_j = -1
        jitter_array = 0.
        jitter_pams = 0.
    elif len(jitter_index) > 0:
        jitter_array = jitter_index * 0.
        n_jitter = int(np.amax(jitter_index) + 1)
        jitter_pams = np.empty(n_jitter)
        for i_j in range(0, n_jitter):
            sel = (jitter_index == i_j)
            jitter_array[sel] = theta[-1-i_j]
            jitter_pams[i_j] = theta[-1-i_j]
    else:
        i_j = 0
        jitter_array = wave_meshgrid*0. + theta[-1]
        jitter_pams = theta[-1]

    """ the second-last value after jitter is always the semi-amplitude of the planet
        the third-last value after jitter is always the planetary radius, if included in the model
    """
    planet_K = theta[-2-i_j]
    planet_R = theta[-3-i_j]

    """ computing interpolated model spectrum
    """
    clv_model = interpolate2d_grid_nocheck(planet_R, clv_rm_radius, clv_rm_grid)

    """ line_array is always structured in this way:
        line_array = np.empty(n_pams, n_lines)
        line_array[0, 0] = wavelength
        line_array[1, 0] = contrast
        line_array[2, 0] = FWHM
        line_array[3, 0] = winds
    """

    n_lines = len(lines_center)
    lines_array = np.empty([4, n_lines])
    i_pams = 0
    for ii in range(0, n_lines):
        lines_array[0, ii] = lines_center[ii]
        lines_array[1, ii] = theta[i_pams]
        i_pams += 1
        lines_array[2, ii] = theta[i_pams]
        i_pams += 1
        lines_array[3, ii] = 0.000

    lines_model = compute_multiple_lines(wave_meshgrid, planet_RVsinusoid, planet_K, lines_array, n_lines)

    flux_res = transmission_spec / clv_model / lines_model - 1.
    ferr_res = transmission_spec_err / clv_model / lines_model

    if return_models:
        return lines_model, clv_model, lines_array, planet_K, planet_R, jitter_pams

    var_dict = {
        'planet_K': planet_K,
        'planet_R': planet_R,
        'FWHM': lines_array[2,:],
        'winds': lines_array[3,:],
    }
    log_prior = 0.
    for key_name, key_vals in priors_dict.items():
        for var_value in np.atleast_1d(var_dict[key_name]):
            log_prior +=  (-(var_value - key_vals[0]) ** 2 / (2 * key_vals[1] ** 2) - 0.5 * np.log(2*np.pi) - np.log(key_vals[1]))

    """ the last value is always the jitter """
    env = 1.0 / (jitter_array ** 2.0 + ferr_res ** 2.0)
    return log_prior -0.5 * (np.size(flux_res) * np.log(2 * np.pi) +
                    np.sum((flux_res) ** 2 * env - np.log(env)))


""" case 15: more than one spectral lines, no winds, all lines have same FWHM """
def logprob_case15(theta,
            boundaries,
            wave_meshgrid,
            transmission_spec,
            transmission_spec_err,
            clv_rm_radius,
            clv_rm_grid,
            planet_RVsinusoid,
            lines_center,
            jitter_index,
            priors_dict,
            return_models=False
            ):

    """ check the boundaries """
    if ((theta < boundaries[:,0]) | (theta > boundaries[:,1])).any():
        return -np.inf

    """ unfolding jitter parameters """
    if jitter_index is None:
        i_j = -1
        jitter_array = 0.
        jitter_pams = 0.
    elif len(jitter_index) > 0:
        jitter_array = jitter_index * 0.
        n_jitter = int(np.amax(jitter_index) + 1)
        jitter_pams = np.empty(n_jitter)
        for i_j in range(0, n_jitter):
            sel = (jitter_index == i_j)
            jitter_array[sel] = theta[-1-i_j]
            jitter_pams[i_j] = theta[-1-i_j]
    else:
        i_j = 0
        jitter_array = wave_meshgrid*0. + theta[-1]
        jitter_pams = theta[-1]

    """ the second-last value after jitter is always the semi-amplitude of the planet
        the third-last value after jitter is always the planetary radius, if included in the model
    """
    planet_K = theta[-2-i_j]
    planet_R = theta[-3-i_j]

    """ computing interpolated model spectrum
    """
    clv_model = interpolate2d_grid_nocheck(planet_R, clv_rm_radius, clv_rm_grid)

    """ line_array is always structured in this way:
        line_array = np.empty(n_pams, n_lines)
        line_array[0, 0] = wavelength
        line_array[1, 0] = contrast
        line_array[2, 0] = FWHM
        line_array[3, 0] = winds
    """

    n_lines = len(lines_center)
    lines_array = np.empty([4, n_lines])
    i_pams = 0
    for ii in range(0, n_lines):
        lines_array[0, ii] = lines_center[ii]
        lines_array[1, ii] = theta[i_pams]
        i_pams += 1
        lines_array[2, ii] = theta[-4-i_j]
        lines_array[3, ii] = 0.000

    lines_model = compute_multiple_lines(wave_meshgrid, planet_RVsinusoid, planet_K, lines_array, n_lines)

    flux_res = transmission_spec / clv_model / lines_model - 1.
    ferr_res = transmission_spec_err / clv_model / lines_model

    if return_models:
        return lines_model, clv_model, lines_array, planet_K, planet_R, jitter_pams

    var_dict = {
        'planet_K': planet_K,
        'planet_R': planet_R,
        'FWHM': lines_array[2,:],
        'winds': lines_array[3,:],
    }
    log_prior = 0.
    for key_name, key_vals in priors_dict.items():
        for var_value in np.atleast_1d(var_dict[key_name]):
            log_prior +=  (-(var_value - key_vals[0]) ** 2 / (2 * key_vals[1] ** 2) - 0.5 * np.log(2*np.pi) - np.log(key_vals[1]))

    """ the last value is always the jitter """
    env = 1.0 / (jitter_array ** 2.0 + ferr_res ** 2.0)
    return log_prior -0.5 * (np.size(flux_res) * np.log(2 * np.pi) +
                    np.sum((flux_res) ** 2 * env - np.log(env)))






""" case 20: more than one spectral lines, no Rp dependance, all line parameters are free and independent """
def logprob_case20(theta,
            boundaries,
            wave_meshgrid,
            transmission_spec,
            transmission_spec_err,
            clv_model,
            planet_RVsinusoid,
            lines_center,
            jitter_index,
            priors_dict,
            return_models=False
            ):

    """ check the boundaries """
    if ((theta < boundaries[:,0]) | (theta > boundaries[:,1])).any():
        return -np.inf

    """ unfolding jitter parameters """
    if jitter_index is None:
        i_j = -1
        jitter_array = 0.
        jitter_pams = 0.
    elif len(jitter_index) > 0:
        jitter_array = jitter_index * 0.
        n_jitter = int(np.amax(jitter_index) + 1)
        jitter_pams = np.empty(n_jitter)
        for i_j in range(0, n_jitter):
            sel = (jitter_index == i_j)
            jitter_array[sel] = theta[-1-i_j]
            jitter_pams[i_j] = theta[-1-i_j]
    else:
        i_j = 0
        jitter_array = wave_meshgrid*0. + theta[-1]
        jitter_pams = theta[-1]

    """ the second-last value after jitter is always the semi-amplitude of the planet
        the third-last value after jitter is always the planetary radius, if included in the model
    """
    planet_K = theta[-2-i_j]

    """ line_array is always structured in this way:
        line_array = np.empty(n_pams, n_lines)
        line_array[0, 0] = wavelength
        line_array[1, 0] = contrast
        line_array[2, 0] = FWHM
        line_array[3, 0] = winds
    """

    n_lines = len(lines_center)
    lines_array = np.empty([4, n_lines])
    i_pams = 0
    for ii in range(0, n_lines):
        lines_array[0, ii] = lines_center[ii]
        lines_array[1, ii] = theta[i_pams]
        i_pams += 1
        lines_array[2, ii] = theta[i_pams]
        i_pams += 1
        lines_array[3, ii] = theta[i_pams]
        i_pams += 1

    lines_model = compute_multiple_lines(wave_meshgrid, planet_RVsinusoid, planet_K, lines_array, n_lines)

    flux_res = transmission_spec / clv_model / lines_model - 1.
    ferr_res = transmission_spec_err / clv_model / lines_model

    if return_models:
        return lines_model, clv_model, lines_array, planet_K, 1.0000, jitter_pams

    var_dict = {
        'planet_K': planet_K,
        'planet_R': 1.00000,
        'FWHM': lines_array[2,:],
        'winds': lines_array[3,:],
    }
    log_prior = 0.
    for key_name, key_vals in priors_dict.items():
        for var_value in np.atleast_1d(var_dict[key_name]):
            log_prior +=  (-(var_value - key_vals[0]) ** 2 / (2 * key_vals[1] ** 2) - 0.5 * np.log(2*np.pi) - np.log(key_vals[1]))

    """ the last value is always the jitter """
    env = 1.0 / (jitter_array ** 2.0 + ferr_res ** 2.0)
    return log_prior -0.5 * (np.size(flux_res) * np.log(2 * np.pi) +
                    np.sum((flux_res) ** 2 * env - np.log(env)))


""" case 21: more than one spectral lines, no Rp dependance, all lines are affected by the same wind """
def logprob_case21(theta,
            boundaries,
            wave_meshgrid,
            transmission_spec,
            transmission_spec_err,
            clv_model,
            planet_RVsinusoid,
            lines_center,
            jitter_index,
            priors_dict,
            return_models=False
            ):

    """ check the boundaries """
    if ((theta < boundaries[:,0]) | (theta > boundaries[:,1])).any():
        return -np.inf

    """ unfolding jitter parameters """
    if jitter_index is None:
        i_j = -1
        jitter_array = 0.
        jitter_pams = 0.
    elif len(jitter_index) > 0:
        jitter_array = jitter_index * 0.
        n_jitter = int(np.amax(jitter_index) + 1)
        jitter_pams = np.empty(n_jitter)
        for i_j in range(0, n_jitter):
            sel = (jitter_index == i_j)
            jitter_array[sel] = theta[-1-i_j]
            jitter_pams[i_j] = theta[-1-i_j]
    else:
        i_j = 0
        jitter_array = wave_meshgrid*0. + theta[-1]
        jitter_pams = theta[-1]

    """ the second-last value after jitter is always the semi-amplitude of the planet
        the third-last value after jitter is always the planetary radius, if included in the model
    """
    planet_K = theta[-2-i_j]

    """ line_array is always structured in this way:
        line_array = np.empty(n_pams, n_lines)
        line_array[0, 0] = wavelength
        line_array[1, 0] = contrast
        line_array[2, 0] = FWHM
        line_array[3, 0] = winds
    """

    n_lines = len(lines_center)
    lines_array = np.empty([4, n_lines])
    i_pams = 0
    for ii in range(0, n_lines):
        lines_array[0, ii] = lines_center[ii]
        lines_array[1, ii] = theta[i_pams]
        i_pams += 1
        lines_array[2, ii] = theta[i_pams]
        i_pams += 1
        lines_array[3, ii] = theta[-3-i_j]

    lines_model = compute_multiple_lines(wave_meshgrid, planet_RVsinusoid, planet_K, lines_array, n_lines)

    flux_res = transmission_spec / clv_model / lines_model - 1.
    ferr_res = transmission_spec_err / clv_model / lines_model

    if return_models:
        return lines_model, clv_model, lines_array, planet_K, 1.0000, jitter_pams

    var_dict = {
        'planet_K': planet_K,
        'planet_R': 1.000000,
        'FWHM': lines_array[2,:],
        'winds': lines_array[3,:],
    }
    log_prior = 0.
    for key_name, key_vals in priors_dict.items():
        for var_value in np.atleast_1d(var_dict[key_name]):
            log_prior +=  (-(var_value - key_vals[0]) ** 2 / (2 * key_vals[1] ** 2) - 0.5 * np.log(2*np.pi) - np.log(key_vals[1]))

    """ the last value is always the jitter """
    env = 1.0 / (jitter_array ** 2.0 + ferr_res ** 2.0)
    return log_prior -0.5 * (np.size(flux_res) * np.log(2 * np.pi) +
                    np.sum((flux_res) ** 2 * env - np.log(env)))


""" case 22: more than one spectral lines, no Rp dependance, all lines have same FWHM """
def logprob_case22(theta,
            boundaries,
            wave_meshgrid,
            transmission_spec,
            transmission_spec_err,
            clv_model,
            planet_RVsinusoid,
            lines_center,
            jitter_index,
            priors_dict,
            return_models=False
            ):

    """ check the boundaries """
    if ((theta < boundaries[:,0]) | (theta > boundaries[:,1])).any():
        return -np.inf

    """ unfolding jitter parameters """
    if jitter_index is None:
        i_j = -1
        jitter_array = 0.
        jitter_pams = 0.
    elif len(jitter_index) > 0:
        jitter_array = jitter_index * 0.
        n_jitter = int(np.amax(jitter_index) + 1)
        jitter_pams = np.empty(n_jitter)
        for i_j in range(0, n_jitter):
            sel = (jitter_index == i_j)
            jitter_array[sel] = theta[-1-i_j]
            jitter_pams[i_j] = theta[-1-i_j]
    else:
        i_j = 0
        jitter_array = wave_meshgrid*0. + theta[-1]
        jitter_pams = theta[-1]

    """ the second-last value after jitter is always the semi-amplitude of the planet
        the third-last value after jitter is always the planetary radius, if included in the model
    """
    planet_K = theta[-2-i_j]

    """ line_array is always structured in this way:
        line_array = np.empty(n_pams, n_lines)
        line_array[0, 0] = wavelength
        line_array[1, 0] = contrast
        line_array[2, 0] = FWHM
        line_array[3, 0] = winds
    """

    n_lines = len(lines_center)
    lines_array = np.empty([4, n_lines])
    i_pams = 0
    for ii in range(0, n_lines):
        lines_array[0, ii] = lines_center[ii]
        lines_array[1, ii] = theta[i_pams]
        i_pams += 1
        lines_array[2, ii] = theta[-3-i_j]
        lines_array[3, ii] = theta[i_pams]
        i_pams += 1

    lines_model = compute_multiple_lines(wave_meshgrid, planet_RVsinusoid, planet_K, lines_array, n_lines)

    flux_res = transmission_spec / clv_model / lines_model - 1.
    ferr_res = transmission_spec_err / clv_model / lines_model

    if return_models:
        return lines_model, clv_model, lines_array, planet_K, 1.0000, jitter_pams

    var_dict = {
        'planet_K': planet_K,
        'planet_R': 1.000000,
        'FWHM': lines_array[2,:],
        'winds': lines_array[3,:],
    }
    log_prior = 0.
    for key_name, key_vals in priors_dict.items():
        for var_value in np.atleast_1d(var_dict[key_name]):
            log_prior +=  (-(var_value - key_vals[0]) ** 2 / (2 * key_vals[1] ** 2) - 0.5 * np.log(2*np.pi) - np.log(key_vals[1]))

    """ the last value is always the jitter """
    env = 1.0 / (jitter_array ** 2.0 + ferr_res ** 2.0)
    return log_prior -0.5 * (np.size(flux_res) * np.log(2 * np.pi) +
                    np.sum((flux_res) ** 2 * env - np.log(env)))


""" case 23: more than one spectral lines, no Rp dependance, all lines are affected by the same wind and have same FWHM """
def logprob_case23(theta,
            boundaries,
            wave_meshgrid,
            transmission_spec,
            transmission_spec_err,
            clv_model,
            planet_RVsinusoid,
            lines_center,
            jitter_index,
            priors_dict,
            return_models=False
            ):

    """ check the boundaries """
    if ((theta < boundaries[:,0]) | (theta > boundaries[:,1])).any():
        return -np.inf

    """ unfolding jitter parameters """
    if jitter_index is None:
        i_j = -1
        jitter_array = 0.
        jitter_pams = 0.
    elif len(jitter_index) > 0:
        jitter_array = jitter_index * 0.
        n_jitter = int(np.amax(jitter_index) + 1)
        jitter_pams = np.empty(n_jitter)
        for i_j in range(0, n_jitter):
            sel = (jitter_index == i_j)
            jitter_array[sel] = theta[-1-i_j]
            jitter_pams[i_j] = theta[-1-i_j]
    else:
        i_j = 0
        jitter_array = wave_meshgrid*0. + theta[-1]
        jitter_pams = theta[-1]

    """ the second-last value after jitter is always the semi-amplitude of the planet
        the third-last value after jitter is always the planetary radius, if included in the model
    """
    planet_K = theta[-2-i_j]

    """ line_array is always structured in this way:
        line_array = np.empty(n_pams, n_lines)
        line_array[0, 0] = wavelength
        line_array[1, 0] = contrast
        line_array[2, 0] = FWHM
        line_array[3, 0] = winds
    """

    n_lines = len(lines_center)
    lines_array = np.empty([4, n_lines])
    i_pams = 0

    for ii in range(0, n_lines):
        lines_array[0, ii] = lines_center[ii]
        lines_array[1, ii] = theta[i_pams]
        i_pams += 1
        lines_array[2, ii] = theta[-4-i_j]
        lines_array[3, ii] = theta[-3-i_j]

    lines_model = compute_multiple_lines(wave_meshgrid, planet_RVsinusoid, planet_K, lines_array, n_lines)

    flux_res = transmission_spec / clv_model / lines_model - 1.
    ferr_res = transmission_spec_err / clv_model / lines_model

    if return_models:
        return lines_model, clv_model, lines_array, planet_K, 1.0000, jitter_pams

    var_dict = {
        'planet_K': planet_K,
        'planet_R': 1.000000,
        'FWHM': lines_array[2,:],
        'winds': lines_array[3,:],
    }
    log_prior = 0.
    for key_name, key_vals in priors_dict.items():
        for var_value in np.atleast_1d(var_dict[key_name]):
            log_prior +=  (-(var_value - key_vals[0]) ** 2 / (2 * key_vals[1] ** 2) - 0.5 * np.log(2*np.pi) - np.log(key_vals[1]))

    """ the last value is always the jitter """
    env = 1.0 / (jitter_array ** 2.0 + ferr_res ** 2.0)
    return log_prior -0.5 * (np.size(flux_res) * np.log(2 * np.pi) +
                    np.sum((flux_res) ** 2 * env - np.log(env)))


""" case 24: more than one spectral lines, no Rp dependance, no winds """
def logprob_case24(theta,
            boundaries,
            wave_meshgrid,
            transmission_spec,
            transmission_spec_err,
            clv_model,
            planet_RVsinusoid,
            lines_center,
            jitter_index,
            priors_dict,
            return_models=False
            ):

    """ check the boundaries """
    if ((theta < boundaries[:,0]) | (theta > boundaries[:,1])).any():
        return -np.inf

    """ unfolding jitter parameters """
    if jitter_index is None:
        i_j = -1
        jitter_array = 0.
        jitter_pams = 0.
    elif len(jitter_index) > 0:
        jitter_array = jitter_index * 0.
        n_jitter = int(np.amax(jitter_index) + 1)
        jitter_pams = np.empty(n_jitter)
        for i_j in range(0, n_jitter):
            sel = (jitter_index == i_j)
            jitter_array[sel] = theta[-1-i_j]
            jitter_pams[i_j] = theta[-1-i_j]
    else:
        i_j = 0
        jitter_array = wave_meshgrid*0. + theta[-1]
        jitter_pams = theta[-1]

    """ the second-last value after jitter is always the semi-amplitude of the planet
        the third-last value after jitter is always the planetary radius, if included in the model
    """
    planet_K = theta[-2-i_j]

    """ line_array is always structured in this way:
        line_array = np.empty(n_pams, n_lines)
        line_array[0, 0] = wavelength
        line_array[1, 0] = contrast
        line_array[2, 0] = FWHM
        line_array[3, 0] = winds
    """

    n_lines = len(lines_center)
    lines_array = np.empty([4, n_lines])
    i_pams = 0
    for ii in range(0, n_lines):
        lines_array[0, ii] = lines_center[ii]
        lines_array[1, ii] = theta[i_pams]
        i_pams += 1
        lines_array[2, ii] = theta[i_pams]
        i_pams += 1
        lines_array[3, ii] = 0.000

    lines_model = compute_multiple_lines(wave_meshgrid, planet_RVsinusoid, planet_K, lines_array, n_lines)

    flux_res = transmission_spec / clv_model / lines_model - 1.
    ferr_res = transmission_spec_err / clv_model / lines_model

    if return_models:
        return lines_model, clv_model, lines_array, planet_K, 1.0000, jitter_pams

    var_dict = {
        'planet_K': planet_K,
        'planet_R': 1.000000,
        'FWHM': lines_array[2,:],
        'winds': lines_array[3,:],
    }
    log_prior = 0.
    for key_name, key_vals in priors_dict.items():
        for var_value in np.atleast_1d(var_dict[key_name]):
            log_prior +=  (-(var_value - key_vals[0]) ** 2 / (2 * key_vals[1] ** 2) - 0.5 * np.log(2*np.pi) - np.log(key_vals[1]))

    """ the last value is always the jitter """
    env = 1.0 / (jitter_array ** 2.0 + ferr_res ** 2.0)
    return log_prior -0.5 * (np.size(flux_res) * np.log(2 * np.pi) +
                    np.sum((flux_res) ** 2 * env - np.log(env)))


""" case 25: more than one spectral lines, no Rp dependance, no winds, all lines have same FWHM """
def logprob_case25(theta,
            boundaries,
            wave_meshgrid,
            transmission_spec,
            transmission_spec_err,
            clv_model,
            planet_RVsinusoid,
            lines_center,
            jitter_index,
            priors_dict,
            return_models=False
            ):

    """ check the boundaries """
    if ((theta < boundaries[:,0]) | (theta > boundaries[:,1])).any():
        return -np.inf

    """ unfolding jitter parameters """
    if jitter_index is None:
        i_j = -1
        jitter_array = 0.
        jitter_pams = 0.
    elif len(jitter_index) > 0:
        jitter_array = jitter_index * 0.
        n_jitter = int(np.amax(jitter_index) + 1)
        jitter_pams = np.empty(n_jitter)
        for i_j in range(0, n_jitter):
            sel = (jitter_index == i_j)
            jitter_array[sel] = theta[-1-i_j]
            jitter_pams[i_j] = theta[-1-i_j]
    else:
        i_j = 0
        jitter_array = wave_meshgrid*0. + theta[-1]
        jitter_pams = theta[-1]

    """ the second-last value after jitter is always the semi-amplitude of the planet
        the third-last value after jitter is always the planetary radius, if included in the model
    """
    planet_K = theta[-2-i_j]

    """ line_array is always structured in this way:
        line_array = np.empty(n_pams, n_lines)
        line_array[0, 0] = wavelength
        line_array[1, 0] = contrast
        line_array[2, 0] = FWHM
        line_array[3, 0] = winds
    """

    n_lines = len(lines_center)
    lines_array = np.empty([4, n_lines])
    i_pams = 0
    for ii in range(0, n_lines):
        lines_array[0, ii] = lines_center[ii]
        lines_array[1, ii] = theta[i_pams]
        i_pams += 1
        lines_array[2, ii] = theta[-3-i_j]
        lines_array[3, ii] = 0.000

    lines_model = compute_multiple_lines(wave_meshgrid, planet_RVsinusoid, planet_K, lines_array, n_lines)

    flux_res = transmission_spec / clv_model / lines_model - 1.
    ferr_res = transmission_spec_err / clv_model / lines_model

    if return_models:
        return lines_model, clv_model, lines_array, planet_K, 1.0000, jitter_pams

    var_dict = {
        'planet_K': planet_K,
        'planet_R': 1.000000,
        'FWHM': lines_array[2,:],
        'winds': lines_array[3,:],
    }
    log_prior = 0.
    for key_name, key_vals in priors_dict.items():
        for var_value in np.atleast_1d(var_dict[key_name]):
            log_prior +=  (-(var_value - key_vals[0]) ** 2 / (2 * key_vals[1] ** 2) - 0.5 * np.log(2*np.pi) - np.log(key_vals[1]))

    """ the last value is always the jitter """
    env = 1.0 / (jitter_array ** 2.0 + ferr_res ** 2.0)
    return log_prior -0.5 * (np.size(flux_res) * np.log(2 * np.pi) +
                    np.sum((flux_res) ** 2 * env - np.log(env)))
