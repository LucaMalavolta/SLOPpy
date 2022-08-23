from __future__ import print_function, division
import numpy as np
from scipy.linalg import lstsq
from scipy.optimize import curve_fit
from sklearn import linear_model, datasets


def berv_telluric_curvefit(xdata, p0, p1, p2):
    return xdata[0] * p0 + xdata[1] * p1 + p2


def berv_linear_curve_fit(airmass, berv, logi_array, sigi_array, n_axis):
    C = []
    pams = np.zeros(3)-0.1
    ltel = np.empty(n_axis)
    shift = np.empty(n_axis)
    zero = np.empty(n_axis)
    airmass_zero = 0. #np.average(airmass)
    berv_zero = 0. #np.average(berv)

    for ii in xrange(0, n_axis):
        popt, pcov = curve_fit(berv_telluric_curvefit,
                               [airmass-airmass_zero, berv-berv_zero],
                               logi_array[:, ii],
                               p0=pams,
                               sigma = sigi_array[:, ii],
                               bounds=([-np.inf, -np.inf, -np.inf], [0.000, np.inf, np.inf]))

        ltel[ii] = popt[0]
        shift[ii] = popt[1]
        zero[ii] = popt[2]

    return ltel, shift, zero


def berv_linear_lstsq(airmass, berv, logi_array):
    A = np.c_[airmass, berv, np.ones(logi_array.shape[0])]
    C, _, _, _ = lstsq(A, logi_array)    # coefficients
    return C[0], C[1], C[2]


def airmass_telluric_curvefit(xdata, p0, p1):
    return xdata * p0 + p1





def airmass_linear_curve_fit_ransac(airmass, logi_array, sigi_array, n_axis):

    pams = np.zeros(2)
    ltel = np.empty(n_axis)
    zero = np.empty(n_axis)
    airmass_zero = np.average(airmass)

    airmass_reshape = (airmass-airmass_zero).reshape(-1,1)

    ransac = linear_model.RANSACRegressor()


    for ii in range(0, n_axis):
        y_zero = np.average(logi_array[:, ii])
        ransac.fit(airmass_reshape, logi_array[:, ii]-y_zero)
        ltel[ii] = ransac.estimator_.coef_[0]
        zero[ii] = ransac.estimator_.intercept_ + y_zero

    return ltel, zero


def airmass_linear_curve_fit(airmass, logi_array, sigi_array, n_axis):
    C = []
    pams = np.zeros(2)
    ltel = np.empty(n_axis)
    zero = np.empty(n_axis)
    airmass_zero = np.average(airmass)

    for ii in range(0, n_axis):
        y_zero = np.average(logi_array[:, ii])
        popt, pcov = curve_fit(airmass_telluric_curvefit,
                               airmass-airmass_zero,
                               logi_array[:, ii]-y_zero,
                               p0=pams,
                               sigma = sigi_array[:, ii],
                               bounds=([-np.inf, -np.inf], [np.inf, np.inf]))

        ltel[ii] = popt[0]
        zero[ii] = popt[1]

    return ltel, zero


def berv_linear_curve_fit_modified(airmass, berv, logi_array, sigi_array, n_axis):
    C = []
    pams = np.zeros(3)
    ltel = np.empty(n_axis)
    shift = np.empty(n_axis)
    zero = np.empty(n_axis)
    airmass_zero = np.average(airmass)
    berv_zero = np.average(berv)

    for ii in range(0, n_axis):
        popt, pcov = curve_fit(berv_telluric_curvefit,
                               [airmass-airmass_zero, berv-berv_zero],
                               logi_array[:, ii],
                               p0=pams,
                               sigma = sigi_array[:, ii],
                               bounds=([-np.inf, -np.inf, -np.inf], [np.inf, np.inf, np.inf]))

        ltel[ii] = popt[0]
        shift[ii] = popt[1]
        zero[ii] = popt[2]

    return ltel, shift, zero



def airmass_linear_lstsq(airmass, logi_array):
    A = np.c_[airmass, np.ones(logi_array.shape[0])]
    C, _, _, _ = lstsq(A, logi_array)    # coefficients
    return C[0], C[1]



