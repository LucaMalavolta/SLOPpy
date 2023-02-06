import numpy as np
from scipy.optimize import curve_fit


def interpolate1d_grid_nocheck(val, arr_1, arr_2):
    # using enumerate() + next() to find index of
    # first element just greater than 0.6
    res = next(x for x, v in enumerate(arr_1)  if v > val)
    interp_out =  (val-arr_1[res-1])/(arr_1[res]-arr_1[res-1])*(arr_2[res,:]-arr_2[res-1,:])  + arr_2[res-1,:]
    return interp_out

def interpolate2d_grid_nocheck(val, arr_1, arr_2):
    # using enumerate() + next() to find index of
    # first element just greater than 0.6
    res = next(x for x, v in enumerate(arr_1)  if v > val)
    interp_out =  (val-arr_1[res-1])/(arr_1[res]-arr_1[res-1])*(arr_2[res,:,:]-arr_2[res-1,:,:])  + arr_2[res-1,:,:]
    return interp_out

def first_derivative(x_arr, y_arr):
    n_points = len(x_arr)
    derivative = np.zeros(n_points, dtype=np.double)

    derivative[1:-1] = (y_arr[2:] - y_arr[:-2]) / (x_arr[2:] - x_arr[:-2])
    derivative[0] = derivative[1]
    derivative[-1] = derivative[-2]

    return derivative
