from __future__ import print_function, division
import numpy as np
from astropy.io import fits
from SLOPpy.subroutines.rebin_subroutines import *
from SLOPpy.subroutines.constants import *
from SLOPpy.instruments.HARPN_DRSv3 import *
from SLOPpy.instruments.HARPS_DRSv3 import *
from SLOPpy.instruments.PEPSI_reduced import *
from SLOPpy.instruments.ESPRESSO import *

def get_calib_data(instrument, night_dict, archive, file_rad, fiber='A', order_selection=None):

    if instrument =='HARPS-N':
        return HARPN_DRSv3_get_calib_data(archive, file_rad, fiber=fiber, order_selection=order_selection)
    elif instrument =='HARPS':
        return HARPS_DRSv3_get_calib_data(archive, file_rad, fiber=fiber, order_selection=order_selection)
    elif instrument in ['PEPSI', 'PEPSI_red', 'PEPSI_blue']:
        return PEPSI_get_calib_data(archive, file_rad, fiber=fiber, order_selection=order_selection)
    elif instrument in ['ESPRESSO', 'ESPRESSO_odd', 'ESPRESSO_even']:
        return ESPRESSO_get_calib_data(archive, file_rad, fiber=fiber, order_selection=order_selection, night_dict=night_dict)
    else:
        raise ValueError("Instrument not supported")


def get_input_data(instrument, night_dict, archive, file_rad, mask, fiber='A', skip_ccf=None, skip_s1d=True, order_selection=None):

    if instrument =='HARPS-N':
        return HARPN_DRSv3_get_input_data(archive, file_rad, mask, fiber=fiber, skip_ccf=skip_ccf, skip_s1d=skip_s1d, order_selection=order_selection)
    elif instrument =='HARPS':
        return HARPS_DRSv3_get_input_data(archive, file_rad, mask, fiber=fiber, skip_ccf=skip_ccf, skip_s1d=skip_s1d, order_selection=order_selection)
    elif instrument in ['PEPSI', 'PEPSI_red', 'PEPSI_blue']:
        return PEPSI_get_input_data(archive, file_rad, mask, fiber=fiber, skip_ccf=skip_ccf, skip_s1d=skip_s1d, order_selection=order_selection)
    elif instrument in ['ESPRESSO', 'ESPRESSO_odd', 'ESPRESSO_even']:
        return ESPRESSO_get_input_data(archive, file_rad, mask, fiber=fiber, skip_ccf=skip_ccf, skip_s1d=skip_s1d, order_selection=order_selection, night_dict=night_dict)
    else:
        raise ValueError("Instrument not supported")
