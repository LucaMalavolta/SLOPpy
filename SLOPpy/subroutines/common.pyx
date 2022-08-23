import sys
import os
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import Normalize
import numpy as np
import scipy.stats as sci_stats
import scipy.optimize as sci_optimize
import scipy.interpolate as sci_int
from astropy.io import fits
import json
import warnings
import pygtc
import re

# Use ordered dictionaries for the observations
from collections import OrderedDict


"""
List of exceptions 
"""
class MissingFileException(Exception):
    pass

class MissingKeywordException(Exception):
    pass

class OutOfRangeException(Exception):
    pass

