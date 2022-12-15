from __future__ import print_function, division
from SLOPpy.subroutines.common import *
from SLOPpy.subroutines.io_subroutines import *
from SLOPpy.transmission_spectrum import *

#__all__ = ['compute_transmission_spectrum_planetRF_iterative',
#           'plot_transmission_spectrum_planetRF_iterative',
#           'compute_transmission_spectrum_stellarRF_iterative',
#           'plot_transmission_spectrum_stellarRF_iterative',
#           'compute_transmission_spectrum_observerRF_iterative',
#           'plot_transmission_spectrum_observerRF_iterative',
#           'compute_transmission_spectrum_iterative',
#           'plot_transmission_spectrum_iterative']


def compute_transmission_spectrum_planetRF(config_in, lines_label):
    compute_transmission_spectrum(config_in, lines_label, reference='planetRF')


def plot_transmission_spectrum_planetRF(config_in, lines_label, night_input, results_input=''):
    plot_transmission_spectrum(config_in, lines_label, night_input, results_input, reference='planetRF')


def compute_transmission_spectrum_stellarRF(config_in, lines_label):
    compute_transmission_spectrum(config_in, lines_label, reference='stellarRF')


def plot_transmission_spectrum_stellarRF(config_in, lines_label, night_input, results_input=''):
    plot_transmission_spectrum(config_in, lines_label, night_input, results_input, reference='stellarRF')


def compute_transmission_spectrum_observerRF(config_in, lines_label):
    compute_transmission_spectrum(config_in, lines_label, reference='observerRF')


def plot_transmission_spectrum_observerRF(config_in, lines_label, night_input, results_input=''):
    plot_transmission_spectrum(config_in, lines_label, night_input, results_input, reference='observerRF')






def compute_transmission_spectrum_planetRF_iterative(config_in, lines_label):

    pca_parameters = from_config_get_pca_parameters(config_in)
    for it in range(0,pca_parameters.get('iterations',5)):
        compute_transmission_spectrum(config_in, lines_label, reference='planetRF', pca_iteration=it)

def compute_transmission_spectrum_stellarRF_iterative(config_in, lines_label):

    pca_parameters = from_config_get_pca_parameters(config_in)
    for it in range(0,pca_parameters.get('iterations',5)):
        compute_transmission_spectrum(config_in, lines_label, reference='stellarRF', pca_iteration=it)

def compute_transmission_spectrum_observerRF_iterative(config_in, lines_label):

    pca_parameters = from_config_get_pca_parameters(config_in)
    for it in range(0,pca_parameters.get('iterations',5)):
        compute_transmission_spectrum(config_in, lines_label, reference='observerRF', pca_iteration=it)

def compute_transmission_spectrum_iterative(config_in, lines_label):

    pca_parameters = from_config_get_pca_parameters(config_in)
    for it in range(0,pca_parameters.get('iterations',5)):
        compute_transmission_spectrum(config_in, lines_label, reference='planetRF', pca_iteration=it)


def  plot_transmission_spectrum_planetRF_iterative(config_in, lines_label, night_input, results_input=''):

    pca_parameters = from_config_get_pca_parameters(config_in)
    for it in range(0,pca_parameters.get('iterations',5)):
        plot_transmission_spectrum(config_in, lines_label, night_input, results_input, reference='planetRF', pca_iteration=it)

def plot_transmission_spectrum_stellarRF_iterative(config_in, lines_label, night_input, results_input=''):

    pca_parameters = from_config_get_pca_parameters(config_in)
    for it in range(0,pca_parameters.get('iterations',5)):
        plot_transmission_spectrum(config_in, lines_label, night_input, results_input, reference='stellarRF', pca_iteration=it)

def plot_transmission_spectrum_observerRF_iterative(config_in, lines_label, night_input, results_input=''):

    pca_parameters = from_config_get_pca_parameters(config_in)
    for it in range(0,pca_parameters.get('iterations',5)):
        plot_transmission_spectrum(config_in, lines_label, night_input, results_input, reference='observerRF', pca_iteration=it)

def plot_transmission_spectrum_iterative(config_in, lines_label, night_input, results_input=''):

    pca_parameters = from_config_get_pca_parameters(config_in)
    for it in range(0,pca_parameters.get('iterations',5)):
        plot_transmission_spectrum(config_in, lines_label, night_input, results_input, reference='planetRF', pca_iteration=it)


