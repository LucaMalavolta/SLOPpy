from __future__ import print_function, division
from SLOPpy.subroutines.common import *
from SLOPpy.subroutines.spectral_subroutines import *
from SLOPpy.subroutines.io_subroutines import *
from SLOPpy.subroutines.fit_subroutines import *
from SLOPpy.subroutines.plot_subroutines import *
from SLOPpy.subroutines.shortcuts import *

__all__ = ["compute_telluric_molecfit_preparation",
           "plot_telluric_molecfit_preparation"]

def compute_telluric_molecfit_preparation(config_in):
    """
    Lazy workaround
    :param config_in:
    :param kwargs:
    :return:
    """

    night_dict = from_config_get_nights(config_in)

    molecfit_dict = from_config_get_molecfit(config_in)

    for night in night_dict:

        try:
            tellprep = load_from_cpickle('telluric_molecfit_preparation', config_in['output'], night)
            continue
        except:
            print()
            print("compute_telluric_molecfit_preparation             Night: ", night)
            print()

        observational_pams = load_from_cpickle('observational_pams', config_in['output'], night)

        tellprep = {
            'work_dir':  config_in['output'] + '_molecfit_' + night,
            'include': {}
        }



        """
        We store all the molecfit files in a subdirectory
        We save the path of the main directory to a temporary file
        """

        os.system('mkdir -p ' + tellprep['work_dir'])
        os.system('mkdir -p ' + tellprep['work_dir'] + '/output/')

        """
        Creation of the include files
        """


        """
        includes_spans_ORF: wavelength ranges _with_ telluric lines, in the ORF
        includes_spans_SRF: wavelength ranges _without_ stellar lines and broadly
            overlapping with telluric ranges, in the SRF
        the two lists must have the same number of columns, with precise correspondence
        """
        tellprep['include']['spans_telluric'] = np.genfromtxt(molecfit_dict['include_telluric'])
        tellprep['include']['spans_stellar_SRF'] = np.genfromtxt(molecfit_dict['include_stellar'])

        #print()
        #print(tellprep['include']['spans_telluric'])
        #print()
        #print(tellprep['include']['spans_stellar_SRF'])

        """ shift the stellar wavelength ranges into ORF """
        tellprep['include']['rv_shift_SRF2ORF'] = -observational_pams['BERV_avg'] + observational_pams['RV_star'][
            'RV_systemic']

        #print()
        #print(tellprep['include']['rv_shift_SRF2ORF'])
        #print(observational_pams['BERV_avg'])
        #print(observational_pams['RV_star']['RV_systemic'])
        #print()
        #print()

        tellprep['include']['spans_stellar'] = tellprep['include']['spans_stellar_SRF']\
                                               * (tellprep['include']['rv_shift_SRF2ORF']
                                                  / (299792458. / 1000.000) + 1.00000)

        #print()
        #print(tellprep['include']['spans_stellar'])

        """ Selecting the overlapping regions between the two lists: we want telluric regions that are not contaminated
            by stellar lines,
        """
        sel_lower = (tellprep['include']['spans_stellar'][:, 0] > tellprep['include']['spans_telluric'][:, 0])
        sel_upper = (tellprep['include']['spans_stellar'][:, 1] < tellprep['include']['spans_telluric'][:, 1])

        """ Final list in the ORF is built"""
        tellprep['include']['selected'] = tellprep['include']['spans_telluric'].copy()
        tellprep['include']['selected'][sel_lower, 0] = tellprep['include']['spans_stellar'][sel_lower, 0]
        tellprep['include']['selected'][sel_upper, 1] = tellprep['include']['spans_stellar'][sel_upper, 1]

        #print()
        #print(tellprep['include']['selected'])

        """ Molecfit line list must be given in vacuum wavelength, even if the stellar spectra is in air wavelength
            conversion from air to vacuum for include file preparation
            where s = 10000 / lambda air and the conversion is: lambda_vac = lambda_air * n.
            http://www.astro.uu.se/valdwiki/Air-to-vacuum%20conversion
        """
        s2 = (10000. / tellprep['include']['selected']) ** 2
        n = 1 + 0.00008336624212083 + 0.02408926869968 / (130.1065924522 - s2) + 0.0001599740894897 / (
                38.92568793293 - s2)
        tellprep['include']['vacuum'] = tellprep['include']['selected'] * n / 10000.

        #fileout = open('./' + tellprep['work_dir'] + '/include_' + night + '.dat', 'w')
        #for i_s, i_e in zip(tellprep['include']['vacuum'][:, 0], tellprep['include']['vacuum'][:, 1]):
        #    fileout.write('{0:12.8f} {1:12.8f}\n'.format(i_s, i_e))
        #fileout.close()


        #quit()
        save_to_cpickle('telluric_molecfit_preparation', tellprep, config_in['output'], night)



def plot_telluric_molecfit_preparation(config_in, night_input=''):
    import matplotlib.pyplot as plt

    night_dict = from_config_get_nights(config_in)
    instrument_dict = from_config_get_instrument(config_in)
    system_dict = from_config_get_system(config_in)

    if night_input == '':
        night_list = night_dict
    else:
        night_list = np.atleast_1d(night_input)

    for night in night_list:
        print("plot_telluric_template                     Night: ", night)
