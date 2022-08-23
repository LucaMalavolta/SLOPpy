from __future__ import print_function, division
import numpy as np
import SLOPpy.subroutines.kepler_exo as kp

__all__ = ["StarParameters", "PlanetParameters"]


G_grav = 6.67428e-11 # Gravitational Constants in SI system [m^3/kg/s^2]
M_sun = 1.9884e30 # Value from TRADES

class StarParameters:
    """
    This class is just a stub with basic properties, to be expanded in th efuture
    E.g.: stellar mass, radius, activity level, etc
    """
    def __init__(self, mass=None, radius=None):
        self.kind = 'star'
        self.mass = mass
        self.radius = radius


class PlanetParameters:
    """
    This class is just a stub with basic properties, to be expanded in the future
    E.g.: give a planetary mass and convert it to expected K, and vice versa
    """
    def __init__(self):
        #array of shape (2), error should be stored here
        self.kind = 'planet'
        self.star_mass = None
        self.star_radius = None

        self.reference_epoch = None
        self.reference_transit = None

        self.period = None
        self.phase = None
        self.eccentricity = None
        self.omega = None

        self.RV_semiamplitude = None
        self.RVplanet_semiamplitude = None

        self.radius = None  # stellar radius in stellar unit
        self.area = None
        self.mass = None    # planet mass (useless for now)
        self.inclination = None #orbital inclination
        self.impact_parameter = None
        self.semimajor_axis = None

        self.gamma = None # RV zero point (including systematic RV of the star)

    def put_reference_epoch(self, Tref):
        self.reference_epoch = Tref

    def put_gamma(self, gamma):
        self.gamma = gamma

    def put_star(self, star):
        self.star_mass = star.mass
        self.star_radius = star.radius

    def put_RVparameters(self, P, K, Tc, e, o, degree=False):

        self.period = P
        self.RV_semiamplitude = K
        self.reference_transit = Tc
        self.eccentricity = e
        # Angles must be in radians
        if degree:
            self.omega = o/180. * np.pi
        else:
            self.omega = o

        self.phase = kp.kepler_Tc2phase_Tref(self.period,
                                             self.reference_transit,
                                             self.eccentricity,
                                             self.omega)

    def put_Transitparameters(self, Rp, i, b, a, degree=True):


        # degree=True by default because inclination is always expressed in degrees
        self.radius = Rp
        self.inclination = i
        self.impact_parameter = b
        self.semimajor_axis = a

    def put_RVplanet(self, Kplanet):
        self.RVplanet_semiamplitude = Kplanet

    def get_RV_kms(self, bjd):
        return self.get_RV_ms(bjd)/1000.000

    def get_RV_ms(self, bjd):
        return kp.kepler_RV_T0P(bjd-self.reference_epoch, self.phase, self.period,
                                self.RV_semiamplitude, self.eccentricity, self.omega)

    def get_RVplanet_kms(self, bjd):
        return self.get_RVplanet_ms(bjd)/1000.000

    def get_RVplanet_ms(self, bjd):
        return (-1.)*kp.kepler_RV_T0P(bjd-self.reference_epoch, self.phase, self.period,
                                self.RVplanet_semiamplitude, self.eccentricity, self.omega)
