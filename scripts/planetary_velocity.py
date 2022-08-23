"""from classes.kepler_exo import *

# Mass of the star HD189733 (in Solar masses)
#Ms = 0.823
Ms = 1.148
# Mass of the planet (in Solar masses)
#Mp = 1.138 / 1.047348644e3
Mp = 0.69 / 1.047348644e3

K1 = kepler_K1(Mp,Ms,3.52474854657,86.59,0.0082)
print K1


## update
"""

import matplotlib.pyplot as plt
import numpy as np
from SLOPpy.subroutines.constants import *
import argparse
from scipy.optimize import fsolve
from SLOPpy.subroutines.kepler_exo import *

def get_mass(M_star2, M_star1, Period, K1, e0):
    # M_star1, M_star2 in solar masses
    # P in days -> Period is converted in seconds in the routine
    # inclination assumed to be 90 degrees
    # Gravitational constant is given in m^3 kg^-1 s^-2
    # output in m/s
    output = K1 - (2. * np.pi * G_grav * M_sun / 86400.0) ** (1.0 / 3.0) * (1.000 / np.sqrt(1.0 - e0 ** 2.0)) * (
                                                                                                                    Period) ** (
                                                                                                                    -1.0 / 3.0) * (
                      M_star2 * (M_star1 + M_star2) ** (-2.0 / 3.0))
    return output


parser = argparse.ArgumentParser(prog='planetary_velocity.py', description='Compute the expected semi-amplitude of the planet')
parser.add_argument('star_mass', type=float, nargs=1, help='Stellar mass [solar]')
parser.add_argument('plan_mass', type=float, nargs=1, help='Planet mass [Jupiter/Earth/K units, default Jupiter]')
parser.add_argument('period', type=float, nargs=1, help='Planetary period [days')
parser.add_argument('inclination', type=float, nargs=1, help='Planetary inclination [degrees]')
parser.add_argument('eccentricity', type=float, nargs=1, help='Planetary eccentricity [pure]')
parser.add_argument('-e', type=float, nargs='?', default=False, const=True, help='Planetary mass in Earth units')
parser.add_argument('-k', type=float, nargs='?', default=False, const=True, help='Planetary mass in m/s')
args = parser.parse_args()


star_mass = args.star_mass[0]
P = args.period[0]
i = args.inclination[0]
e = args.eccentricity[0]
if args.e and args.k:
    print('Either -k or -e, not both!')
    quit()

planet_mass = args.plan_mass[0]
if args.e:
    planet_mass *= Mears
elif args.k:
    x0 = Mjups
    K_input = planet_mass
    planet_mass = fsolve(get_mass, x0, args=(star_mass, P, K_input, e))
else:
    planet_mass *= Mjups

print(kepler_K1(planet_mass, star_mass, P, i, e))


#sampler = args.sample[0]
#file_conf = args.config_file[0]
