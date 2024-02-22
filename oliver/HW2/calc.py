from calculator_functions import case1_solver, calc_bubble_point, calc_dew_point
from shortcut_flash import shortcut_flash
from absorber_functions import absorber_shortcut, Absorber, absorber_bubble_point
from distialltion_functions import distillation_column
import numpy as np


class Epsilon:
    def __init__(self, value, position) -> None:
        self.value = value
        self.position = position



# water, ethyl alcohol, ethylene, diethyl-ether (order)
# Define antoine constants
antoine_coeffs = np.array(
    [
        [8.10765, 1750.286, 235.0], # water
        [8.04494, 1554.3, 222.65], # ethanol
        [6.74756, 585, 255], # ethylene
        [7.4021, 1391.4, 273.16], # DEE
    ]
)
antoine_coeffs[:, :2] *= np.log(10) # Convert A, B from log base 10 to log base e.
antoine_coeffs[:, 2] -= 273.15 # Convert C coefficient from C to K.

print('W, EA, ET, DEE')
# Define split fraction for diethyl-ether
epsilon = Epsilon(0.5, 4)

#Print split fractions for top and bottom streams exiting flash.
print(1- shortcut_flash(51379.22, 310, antoine_coeffs, epsilon))
print(shortcut_flash(51379.22, 310, antoine_coeffs, epsilon))

#Define split fraction of ethylene for absorber
epsilon = Epsilon(0.99, 2)
print(absorber_shortcut(51004.2, 310, epsilon, 10, antoine_coeffs))


# Now do flash calculations and solve for temperatures in each unit operation

epsilon = Epsilon(0.5, 4)
fk_1 = np.array([766.2, 96.29, 1277, 60.50])
print ('---- flash calculations ------')
case1_solver(epsilon, antoine_coeffs, 390, 51379.22, fk_1, 'P' )
u31 = np.array([44.49, 11.55, 1257.845,1.21 ])
u32 = np.array([779.38, 84.74, 19.115, 1.21])
calc_bubble_point(51379.22, 350, u31, antoine_coeffs, 'P')
calc_dew_point(51379.22, 350, u32, antoine_coeffs, 'P', 5, 300)

print('--- Absorber Calculations-------')

#calculate absorber
u31 = np.array([44.49, 11.55, 1257.845, 1.21])
absorber_results = Absorber(3, 0.979, 68*750, 310, u31, antoine_coeffs, 1.4 )

u42 = absorber_results[2]
absorber_bubble_point(68*750, 340, antoine_coeffs, absorber_results[2], 0, 3)
calc_dew_point(68*750, 340, absorber_results[1], antoine_coeffs, 'P')
print('L0 flows of absorber', absorber_results[0])
print('mu41 of absorber (vap out)', absorber_results[1])
print('mu42 of absorber (liq out)', absorber_results[2])

#calculate mixer exit temperature
print ('--- mixer calculations ---')
u6 = np.array([852.27, 84.74, 38, 2.132])

calc_bubble_point(750*68, 350, u6, antoine_coeffs, 'P' )

#first distillation column (de-watering)
print('---first distillation column calculations---')
eps_lk = Epsilon(0.995, 1)
eps_hk = Epsilon(0.9, 0)
distillation_column(u6, 528, 68*750, eps_lk, eps_hk, antoine_coeffs, np.array([0.1, 0.995, 1, 1]), 17.56 * 750, 18.06 * 750)

#Second distillation column (de-ethering)
print('---second distillation column calculations---')
eps_lk = Epsilon(0.995, 3)
eps_hk = Epsilon(0.995, 1)
u71 = np.array([85.23, 38, 95.69, 2.132])

distillation_column(u71, 327.40, 17.56*750, eps_lk, eps_hk, antoine_coeffs, np.array([0, 0.005, 1, 0.995]), 10.7 * 750, 11.2 * 750 )

#Finishing distillation column
print('---final distillation column calculations---')

eps_lk = Epsilon(0.995, 1)
eps_hk = Epsilon(0.811, 0)
u82 = np.array([85.23, 95.21, 0, 0.011])

distillation_column(u82, 449, 17*11.2, eps_lk, eps_hk, antoine_coeffs, np.array([0.189, 0.995, 0, 1]), 1*750, 1.5*750)