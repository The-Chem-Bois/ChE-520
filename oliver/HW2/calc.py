from calculator_functions import case1_solver, calc_bubble_point, calc_dew_point
from shortcut_flash import shortcut_flash
from absorber_functions import absorber_shortcut, Absorber, absorber_bubble_point
from distialltion_functions import distillation_column
import numpy as np
# Reactor

# EL
eta = 0.07
K = 0.2
muEL1 = 96
muEL2 = (1 - eta) * muEL1

# W
muW1 = 0.6 * muEL1
muW2 = muW1 - eta * muEL1

# EA
muEA1 = 0
muEA2 = eta * muEL1 + muEA1

# DEE
muDEE2 = (muEA2) ** 2 / (muW2 * K)

# print("EL", muEL1, muEL2, "\nW", muW1, muW2, "\nEA", muEA1, muEA2, "\nDEE", muDEE2)


class Epsilon:
    def __init__(self, value, position) -> None:
        self.value = value
        self.position = position


epsilon = Epsilon(0.5, 4)

# water, ethyl alcohol, ethylene, diethyl-ether

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
# fk = np.array([50.88, 6.72, 89.28, 4.44])
print('W, EA, ET, DEE')
# print(1- shortcut_flash(51379.22, 310, antoine_coeffs, epsilon))
# print(shortcut_flash(51379.22, 310, antoine_coeffs, epsilon))

# epsilon = Epsilon(0.99, 2)
# print(absorber_shortcut(51004.2, 310, epsilon, 10, antoine_coeffs))

epsilon = Epsilon(0.5, 4)
fk_1 = np.array([766.2, 96.29, 1277, 60.50])

case1_solver(epsilon, antoine_coeffs, 390, 51379.22, fk_1, 'P' )
calc_bubble_point(51379.22, 350, fk_1, antoine_coeffs, 'P')

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
u6 = np.array([852.27, 84.74, 38, 2.132])

calc_bubble_point(750*68, 350, u6, antoine_coeffs, 'P' )

#first distillation column (de-watering)
eps_lk = Epsilon(0.995, 1)
eps_hk = Epsilon(0.9, 0)
distillation_column(u6, 528, 68*750, eps_lk, eps_hk, antoine_coeffs, np.array([0.1, 0.995, 1, 1]), 17.56 * 750, 18.06 * 750)