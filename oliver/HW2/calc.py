from calculator_functions import case1_solver
from shortcut_flash import shortcut_flash
from absorber_functions import absorber_shortcut
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

# print(shortcut_flash(51379.22, 310, antoine_coeffs, epsilon))

epsilon = Epsilon(0.99, 2)
print(absorber_shortcut(51004.2, 310, epsilon, 10, antoine_coeffs))

