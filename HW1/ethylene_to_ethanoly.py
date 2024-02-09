'''
Oliver Erdmann
20843970
'''
import numpy as np;
from calculator_functions import case1_solver, calc_bubble_point, calc_dew_point;

class Epsilon: # Using this class to create my Epsilon object with a value and key position (index).
    def __init__(self, value: float, position: int) -> None:
        self.value = value
        self.position = position


# Define Constants

Pressure = 51375 # mmHg = 68.5 bar
epsilon = Epsilon(0.5, 4) # This is split fraction

# Define our input feed and keep everything in consistent order
# water, ethyl alcohol, ethylene, diethyl-ether, methane, propylene, isopropyl-alcohol

antoine_coeffs = np.array([
    [8.10765, 1750.286, 235.0],
    [8.04494, 1554.3, 222.65 ],
    [6.74756, 585, 255],
    [7.4021, 1391.4, 273.16],
    [6.61184, 389.93, 266 ],
    [6.81960, 785, 247],

    [6.66040, 813.055, 132.93]
]);

antoine_coeffs[:, :2] *= np.log(10) # Convert A, B from log base 10 to log base e.
antoine_coeffs[:, 2] -= 273.15 # Convert C coefficient from C to K.

#flow rates
u2 = np.array([680.72, 90.79, 1198.77, 2.421, 200, 266.71, 1.8802 ]);

# Execute case 1 solver
flash_results = case1_solver(epsilon, antoine_coeffs, 313.15, Pressure, u2, specification='P', tolerance=0.06);
print(flash_results)
#Execute bubble point solver
calc_bubble_point(Pressure, 310, u2, antoine_coeffs, specification='P')
# #Execute dew point solver
calc_dew_point(Pressure, 390, u2, antoine_coeffs, specification='P', tol=0.1 )