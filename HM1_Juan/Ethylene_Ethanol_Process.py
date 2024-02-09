'''
Juan De Leon
ID: 20765671
'''
import numpy as np
from Hmw1 import bubble_point, dew_point, case1_solver

names = ["methane", "ethylene", "propylene", "diethyl-ether", "ethanol", "isopropanol", "water"]

#Antoine's Coefficients for the process components
coefficients = np.array([[6.61184, 389.93, 266.0],\
                      [6.74756, 585.,255.],\
                      [6.81960, 785., 247.],\
                      [7.4021,1391.4, 273.16],\
                      [8.04494,1554.3,222.65],\
                      [6.66040, 813.055, 132.96],\
                      [8.10765,1750.286,235.]])

# convert from log10->ln
coefficients[:, 0] *= np.log(10)
coefficients[:, 1] *= np.log(10)

# convert Celsius->Kelvin
coefficients[:, 2] -= 273.15

fk = np.array([[200, 1198.77,266.71,2.421,90.79,1.8802,680.72]]) #Array of componet flowrates at input

bp_guess = 310 #Temperature (K) guess since it is above the boiling point
pressure = 51379.247 #mmHg

#Question 1
#Calculates bubble point temperature for the inlet using bubble point function
bp_temp = bubble_point(bp_guess, pressure, coefficients, fk, output = "T" )
print(bp_temp)

dw_guess = 390 #Temperature (K) guess since it is above the melting point
#Calculates dew point temperature for the inlet using dew point function
dw_temp = dew_point(dw_guess, pressure, coefficients, fk, output = "T")
print(dw_temp)

#Question 2
#Calculates the flash temperature for the given problem using case 1 method
flash_temp = case1_solver(fk, coefficients, 313, pressure, 0.5, 3, fixed="P")
print(flash_temp)