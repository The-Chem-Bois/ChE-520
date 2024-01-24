#CHE 520 Homework 1 Ethylene to Ethanol Process
#Christina Mohan
#20844467

#Import numpy to access required functions
import numpy as np
#Import the functions from the other file
from HW1_Functions_CM import *

#Define the Antoine Equation parameters from Table 1.3
names = ["methane", "ethylene", "propylene", "diethyl-ether", "ethanol", "isopropanol", "water"]
coefficients = np.array([[6.61184, 389.93, 266.0],\
                      [6.74756, 585.,255.],\
                      [6.81960, 785., 247.],\
                      [7.4021,1391.4, 273.16],\
                      [8.04494,1554.3,222.65],\
                      [6.66040, 813.055, 132.96],\
                      [8.10765,1750.286,235.]])

# convert from log10->ln
coefficients[:, 0] *= np.log(10.)
coefficients[:, 1] *= np.log(10.)

# convert Celsius->Kelvin
coefficients[:, 2] -= 273.15

#Create arrays for each Antoine parameter
A = coefficients[:,0]
B = coefficients[:,1]
C = coefficients[:,2]

#Define the molar flowrates for each component in gmol/s
mu = np.array([200, 1198.77,266.71,2.421,90.79,1.8802,680.72])

n = 3 #Let DEE as the key component 

P = 68.5*750.062 #Pressure in mmHg
T_guess_bubble = 310 #Initial temperature guess in K; set to be slightly above the boiling point of DEE to ensure it vaporizes
T_guess_dew = 390 #Initial temperature guess in K; set to be slightly above the melting point of DEE to ensure it vaporizes

T_bubble = bubble_point(P, T_guess_bubble, mu, A, B, C, n, Find_T = True) #Calculate the bubble point

#T_dew = dew_point(P, T_guess_dew, mu, A, B, C, n, Find_T = True) #Calculate the dew point

eps = 0.5

#c1_flash = case1_flash(eps, P, T_guess, mu, A, B, C, n)
#print(c1_flash)

