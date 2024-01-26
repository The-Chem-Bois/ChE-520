from numpy import *

from Hmw1 import antoine

names = ["methane", "ethylene", "propylene", "diethyl-ether", "ethanol", "isopropanol", "water"]
coefficients = array([[6.61184, 389.93, 266.0],\
                      [6.74756, 585.,255.],\
                      [6.81960, 785., 247.],\
                      [7.4021,1391.4, 273.16],\
                      [8.04494,1554.3,222.65],\
                      [6.66040, 813.055, 132.96],\
                      [8.10765,1750.286,235.]])

# convert from log10->ln
coefficients[:, 0] *= log(10.)
coefficients[:, 1] *= log(10.)

# convert Celsius->Kelvin
coefficients[:, 2] -= 273.15

# convert from mmHG to bar
coefficients[:, 0:2] -= log(750.0638)

print(antoine(1., coefficients, output="T"))
print(antoine(310., coefficients, output="P"))
