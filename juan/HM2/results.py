'''
Juan De Leon
ID: 20765671
'''
import numpy as np
import units as un
import operators as op

#Antoine Coefficients for the process elements
antoine_coeffs = np.array(
    [[8.10765, 1750.286, 235.0], # water W
    [8.04494, 1554.3, 222.65], # ethanol EA
    [6.74756, 585, 255], # ethylene EL
    [7.4021, 1391.4, 273.16],] # Diethyl ether DEE
)
antoine_coeffs[:, :2] *= np.log(10) # Convert A, B from log base 10 to log base e.
antoine_coeffs[:, 2] -= 273.15 # Convert C coefficient from C to K.

print("Order: W, EA, EL, DEE")
 
#Get the split fractions from the flash unit
#Using 68.5 bar, 310K and 0.5 split fraction DEE
#Vapor Stream u31 (top) split fractions
print("Flash Split Fractions")
print(un.short_flash(310, 51379.22, antoine_coeffs, 0.5, 4))

#Liquid Stream u32 (bottom) split fractions
print(1-un.short_flash(310, 51379.22, antoine_coeffs, 0.5, 4))

#Get split fractions from the Absorber unit
#Using 68 bar, 310K, 0.99 recovery of Ethanol and AE = 10
print("Absorber unit split fractions")
print(un.short_absorber(310, 51004.2, antoine_coeffs, 0.99, 2, 10))

#Use Case1 Solver to find the Operating Temperature
#Uses u2 stream
print("Flash Unit Calculations")
u2 = np.array([[726.86, 96.6, 1275.43, 2.57]])
print(op.case1_solver(u2, antoine_coeffs, 310, 51379.22, 0.5, 4, "P"))

#Calculates the flash bubble point using the vapor stream u31
u31 = np.array([[39.25, 11.79, 1256.3, 1.29]])
print(op.bubble_point(310, 51379.22, antoine_coeffs, u31, "T" ))

#Calculates the flash dew point using the liquid stream u32
u32 = np.array([[687.61, 84.81, 19.13, 1.29]])
print(op.dew_point(310, 51379.22, antoine_coeffs, u32, "T"))

#Calculates the absorber unknowns
#Using 310K, 68 bar, 0.99 recovery of Diethyl ether, AE = 1.4, and the composition from u31
print("Absorber Unit Calculations")
u03, u41, u42 = un.full_absorber(310, 51004.2, antoine_coeffs, 0.99, 4, 1.4, u31)
print(u03, u41, u42)

#Calculates the absorber bubble point using the liquid stream u42
print(op.bubble_point(310, 51379.22, antoine_coeffs, u42, "T" ))

#Calculates the absorber dew point using the vapor stream u41
print(op.dew_point(310, 51379.22, antoine_coeffs, u41, "T" ))

#Calculates the bubble point of the mixer (unit 6)
print("Mixer Calculations")
u6 = np.array([[755.12, 96.48, 45.51, 2.27]])
print(op.bubble_point(310, 51004.2, antoine_coeffs, u6, "T" ))

#Calculates the de-watering column bubble point, and dewpoints using the vapor stream u71 and liquid stream u72
print("Dewatering Column Calculations")
u71 = np.array([[75.51, 96, 45.51, 2.27]])
u72 = np.array([[679.61, 0.48, 0, 0]])
print(op.bubble_point(310, 13171.081, antoine_coeffs, u71, "T" ))
print(op.dew_point(521.01, 13546.112, antoine_coeffs, u72, "T" ))

#Calculates the de-ethering column dewpoint using liquid stream u82
print("De-ethering Column Calculations")
u82 = np.array([[75.51, 95.52, 0, 0.0023]])
print(op.dew_point(318.3, 8400.69, antoine_coeffs, u82, "T" ))

#Calculates the finishing column bubble point, and dewpoints using the vapor stream u91 and liquid stream u92
print("Dewatering Column Calculations")
u91 = np.array([[13.78, 95.04, 0, 0.0023]])
u92 = np.array([[61.73, 0.48, 0, 0]])
print(op.bubble_point(310, 750.062, antoine_coeffs, u91, "T" ))
print(op.dew_point(310, 1125.09, antoine_coeffs, u92, "T" ))