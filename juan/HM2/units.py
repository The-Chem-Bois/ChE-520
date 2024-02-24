'''
Juan De Leon
ID: 20765671
'''
import operators as op
import numpy as np
import math


#Shortcut Method to find Split Fractions in the Flash Unit
def short_flash(temperature, pressure, coefficients, xi_n, key):

    key_position = key-1
    p_vaps = op.antoine(temperature, coefficients,"P")
    K_k = p_vaps/pressure 
    alpha = K_k/K_k[key_position]

    xi= alpha*xi_n/(1+(alpha-1)*xi_n)
    
    return xi

#Shortcut Method to find Split Fractions in the Absorber Unit
def short_absorber(temperature, pressure, coefficients, r, key, Ae):

    key_position = key-1
    p_vaps = op.antoine(temperature, coefficients, "P")
    K_k = p_vaps/pressure
    alpha = K_k/K_k[key_position]

    N = np.log((r-Ae)/(Ae*(r-1)))/ np.log(Ae)

    A_k = Ae/alpha
    B_kN = (1-A_k**(N+1))/(1-A_k)

    Xi_kV = B_kN**-1
    Xi_kL = 1-Xi_kV

    return Xi_kV, Xi_kL

#Full absorber method to calculate the input and output streams as welll as the number of equilibrium stages
def full_absorber(temperature, pressure, coefficients, r, key, Ae, V_kn1):

    key_position = key-1
    p_vaps = op.antoine(temperature, coefficients,"P")
    K_k = p_vaps/pressure 
    K_n = K_k[key_position]
    alpha = K_k/K_n
    
    L_0= np.zeros(4)
    L_0[0] = Ae*(np.sum(V_kn1))*K_n

    N = np.log((r-Ae)/(Ae*(r-1)))/(np.log(Ae))
    
    
    A_k = Ae/alpha
    B_kN = (1-A_k**(N+1))/(1-A_k)
    B_kN1 = (1-A_k**(N))/(1-A_k)

    Xi_kV = B_kN**-1
    Xi_kL = (B_kN1)/B_kN
    
    V_k1 = (Xi_kV*V_kn1)+(Xi_kL*L_0)
    L_kN = (1-Xi_kV)*V_kn1+(1-Xi_kL)*L_0

    return L_0, V_k1, L_kN

if __name__ == "__main__":
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

    # print(antoine_coeffs)
    print("Order: W, EA, EL, DEE")
    # print(short_flash(51379.22, 310, antoine_coeffs, 0.5, 4))
    # print(short_absorber(51004.2, 310, antoine_coeffs, 0.99, 2, 10))
    u31 = np.array([44.49, 11.55, 1257.845, 1.21])
    print(full_absorber(310, 51004.2, antoine_coeffs, 0.979, 4,1.4,u31))

