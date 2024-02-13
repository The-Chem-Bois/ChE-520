'''
Oliver Erdmann
20843970
'''
import numpy as np
from calculator_functions import get_vapor_pressure, calc_bubble_point


def Absorber(n, r,  P, T, V_n1, antoine_coeffs):
    '''
    n: index of keycomponent (including 0)
    r: fix recovery (typically 0.99)
    P: Pressure
    T: Solvent temperature
    V_n1: Bottom flow rates in the absorber of each component.
    antoine_coeffs: Includes parameter A, B, and C for all keycomponents.
    '''
    Ae = 1.4
    p_vap = get_vapor_pressure(antoine_coeffs, T)

    #calculate L0
    L0 = 1.4*V_n1 * p_vap/P
    K = p_vap/P
    alpha = K/K[n]

    # Calculate number of stages via Kremser equation
    N = (np.log((r*V_n1+L0[n] - Ae*V_n1)/(L0[n]-Ae*(1-r)*V_n1)))/np.log(Ae)

    #calculate abosrption factor for all components
    A_k = 1.4/alpha
    B_k1 = (1-A_k**(N+1))/(1-A_k)
    B_k = (1-A_k**(N))/(1-A_k)

    #calculate flows
    vk = V_n1/B_k + B_k1/B_k*L0
    lk = (1-B_k1/B_k)*L0 + (1-1/B_k)*V_n1

    x = lk/np.sum(lk)
    y = vk/np.sum(vk)

    bubble = calc_bubble_point(P, T, V_n1, antoines, 'P')

    eps_v = B_k**(-1)
    eps_l = B_k1/B_k

    return (eps_v, eps_l)


def absorber_shortcut(P, T, eps, AE, antoine_coeffs):
    '''
    Get the split fractions for the absorber unit operation

    P: Pressure in mmHg
    T: Temperature in Kelvin
    eps: Epsilon value with index position and value
    AE: some parameter which needs to be specified for absorber
    antoine_coeffs: Array of anotine coefficients for each component

    '''
    
    p_vaps = get_vapor_pressure(antoine_coeffs, T);

    K = p_vaps/P
    n = eps.position -1
    epsilon = r = eps.value

    alpha = K/K[n]

    N = np.log((r-AE)/(AE*(r-1)))/np.log(AE)

    Ak = AE/(alpha)
    Bk = (1-Ak**(N+1))/(1-Ak)

    Eps_v = Bk**(-1)
    Eps_l = 1- Eps_v

    return (Eps_v, Eps_l)


if __name__ == "__main__":

    pressure = 7500 # mmhg
    temp = 300 # kelvin
    antoines = np.array([[5.40221, 1838.675, -31.737], [4.42448, 1312.253, -32.445]])

    results = Absorber(1, .95, pressure, temp, 11, antoines)
    print(results)