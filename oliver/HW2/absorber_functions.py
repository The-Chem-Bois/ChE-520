'''
Oliver Erdmann
20843970
'''
import numpy as np
from calculator_functions import get_vapor_pressure, calc_bubble_point
import math


def Absorber(n, r,  P, T, V_n1, antoine_coeffs, Ae):
    '''
    n: index of keycomponent (including 0)
    r: fix recovery (typically 0.99)
    P: Pressure
    T: Solvent temperature
    V_n1: Bottom flow rates in the absorber of each component.
    antoine_coeffs: Includes parameter A, B, and C for all keycomponents.
    '''
    p_vap = get_vapor_pressure(antoine_coeffs, T)

    #calculate L0
    K = p_vap/P
    alpha = K/K[n]
    L0 = np.zeros(4)
    L0[0] = Ae * np.sum(V_n1) * K[n]

    # Calculate number of stages via Kremser equation
    N = np.log((r-Ae)/(Ae*(r-1)))/np.log(Ae)
    N = math.ceil(N)

    #calculate abosrption factor for all components
    A_k = Ae/alpha
    B_kN = (1-A_k**(N+1))/(1-A_k)
    B_kN1 = (1-A_k**(N))/(1-A_k)

    # B_k = (1-A_k**(N))/(1-A_k)

    calc_bubble_point(P, T, V_n1, antoine_coeffs, 'P')

    eps_v = B_kN**(-1)
    eps_l = (B_kN1)/B_kN

    V1 = eps_v*V_n1 + eps_l*L0
    LN = (1-eps_v)*V_n1 + (1-eps_l)*L0


    # breakpoint()

    return (L0, V1, LN, N, eps_v, eps_l)

def absorber_shortcut(P, T, eps, AE, antoine_coeffs):
    '''
    Get the split fractions for the absorber unit operation

    P: Pressure in mmHg
    T: Temperature in Kelvin
    eps: Epsilon value with index position and value
    AE: some parameter which needs to be specified for absorber
    antoine_coeffs: Array of anotine coefficients for each component

    '''
    #Get vapor pressures using Antoine's equation for each component
    p_vaps = get_vapor_pressure(antoine_coeffs, T);
    #Calculate K, specify postion, and specify r
    K = p_vaps/P
    n = eps.position -1
    epsilon = r = eps.value

    #Calculate relative volatility
    alpha = K/K[n]
    # Calculate N number of stages via Kremser Equation
    N = np.log((r-AE)/(AE*(r-1)))/np.log(AE)

    Ak = AE/(alpha)
    Bk = (1-Ak**(N+1))/(1-Ak)
    # Find the split fractions in the vapor and liquid phase
    Eps_v = Bk**(-1)
    Eps_l = 1- Eps_v

    return (Eps_v, Eps_l)