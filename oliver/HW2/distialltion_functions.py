import numpy as np
from calculator_functions import get_vapor_pressure, calc_bubble_point, calc_dew_point

def distillation_column(fk, T_in, P_in, eps_lk, eps_hk, antoine_coeffs, eps, P_cond, P_reb):
    '''
    Function to determine number of stages, and reboiler and condensor temperatures.

    fk: Inlet flowrates of each component: ARRAY
    T_in: Incoming temperature of the stream in Kelvin: FLOAT
    P_in: Pressure of the inlet stream in mmHg: FLOAT
    eps_lk: Object containing split fraction (value) and index (position) of light key: EPSILON Class
    eps_hk: Object containing split fraction (value) and index (postion) of the heavy key: EPSILON CLASS
    antoine_coeffs: Antoine Coefficients of each key component: ARRAY
    eps: Desired split fractions in the vapor phase for each key component: ARRAY
    P_cond: Pressure of total condensor in mmHg: FLOAT
    P_reb: Pressure of reboiler in mmHg: FLOAT
    --- RETURNS ---

    '''
    # get indicies of heavy key and light key
    n_hk = eps_hk.position
    n_lk = eps_lk.position

    # get split fraction values of heavy key and light key
    val_hk = eps_hk.value
    val_lk = eps_lk.value

    #Get relative volatility of LK in reference to HK
    p_vaps = get_vapor_pressure(antoine_coeffs, T_in)
    
    alpha_LK_HK = (p_vaps[n_lk])/(p_vaps[n_hk])
    # Calculate number of stages: N
    N = np.log((val_lk*val_hk)/((1-val_hk)*(1-val_lk)))/np.log(alpha_LK_HK)


    tops = eps*fk
    bottoms = (1-eps)*fk

    #Calculate dew point for bottoms, and bubble point for tops
    calc_bubble_point(P_cond, T_in, tops, antoine_coeffs, 'P', 0.01, 200 )
    calc_dew_point(P_reb, T_in, bottoms, antoine_coeffs, 'P')

    print(f'{N} number of stages')

    return (N)