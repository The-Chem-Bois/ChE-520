'''
Oliver Erdmann
20843970
Functions written in Python3
'''
import numpy as np

def get_vapor_pressure(antoine_coeffs, T):
    # calculates vapor pressure of all components provided array of antoine coefficients and temperature.
    P_vaps = np.zeros(len(antoine_coeffs))

    for i in range(0, len(antoine_coeffs)):
        # Get antoine coefficients and calculate vapor pressures
        A, B, C = antoine_coeffs[i]
        P_vap = np.exp(A - B/(C+T))
        P_vaps[i] = P_vap
    return P_vaps

def calc_relative_volatility(epsilon_or_phi, P, fk, P_vaps, parameter: str):
    # Calculates epsilon values, relative volatilities, vapor and liquid compositions, vapor flow, liquid flow, total flow, and average volatility
    # Parameter can be 'epsilon' or 'phi' as strings.

    K_n = P_vaps/P # k values (P vapor pressure)/ (Total Pressure) for all key components


    eps_n = None # initialize variable eps_n
    alpha_k = None # initialize variable alpha_k for volatility of each component
    
    if parameter == 'epsilon':
        alpha_k = K_n / K_n[epsilon_or_phi.position - 1] # calculates array of relative volatilities for each component relative to the epsilon component provided
        eps_n = (alpha_k * epsilon_or_phi.value)/(1 + (alpha_k - 1)* epsilon_or_phi.value) # calculate epsilon for each component
    elif parameter == 'phi':
        theta = K_n * epsilon_or_phi / (1- epsilon_or_phi)
        eps_n = theta/(1 + theta)
        # By default, select the first key
        alpha_k = K_n/ K_n[0]

    # Calculate Vapor and Liquid flow rates (total and for each component)
    V_k = eps_n*fk
    L_k = (1-eps_n)*fk
    V = sum(V_k)
    L = sum(L_k)
    Total_Flow = V + L
    # Calculate compositions in liquid (x) and vapor (y) phases.
    x_k = L_k/L
    y_k = V_k/V

    alpha_bar = sum(x_k * alpha_k) # calculate average volatility

    return (K_n, alpha_k, eps_n, V_k, L_k, V, L, Total_Flow, x_k, y_k, alpha_bar)

def calc_bubble_dew_T(alpha_k, alpha_bar, antoine_coeffs, index ): 
    #Calcuate vapor pressure based in terms of alpha calculated on majority liquid key component
    P_k_vap = alpha_k[index]/alpha_bar * P
    # Rearrange antoine's equation and solve for T
    A,B,C = antoine_coeffs[index]
    T_calc = B/(A - np.log(P_k_vap)) - C

    return T_calc

def calc_bubble_dew_P(alpha_k, alpha_bar, P_vaps, index):
     
     P_calc = alpha_bar/alpha_k[index] * P_vaps[index]

     return P_calc


def case1_solver(epsilon, antoine_coeffs, T: float, P: float, fk, specification: str, tolerance = 0.5, maxiter= 50):
    '''
    Use function when epsilon is defined, and if T or P is defined.

    epsilon: An object with property value for epsilon as well as its key position. (i.e. epsilon for component 2 would therefore be 2).
    antoine_coeffs: Array containing arrays of each key's anotine coefficients. (Order matters, requires 3 coefficients A, B, C)
    T: The initial guessed or given temperature in Kelvin.
    P: The initial guessed or given pressure in mmHg.
    fk: Array of the fk value for each component in mol/s. (Order matters).
    specification: String of 'P' for pressure, or 'T' for temperature that indicated which parameter is specified.
    tolerance: Acceptable tolerance, default is 0.5.
    maxiter: The maximum number of iterations, default is 50.
    '''

    iterations = 1

    while iterations < maxiter:

        #Get an array of vapor pressures for each key component
        P_vaps = get_vapor_pressure(antoine_coeffs, T)
        # Calculate volatility, epsilon, and liquid composition values
        K_n, alpha_k, eps_n, V_k, L_k, V, L, Total_Flow, x_k, y_k, alpha_bar = calc_relative_volatility(epsilon, P, fk, P_vaps, 'epsilon')

        majority_index = np.argmax(x_k) # get the index of the component that has the highest liquid composition


        if specification == 'P':
            # Running this code block if T was guessed
            T_calc = calc_bubble_dew_T(alpha_k, alpha_bar, antoine_coeffs, majority_index)

            if (np.abs(T_calc - T) <= tolerance):
                print(f'Converged after {iterations} iterations! Epsilons of each component: {eps_n}, temperature: {T} Kelvin')
                break;
            else:
                T = (T + T_calc)/ 2
                iterations += 1

        elif specification == 'T':
            # Running this code block if P was guessed
            P_calc = calc_bubble_dew_P(alpha_k, alpha_bar, P_vaps, majority_index)
        
            if (np.abs(P_calc - P) <= tolerance):
                print(f'Converged after {iterations} iterations! Epsilons of each component: {eps_n}, pressure: {P} mmHg')
                break;
            else:
                P = (P + P_calc)/2
                iterations += 1

        if iterations >= maxiter:
            print(f'Failed to converge')
            break;

def case2_solver(epsilon, antoine_coeffs, T: float, P: float, fk, tolerance = 0.5, maxiter = 50):
    '''
    For a specified T and P, choose a key component n and specify its epsilon

    epsilon: An object with a guessed property value for epsilon as well as its key position (n). (properties: value, position)
    antoine_coeffs: Array containing arrays of each key's anotine coefficients. (Order matters, requires 3 coefficients A, B, C)
    T: Temperature in Kelvin.
    P: Pressure in mmHg.
    fk: Array of the fk value for each component in mol/s. (Order matters).
    tolerance: Acceptable tolerance, default is 0.5.
    maxiter: The maximum number of iterations, default is 50.
    '''

    P_vaps = np.zeros(len(fk)) # create an array for storing all key component vapor pressures
    K_n = np.zeros(len(fk)) # create an array for storing all k values (P vapor pressure)/ (Total Pressure) for all key components
    iterations = 1

    while iterations < maxiter:

        #Get an array of vapor pressures for each key component
        P_vaps = get_vapor_pressure(antoine_coeffs, T)

        # Calculate K, epsilon, and flow rates (vapor and total)
        K_n, alpha_k, eps_n, V_k, L_k, V, L, Total_Flow, x_k, y_k, alpha_bar = calc_relative_volatility(epsilon, P, fk, P_vaps)

        phi = V/Total_Flow
        

        theta = K_n[epsilon.position - 1]*phi/(1-phi)
        epsilon_prime = theta/(1+theta)

        if (np.abs(epsilon.value - epsilon_prime) <= tolerance):
            print(f"Converged after {iterations} iterations. The following epsilons for each component are: {eps_n}")
            break;
        else:
            epsilon.value = (epsilon.value + epsilon_prime)/2
            iterations += 1
        if iterations >= maxiter:
            print("failed to converge")
            break;
    
def case3_solver(phi, antoine_coeffs, T: float, P: float, fk, specification: str, tolerance = 0.01, maxiter = 50):
    '''
    For a specified T and P, choose a key component n and specify its epsilon

    phi: An object with a guessed property value for epsilon as well as its key position (n). (properties: value, position)
    antoine_coeffs: Array containing arrays of each key's anotine coefficients. (Order matters, requires 3 coefficients A, B, C)
    T: Temperature in Kelvin.
    P: Pressure in mmHg.
    fk: Array of the fk value for each component in mol/s. (Order matters).
    tolerance: Acceptable tolerance, default is 0.5.
    maxiter: The maximum number of iterations, default is 50.
    '''
    iterations = 1

    while iterations < maxiter:

        P_vaps = get_vapor_pressure(antoine_coeffs, T);
        K_n, alpha_k, eps_n, V_k, L_k, V, L, Total_Flow, x_k, y_k, alpha_bar = calc_relative_volatility(phi, P, fk, P_vaps, parameter='phi' )

        majority_index = np.argmax(x_k) # get the index of the component that has the highest liquid composition

        if specification == 'P':
            # Running this code block if T was guessed
            T_calc = calc_bubble_dew_T(alpha_k, alpha_bar, antoine_coeffs, majority_index)

            if (np.abs(T_calc - T) <= tolerance):
                print(f'Converged after {iterations} iterations! Epsilons of each component: {eps_n}, temperature: {T} Kelvin')
                break;
            else:
                T = (T + T_calc)/ 2
                iterations += 1

        elif specification == 'T':
            # Running this code block if P was guessed
            P_calc = calc_bubble_dew_P(alpha_k, alpha_bar, P_vaps, majority_index)
        
            if (np.abs(P_calc - P) <= tolerance):
                print(f'Converged after {iterations} iterations! Epsilons of each component: {eps_n}, pressure: {P} mmHg')
                break;
            else:
                P = (P + P_calc)/2
                iterations += 1

        if iterations >= maxiter:
            print(f'Failed to converge')
            break;
        

    


if __name__ == "__main__":

    class Epsilon:
        def __init__(self, value: float, position: int) -> None:
            self.value = value
            self.position = position
    
    #case 1 we are providing Eps_2 = 0.8 and P = 1 bar
    eps1 = Epsilon(0.9, 2 )
    P = 750 ## 1 bar = 750 mmHg
    T = 390 ## Initial guess
    fk = np.array([30, 50, 40]);
    phi = 0.5
    antoine_coeffs = np.array([[15.9008, 2788.51, -52.34], [16.0137, 3096.52, -53.67], [16.1156, 3395.57, -59.44]]);

    # case1_solver(eps1, antoine_coeffs, T, P, fk, 'P', 0.01);
    # case2_solver(eps1, antoine_coeffs, 385, 750, fk, tolerance=0.001, maxiter=100)\
    case3_solver(phi, antoine_coeffs, 390, P, fk, 'P');

''' 
Questions to ask Nasser

1. Do units of fk matter? Or do we only care that they are consistent with each other?
2. Extra examples so I could verify my functions?
3. What is fk? What is K?
4. Bubble point and dew point equations explanation
5. What is V_k and L_k?
6. !!! For case3, when calculating alpha_k_n (volatility) for each component. W.r.t. to which component am I calculating each volatility? (Choose)
7. What are we returning for "flash calculations" ? 

'''