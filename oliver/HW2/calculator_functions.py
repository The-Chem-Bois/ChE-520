'''
Oliver Erdmann
20843970
Functions written in Python3
'''
import numpy as np

def get_vapor_pressure(antoine_coeffs, T):
    '''
    Gets vapor pressures for any number of key components based on the length of antoine coefficients at a temperature T.

    antoine_coeffs: array of A,B,C antoine coefficients for each key component. Supports units of Kelvin and natural log calculations.
    T: Temperature of system in Kelvin.
    --Returns--
    vapor pressures
    '''
    # Initialize an array with enough space to store the vapor pressure of each component.
    P_vaps = np.zeros(len(antoine_coeffs))

    for i in range(0, len(antoine_coeffs)):
        # Get antoine coefficients and calculate vapor pressures
        A, B, C = antoine_coeffs[i]
        P_vap = np.exp(A - B/(C+T))
        # Assign vapor pressure to spot in array
        P_vaps[i] = P_vap
    return P_vaps

def calc_relative_volatility(epsilon_or_phi, P, fk, P_vaps, parameter: str, key_component = 0):
    '''
    Calculates and returns the epsilon values, relative volatilties, vapor and liquid compositions, vapor flow, liquid flow, total flow, and relative volatilities
    and average volatilities.

    epsilon_or_phi: Accepts a split fraction value (epsilon) or Phi value. Float.
    P: Pressure of system in mmHg. Float.
    fk: Flow rates in feed of each key component. Make sure units are consistent. Type is Array.
    P_vaps: Vapor pressure of each key component in mmHg. Type is array.
    parameter: A string specifying 'epsilon' for split fraction or 'phi' if phi was specified.
    key_component: Specify index of majority key. If no key specified will default to 0, which may sometimes work.

    --- Returns ---
    In order:
    K values, relative volatilities, split fractions, vapor flows, liquid flow, total vapor flow, total liquid flow, total flow, 
    liquid compositions, vapor compositions, average relative volatility.
    '''
    

    K_n = P_vaps/P # k values (P vapor pressure)/ (Total Pressure) for all key components.


    eps_n = None # initialize variable eps_n
    alpha_k = None # initialize variable alpha_k for volatility of each component.
    
    if parameter == 'epsilon':
        alpha_k = K_n / K_n[epsilon_or_phi.position - 1] # calculates array of relative volatilities for each component relative to the epsilon component provided
        eps_n = (alpha_k * epsilon_or_phi.value)/(1 + (alpha_k - 1)* epsilon_or_phi.value) # calculate epsilon for each component
    elif parameter == 'phi':
        theta = K_n * epsilon_or_phi / (1- epsilon_or_phi)
        eps_n = theta/(1 + theta)
        # By default, select the first key
        alpha_k = K_n/ K_n[key_component]

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

def check_T(P, alpha_k, alpha_bar, antoine_coeffs, index ):
    '''
    Calculates the temperature based on relative volatilities and pressure. Reverses antoine's equation to get temperature.

    P: Pressure in mmHg
    alpha_k: Relative volatilities of each key component as an array
    alpha_bar: Average relative volatitiles as float
    antoine_coeffs: Array of A, B, and C antoine coefficient values. Support K and natural log.
    index: The index of the majority key component in the liquid phase. Integer.
    --returns--
    temperature in Kelvin
    '''
    #Calcuate vapor pressure based in terms of alpha calculated on majority liquid key component
    P_k_vap = alpha_k[index]/alpha_bar * P
    # Rearrange antoine's equation and solve for T
    A,B,C = antoine_coeffs[index]
    T_calc = B/(A - np.log(P_k_vap)) - C

    return T_calc

def check_P(alpha_k, alpha_bar, P_vaps, index):
     '''
     Function calcualtes the pressure based on relative volatility and vapor pressure.

     alpha_k: Relative volatitiles of each component as an array.
     alpha_bar: Average relative volatility as float.
     P_vaps: The vapor pressure of each key component in mmHg as an array.
     index: The index of the majority key component in the liquid phase. Integer.
     --returns--
     Pressure in mmHg
     '''
     
     P_calc = alpha_bar/alpha_k[index] * P_vaps[index]

     return P_calc

def calc_bubble_point (P, T, fk, antoines, specification, tol = 0.01, maxiter = 50):
    '''
    Calculates either the Pressure or Temperature in the extreme case of a bubble point. Prints result, no return.

    P: Provide Pressure or guessed pressure in mmHg
    T: Provide Temperature or guessed temperature in K
    fk: Array of feed flowrate for each component
    specification: A string value of 'T' or 'P' indicating whether the specified parameter is temperature or pressure respectively.
    antoines: Array of antoine constants (A, B, C) for each component
    zk: Array of mole fractions of each component in feed
    key_component: Provide integer of most abundant component in feed.
    tolerance: default is 0.01 - float
    maxiter: maximum number of iteratiions, default is 50 - interger
    '''
    iterations = 1

    while iterations < maxiter:
        epsilon = 0 # for bubble point flash calcualtion, epsilon (split fraction) will be 0.
        lk = fk # liquid flowrate of component is equal to feed flowrate of component
        zk = fk/sum(fk)
        xk = zk # liquid mole fraction of component is equal to mole fraction of component in feed.

        n = np.argmax(fk) # sets the index of the most abumdant component in the feed.

        P_vaps = get_vapor_pressure(antoines, T)
        K_k = P_vaps/P #get K values for all components and use it to calculate their vapor pressures

        alpha_k = K_k/K_k[n] # get relative volaities for each component relative to most abundant one in feed.
        alpha_bar = sum(xk * alpha_k) # calculates average volatility here

        if specification == 'P':
            # If Pressure is sepcifed this block runs. Calculates temperature and compares to guessed T.
            P_k_vap = alpha_bar**-1 * P

            A,B,C = antoines[n]
            T_calc = B/(A-np.log(P_k_vap)) - C
            iterations += 1
            # breakpoint()
            if (np.abs(T_calc - T) <= tol):
                print (f'Bubble point found at {T_calc} Kelvin after {iterations} iterations!')
                break;
            elif (iterations >= maxiter):
                print('Did not converge')
                break;

            T = (T_calc + T)/2

        elif specification == 'T':
            # If Temperature is specified, this block runs. Calculates pressure and compares to guessed P.
            P_calc = alpha_bar * P_vaps[n]
            iterations += 1

            if (np.abs(P_calc - P) <= tol):
                print(f'Bubble point pressure found at {P_calc} mmHg after {iterations} iterations!')
                break;
            elif (iterations >= maxiter):
                print('Did not converge')
                break;
            P = (P_calc + P)/2

def calc_dew_point (P, T, fk, antoines, specification, tol = 0.5, maxiter = 50 ):

    '''
    Calculates either the pressure (mmHg) or temperature (Kelvin) of an extreme dew point. Prints result, does not return.

    P: Pressure in mmHg, either specified or guessed. Float.
    T: Temperature in Kelvin, either specified or guessed. Float.
    fk: Flowrates of each component in feed. Consistent units. Array.
    antoines: Antoince coefficients of A,B,C of each component in an array. Supports kelvin and natural log calculations.
    specification: A string value of 'T' or 'P' indicating whether temperature or pressure is specified respectively as a parameter.
    tol: float of tolerance. Default is 0.5
    maxiter: maximum number of iterations. Default is 50. Integer
    '''
    iterations = 1


    while iterations < maxiter:
        # Define conditions for dew point
        epsilons = np.ones(len(fk)) # split fractions
        vk = fk
        yk = zk = fk/sum(fk)
        # get vapor pressures and index of majority key in the feed
        P_vaps = get_vapor_pressure(antoines, T)
        n = np.argmax(fk) 
        # calculate k and relative volatilities
        K_k = P_vaps/P
        alpha_k = K_k/K_k[n]

        if specification == 'P':
            #if specification is pressure, run this block. Calculates for Temperature by calculating vapor pressure, and reversing antoines equation for T.
            P_k_vap = P * sum(yk/alpha_k)
            A,B,C = antoines[n]
            T_calc = B/(A-np.log(P_k_vap)) - C
            iterations += 1
            if (np.abs(T_calc - T) <= tol):
                print (f'Dew point found at {T_calc} Kelvin after {iterations} iterations!')
                break;
            elif (iterations >= maxiter):
                print('Did not converge')
                break;
            # If did not converge but loop is still continuing, update T.
            T = (T_calc + T)/2
        
        elif specification == 'T':
            # Block runs if specification is Temperature. Calculates for Pressure and compares for convergence.
            P_calc = P_vaps[n] * (sum(yk/alpha_k)**-1)
            iterations += 1

            if (np.abs(P_calc - P) <= tol):
                print(f'Dew point found at pressure {P_calc} mmHg after {iterations} iterations!')
                break;
            elif (iterations >= maxiter):
                print('Did not converge')
                break;
            # Update P if convergence is not attained but loop continues.
            P = (P_calc + P)/2




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
    --returns--
    relative volatilities, average volatility, split fractions, Pressure, Temperature, vapor flow rate, liquid flow rate, total flow rate,
    liquid fractions, vapor fractions.
    '''

    iterations = 1

    while iterations < maxiter: # Run loop while iterations is less than maxiter that was specified

        #Get an array of vapor pressures for each key component
        P_vaps = get_vapor_pressure(antoine_coeffs, T)
        # Calculate volatility, epsilon, and liquid composition values
        K_n, alpha_k, eps_n, V_k, L_k, V, L, Total_Flow, x_k, y_k, alpha_bar = calc_relative_volatility(epsilon, P, fk, P_vaps, 'epsilon')

        majority_index = np.argmax(x_k) # get the index of the component that has the highest liquid composition


        if specification == 'P':
            # Running this code block if T was guessed
            T_calc = check_T(P, alpha_k, alpha_bar, antoine_coeffs, majority_index)

            if (np.abs(T_calc - T) <= tolerance): # check if it meets our tolerance criteria
                print(f'Converged after {iterations} iterations! Epsilons of each component: {eps_n}, temperature: {T} Kelvin')
                return(alpha_k, alpha_bar, eps_n, P, T, V, L, Total_Flow, x_k, y_k )
            else:
                # Update T and add to iteration count if did not meet tolerance criteria
                T = (T + T_calc)/ 2
                iterations += 1

        elif specification == 'T':
            # Running this code block if P was guessed
            P_calc = check_P(alpha_k, alpha_bar, P_vaps, majority_index)
        
            if (np.abs(P_calc - P) <= tolerance):
                print(f'Converged after {iterations} iterations! Epsilons of each component: {eps_n}, pressure: {P} mmHg')
                return(alpha_k, alpha_bar, eps_n, P, T, V, L, Total_Flow, x_k, y_k )
            else:
            # Update T and add to iteration count if did not meet tolerance criteria
                P = (P + P_calc)/2
                iterations += 1

        if iterations >= maxiter:
            print(f'Failed to converge')
            break;

def case2_solver(epsilon, antoine_coeffs, T: float, P: float, fk, tolerance = 0.5, maxiter = 50):
    '''
    For a specified T and P, choose a key component n and specify its split fraction.

    epsilon: An object with a guessed property value for epsilon as well as its key position (n). (properties: value, position)
    antoine_coeffs: Array containing arrays of each key's anotine coefficients. (Order matters, requires 3 coefficients A, B, C)
    T: Temperature in Kelvin.
    P: Pressure in mmHg.
    fk: Array of the fk value for each component in mol/s. (Order matters).
    tolerance: Acceptable tolerance, default is 0.5.
    maxiter: The maximum number of iterations, default is 50.
    --returns--
    relative volatitilies, average relative volatility, split fractions, Vapor flow rate, liquid flow rate, total flow rate,
    liquid fractions, vapor fractions.
    '''

    P_vaps = np.zeros(len(fk)) # create an array for storing all key component vapor pressures
    K_n = np.zeros(len(fk)) # create an array for storing all k values (P vapor pressure)/ (Total Pressure) for all key components
    iterations = 1

    while iterations < maxiter:

        #Get an array of vapor pressures for each key component
        P_vaps = get_vapor_pressure(antoine_coeffs, T)

        # Calculate K, epsilon, and flow rates (vapor and total)
        K_n, alpha_k, eps_n, V_k, L_k, V, L, Total_Flow, x_k, y_k, alpha_bar = calc_relative_volatility(epsilon, P, fk, P_vaps, 'epsilon')
        # Calculate phi, theta, and epsilon
        phi = V/Total_Flow
        

        theta = K_n[epsilon.position - 1]*phi/(1-phi)
        epsilon_prime = theta/(1+theta)
        # Check if the epsilon difference meet the specified tolerance
        if (np.abs(epsilon.value - epsilon_prime) <= tolerance):
            print(f"Converged after {iterations} iterations. The following epsilons for each component are: {eps_n}")
            return (alpha_k, alpha_bar, eps_n, V, L, Total_Flow, x_k, y_k)
        else:
            #If did not converge, update epsilon object and add iterations
            epsilon.value = (epsilon.value + epsilon_prime)/2
            iterations += 1
        if iterations >= maxiter:
            print("failed to converge")
            break;
    
def case3_solver(phi, antoine_coeffs, key_component: int, T: float, P: float, fk, specification: str, tolerance = 0.01, maxiter = 50):
    '''
    Solve flash conditions with a specified phi and T or P.

    phi: vapor flowrate to total flowrate ratio. Float
    antoine_coeffs: Array of antoine coefficients A,B,C for each key component.
    key_component: Integer. Specify index + 1 of keycomponent.
    T: Temperature, specified or guessed. In Kelvin. Float.
    P: Pressure, specified or guessed. In mmHg. Float.
    fk: Flow rates in feed, consistent unit. Array of flowrate for each key component.
    specification: A string of 'P' or 'T' specifying which parameter is specified/fixed for pressure or temperature respectively.
    tolerance: Acceptable tolerance. Float. Default is 0.01.
    maxiter: Maximum number of iterations. Default is 50. Integer.
    -- returns --
    relative volatilities, average relative volatility, vapor pressures, split fractions, Pressure, Temperature, Vapor flowrate, liquid flowrate,
    liquid compositions, vapor compositions.
    '''
    iterations = 1

    while iterations < maxiter:

        P_vaps = get_vapor_pressure(antoine_coeffs, T);
        K_n, alpha_k, eps_n, V_k, L_k, V, L, Total_Flow, x_k, y_k, alpha_bar = calc_relative_volatility(phi, P, fk, P_vaps, parameter='phi', key_component=key_component )

        majority_index = np.argmax(x_k) # get the index of the component that has the highest liquid composition

        if specification == 'P':
            # Running this code block if T was guessed
            T_calc = check_T(P, alpha_k, alpha_bar, antoine_coeffs, majority_index)

            if (np.abs(T_calc - T) <= tolerance):
                print(f'Converged after {iterations} iterations! Epsilons of each component: {eps_n}, temperature: {T} Kelvin')
                return (alpha_k, alpha_bar, P_vaps, eps_n, P, T, V, L, Total_Flow, x_k, y_k)
            else:
                T = (T + T_calc)/ 2
                iterations += 1

        elif specification == 'T':
            # Running this code block if P was guessed
            P_calc = check_P(alpha_k, alpha_bar, P_vaps, majority_index)
        
            if (np.abs(P_calc - P) <= tolerance):
                print(f'Converged after {iterations} iterations! Epsilons of each component: {eps_n}, pressure: {P} mmHg')
                return (alpha_k, alpha_bar, P_vaps, eps_n, P, T, V, L, Total_Flow, x_k, y_k)
            else:
                P = (P + P_calc)/2
                iterations += 1

        if iterations >= maxiter:
            print(f'Failed to converge')
            break;

# Test code cases down here
if __name__ == "__main__":

    class Epsilon:
        def __init__(self, value, position) -> None:
            self.value = value
            self.position = position
    
    eps = Epsilon(0.39, 2)
    P = 750
    T = 385

    antoines = np.array([[15.9008, 2788.51, -52.34], [16.0137, 3096.52, -53.67], [16.1156, 3395.57, -59.44]]);
    fk = np.array([30, 50, 40])

    # case1_solver(eps, antoines, T, P, fk, 'P', 0.1, 50 )
    # case2_solver(eps, antoines, T, P, fk, tolerance= 0.0005, maxiter=200 )
    case3_solver(0.5, antoines, 2, T, P, fk, 'P' )