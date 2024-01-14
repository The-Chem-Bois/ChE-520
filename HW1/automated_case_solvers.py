'''
Oliver Erdmann
20843970
Functions written in Python3
'''
import numpy as np

def case_1_solver(epsilon, antoine_coeffs, T: float, P: float, fk: float, specification: str, tolerance = 0.5, maxiter= 50):
    '''
    epsilon: An object with property value for epsilon as well as its key position. (i.e. epsilon for component 2 would therefore be 2).
    antoine_coeffs: Array containing arrays of each key's anotine coefficients. (Order matters, requires 3 coefficients A, B, C)
    T: The initial guessed or given temperature in Kelvin.
    P: The initial guessed or given pressure in mmHg.
    fk: Array of the fk value for each component in mol/s. (Order matters).
    specification: String of 'P' for pressure, or 'T' for temperature that indicated which parameter is specified.
    tolerance: Acceptable tolerance, default is 0.5.
    maxiter: The maximum number of iterations, default is 50.
    '''

    P_vaps = np.zeros(len(fk)) # create an array for storing all key component vapor pressures
    K_n = np.zeros(len(fk)) # create an array for storing all k values (P vapor pressure)/ (Total Pressure) for all key components
    iterations = 0

    while iterations < maxiter:

        for i in range(0, len(fk)):

            #Calculate vapor pressures
            A, B, C = antoine_coeffs[i] # get coefficients
            P_vap = np.exp(A - B/(C + T))
            P_vaps[i] = P_vap

            # Start to calculate relative volatilities (alpha)
            k = P_vap/P
            K_n[i] = k

        alpha_k = K_n / K_n[epsilon.position - 1] # calculates array of relative volatilities for each component relative to the epsilon component provided

        #Calculate epsilon for each component
        eps_n = (alpha_k * epsilon.value)/(1 + (alpha_k - 1)* epsilon.value)
        V_k = eps_n*fk
        L_k = (1-eps_n)*fk
        V = sum(V_k)
        L = sum(L_k)
        Total_Flow = V + L
        
        x_k = L_k/L
        y_k = V_k/V

        max_index = np.argmax(x_k) # get the index of the component that has the highest liquid composition

        alpha_bar = sum(x_k * alpha_k) # calculate average volatility

        if specification == 'P':
            # Running this code block if T was guessed
            P_k_vap = alpha_k[max_index]/alpha_bar * P
            # Rearrange antoine's equation and solve for T
            A,B,C = antoine_coeffs[max_index]
            T_calc = B/(A - np.log(P_k_vap)) - C

            if (np.abs(T_calc - T) <= tolerance):
                print(f'Converged after {iterations} iterations! Epsilons of each component: {eps_n}, temperature: {T} Kelvin')
                break;
            else:
                T = (T + T_calc)/ 2
                iterations += 1

        elif specification == 'T':
            pass
            # Running this code block if P was guessed
            P_calc = alpha_bar/alpha_k[max_index] * P_vaps[max_index]

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
    eps1 = Epsilon(0.8, 2 )
    P = 750 ## 1 bar = 750 mmHg
    T = 390 ## Initial guess
    fk = np.array([30, 50, 40]);
    antoine_coeffs = np.array([[15.9008, 2788.51, -52.34], [16.0137, 3096.52, -53.67], [16.1156, 3395.57, -59.44]]);

    case_1_solver(eps1, antoine_coeffs, T, P, fk, 'P', 0.01);

''' 
Questions to ask Nasser

1. Do units of fk matter? Or do we only care that they are consistent with each other?
2. Extra examples so I could verify my functions?

'''