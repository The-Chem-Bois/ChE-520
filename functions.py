"""These are some useful python functions that can help with automated calculations throughout the course"""
import numpy as np


# For case 1 we are given eps2 = 0.9, P = 1bar
# I wrote this function in a javascript mindset, so it's very redundant and not optimized. It still works however :)
def Case1Solver(
    eps,
    antoine_coeffs,
    T,
    P,
    fk,
    Specification,
    tolerance = 0.05,
    maxiter=50,
):
    """
    Split fraction is known (eps) and T (or P) is specified.
    Eps: is your epsilon, provide it as an object with an index value for which key element it corresponds to {eps: value, index: int}.
    antoine_coeffs: Array of Antoine Coefficients are provided by: A, B, C as an object. Order here matters. (must correspond to index set for eps)
    T: Temp in Kelvin
    P: Pressure in mmHg
    fk: Provide fk factors (order matters!). Just an array of floats.
    Specification: 'T' or 'P'
    maxiter: maximum iterations, default is 50.

    """

    # First make a guess for T or P for the provided/specified eps (user provides this)

    iteration = 0
    # Next step is to evaluate K_n, alpha_k_n
    while iteration < maxiter:
        component_index = eps["index"]
        epsilon = eps["eps"]
        P_vaps = {}
        K_n = {}
        alpha_k_n = {}
        epsilon_k = {}
        v_k = {}
        l_k = {}
        x_k = {}
        y_k = {}
        alpha_prod_list = []
        T_calc = None
        P_calc = None
        alpha_bar = None
        L = None
        V = None


        for index, coeffs in enumerate(antoine_coeffs, start=1):
            A = coeffs["A"]
            B = coeffs["B"]
            C = coeffs["C"]
            P_vap = np.exp(A - B / (T + C))
            P_vaps[index] = P_vap

        for index, P_vap in enumerate(P_vaps.values(), start=1):
            K = P_vap / P
            K_n[index] = K
            alpha = P_vap / P_vaps[component_index]
            alpha_k_n[index] = alpha
            # must find epsilon_k
            eps_k = (alpha * epsilon) / (1 + (alpha - 1) * epsilon)
            epsilon_k[index] = eps_k
            v = eps_k * fk[index - 1]
            v_k[index] = v
            l = (1 - eps_k) * fk[index - 1]
            l_k[index] = l

        L = sum(l_k.values()) # is in kmol/hr if all units provided are correct
        V = sum(v_k.values())

        for i in range (1 , len(v_k)+ 1):
            x = l_k[i]/L
            x_k[i] = x
            y = v_k[i]/V
            y_k[i] = y
        
        Total_Flow = L + V

        for i in range (1, len(x_k) + 1):
            alpha_prod = x_k[i] *  alpha_k_n[i]
            alpha_prod_list.append(alpha_prod)

        alpha_bar = sum(alpha_prod_list)

        #get the key for the maximum liquid composition
        max_key, max_value = max(x_k.items(), key=lambda x: x[1])

        if Specification == "T":
            P_k_vap = alpha_k_n[max_key]/alpha_bar * P

            ant_coef = antoine_coeffs[max_key - 1]
            T_calc = ant_coef["B"]/(ant_coef["A"] - np.log(P_k_vap) ) - ant_coef["C"]
            if np.abs(T_calc - T ) <= tolerance:
                print(f"Converged, found temperature {T_calc} K, after {iteration} iterations")
                break;
            else:
                T = (T_calc + T)/2
                iteration += 1

        elif Specification == "P":
            ant_coef = antoine_coeffs[max_key -1]
            P_k_vap = np.exp(ant_coef["A"] + ant_coef["B"]/(ant_coef["C"] + T))

            if np.abs(P_k_vap - P) <= tolerance:
                print(f"Converged, found pressure {P_k_vap} K, after {iteration} iterations")
                break;
            else:
                P = (P_k_vap + P)/2
                iteration += 1
        
        
        if iteration >= maxiter:
            print(f'No result found, diverged at {iteration} iterations. Most recent calculated temperature: {T_calc} K')
            break;

if __name__ == "__main__":
    eps2 = {"eps": 0.9, "index": 2}
    coeffs = [
        {"A": 15.9008, "B": 2788.51, "C": -52.34},
        {"A": 16.0137, "B": 3096.52, "C": -53.67},
        {"A": 16.1156, "B": 3395.7, "C": -59.44},
    ]
    fk = [30, 50, 40]
    T = 390
    P = 750
    Spec = "T"

    Case1Solver(eps2, coeffs, T, P, fk, Spec)
