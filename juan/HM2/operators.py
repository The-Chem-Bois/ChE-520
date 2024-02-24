'''
Juan De Leon
ID: 20765671
'''
import numpy as np

def antoine(value, coefficients, output = None):
  

    A = coefficients[:, 0]
    B = coefficients[:, 1]
    C = coefficients[:, 2]
   
    if output == "P":
        for i in coefficients:
            p_vap = np.exp((A-(B/(C+value))))
        return p_vap
    
    elif output =="T":
        for i in coefficients:
            temperature = (B /(A- np.log(value)))-C
        return temperature
    else:
        print ("Error")
        return

def relative_volatiity(fk, coefficients, temp, pressure, key, xi_n = None, phi = None,):

    key_position = key-1
    p_vaps = antoine(temp, coefficients, "P")
    K_k = p_vaps/pressure
    K_n = K_k[key_position]    
    alpha_k = K_k/K_n
   
    if xi_n == None:
        theta = (K_n*phi)/(1-phi)
        xi_n = theta/(1+theta)
     
    xi_old = alpha_k*xi_n/(1+(alpha_k-1)*xi_n)
    xi = np.delete(xi_old, key_position) 
    xi_k = np.insert(xi,key_position, xi_n)  

    v_k = xi_k*fk
    l_k = (1-xi_k)*fk
    V = np.sum(v_k)
    L = np.sum(l_k)
    F = L+V
    y_k = v_k/V
    x_k = l_k/L

    average_alpha = np.sum(x_k*alpha_k)
    
    return alpha_k, average_alpha, x_k, V, F, K_n

def bubble_point(temperature, pressure, coefficients, fk, output = None):
    iterations = 1
    tol = 0.01
    maxiter = 50
    while iterations < maxiter:

        z_k = fk/np.sum(fk)
        x_k = z_k
        key_i = np.argmax(fk)

        p_vaps = antoine(temperature, coefficients, 'P')
        K_k = p_vaps/pressure
        K_n = K_k[key_i]
        alpha_k = K_k/K_n
        alpha_avg = np.sum(x_k*alpha_k)

        if output == "P":
            bp_pressure = alpha_avg*p_vaps[key_i]
            iterations +=1
            if (np.abs(bp_pressure-pressure) <= tol):
                return bp_pressure
            else:
                pressure = (bp_pressure + pressure)/2

        elif output == "T":
            n_vap = alpha_avg**-1*pressure
            bp_temp = antoine(n_vap, np.array([coefficients[key_i]]), output)
            iterations += 1
            if(np.abs(bp_temp-temperature) <=tol):
                return bp_temp
            else:
                temperature = (bp_temp+temperature)/2
    
    print('Did not converge')
    return None
    

def dew_point(temperature, pressure, coefficients, fk, output = None):
    iterations = 1
    tol = 0.01
    maxiter = 50

    while iterations < maxiter:

        z_k = fk/np.sum(fk)
        key_i = np.argmax(fk)

        p_vaps = antoine(temperature, coefficients, "P")
        K_k = p_vaps/pressure
        K_n = K_k[key_i]
        alpha_k = K_k/K_n

        suma = np.sum(z_k/alpha_k)
        
        if output =='P':
            dw_pressure = p_vaps[key_i]*suma**-1
            iterations += 1
            if (np.abs(dw_pressure-pressure) <=tol):
                return dw_pressure
            else:
                pressure = (dw_pressure + pressure)/2
       
        elif output == "T":
            n_vap = pressure * suma
            dw_temp = antoine(n_vap, np.array([coefficients[key_i]]), output)
            iterations += 1
            if(np.abs(dw_temp-temperature) <=tol):
                return dw_temp
            else:
                temperature = (dw_temp+temperature)/2
    print("Did not converge")

def case1_solver(fk, coefficients, temperature, pressure, xi_n, key, fixed = None):
 
    iterations = 1
    tol = 0.05
    maxiter = 50
    
    while (iterations <= maxiter):

        alpha_k, average_alpha, x_k, V, F, K_nk = relative_volatiity(fk, coefficients, temperature, pressure, key, xi_n)
        heavy_i = np.argmax(x_k)-1
        alpha_n = alpha_k[heavy_i]

        if fixed == "P":
            p_vap_k = (alpha_n/average_alpha)*pressure
            fix_temp = antoine(p_vap_k,np.array([coefficients[heavy_i]]),"T")
            
            iterations += 1
            if (np.abs(fix_temp-temperature)<= tol):
                return fix_temp
            else:
                temperature = (fix_temp+temperature)/2
       
        if fixed == "T":
            p_vap_k = antoine(temperature, np.array([coefficients[heavy_i]]), "P")
            fix_pressure = (average_alpha/alpha_n)*p_vap_k
            
            iterations += 1
            if (np.abs(fix_pressure-pressure)<= tol):
                return fix_pressure
            else:
                pressure = (fix_pressure+pressure)/2
       
    print("Did not converge")
    return

def case2_solver (fk, coefficients, temperature, pressure, xi_n, key):

    iterations = 1
    tol = 0.05
    maxiter = 50
    
    while (iterations <= maxiter):

        alpha_k, average_alpha, x_k, V, F, K_n = relative_volatiity(fk, coefficients, temperature, pressure, key, xi_n)
        
        phi = V/F 
        theta = (K_n*phi)/(1-phi)
        xi_new = theta/(1+theta)
        iterations += 1
        if (np.abs(xi_new-xi_n) <= tol):
            return xi_new
        else:
            xi_n = (xi_new + xi_n)/2

    print("Did not converge")
    return None


def case3_solver(fk, coefficients, temperature, pressure, key, phi, fixed = None):

    iterations = 1
    tol = 0.05
    maxiter = 50
    
    while (iterations <= maxiter):

        alpha_k, average_alpha, x_k, V, F, K_n = relative_volatiity(fk, coefficients, temperature, pressure, key, phi)
        heavy_i = np.argmax(x_k)-1
        alpha_n = alpha_k[heavy_i]

        if fixed == "P":
            p_vap_k = (alpha_n/average_alpha)*pressure
            fix_temp = antoine(p_vap_k,np.array([coefficients[heavy_i]]),"T")
            
            iterations += 1
            if (np.abs(fix_temp-temperature)<= tol):
                return fix_temp
            else:
                temperature = (fix_temp+temperature)/2
       
        if fixed == "T":
            p_vap_k = antoine(temperature, np.array([coefficients[heavy_i]]), "P")
            fix_pressure = (average_alpha/alpha_n)*p_vap_k
            
            iterations += 1
            if (np.abs(fix_pressure-pressure)<= tol):
                return fix_pressure
            else:
                pressure = (fix_pressure+pressure)/2
    print("Did not converge")
    return None




if __name__ == "__main__":
    '''
    Testing
    
    coefficients = np.array([[15.9008,2788.51,-52.34],[16.0137,3096.52,-53.67],[16.1156,3395.57,-59.44]]) 
    fk = np.array([[30,50,40]])
    print(antoine(326.40,coefficients, output="P"))
    


    print(relative_volatiity(fk,coefficients,380,750,2,phi=0.5))
    print(bubble_point(310, 51379.247, coefficients, fk, output="T"))
    print(antoine(310, coefficients, output="T"))
    print(dew_point(390, 750, coefficients, fk, output="T"))
    print(antoine(390, np.array(coefficients[1]), output="T"))
    print(case1_solver(fk, coefficients, 391.5915, 1050,0.8,2,"T"))
    print(case2_solver(fk, coefficients, 385, 750,0.7,2))
    print(case3_solver(fk, coefficients, 380, 750,2,0.5,"P"))
    
    '''
    


    
