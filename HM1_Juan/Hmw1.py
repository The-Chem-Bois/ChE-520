'''
Juan De Leon
ID: 20765671
'''
import numpy as np
#Task 1

def antoine(value, coefficients, output = None):
    # convert from log10->ln
    # coefficients[:, 0] *= np.log(10.)
    # coefficients[:, 1] *= np.log(10.)

    # convert Celsius->Kelvin
    # coefficients[:, 2] -= 273.15

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

#Task 2
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
    
    return alpha_k, average_alpha

def bubble_point(temperature, pressure, coefficients, fk, output = None):
    iterations = 1
    tol = 0.01
    maxiter = 50
    while iterations < maxiter:

        xi_k = 0
        l_k = fk
        z_k = fk/np.sum(fk)
        x_k = z_k
        key_i = np.argmax(fk)

        if output == "P":
            p_vaps = antoine(temperature, coefficients, output)
            K_k = p_vaps/pressure
            K_n = K_k[key_i]

            alpha_k = K_k/K_n
            alpha_avg = np.sum(x_k*alpha_k)
            bp_pressure = alpha_avg*p_vaps[key_i]
            iterations +=1
            if ((bp_pressure-pressure) <= tol):
                return bp_pressure
            else:
                pressure = (bp_pressure + pressure)/2
            
        else:
            pass
    
    print('Did not converge')
    return None
    
    # elif output == T:
    #     iterations = 1
    #     tol = 0.01
    #     maxiter = 50
    #     while iterations < maxiter:
    # return
    
    # else:
    #     print("Did not converge")




    #Task 4



    #Task 5

#  def case1_solver(fk, coefficients, pressure, temp, xi_n, key):

#     p_vaps = antoine(66)
#     return








if __name__ == "__main__":

    coefficients = np.array([[15.9008,2788.51,-52.34],[16.0137,3096.52,-53.67],[16.1156,3395.57,-59.44]]) 
    fk = np.array([[30,50,40]])
    antoine(380,coefficients, output="P")
    

    # print(relative_volatiity(fk,coefficients,390,750,2,0.8))
    print(bubble_point(390, 750, coefficients, fk, output="P"))
    
