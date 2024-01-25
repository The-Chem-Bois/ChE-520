#CHE 520 Homework 1 Computer Augmentation
#Christina Mohan
#20844467

#Import numpy to access required functions
import numpy as np

#Define a function to compute the vapour pressure or temperature
def Antoine(A, B, C, T_or_P_val, Find_T = False):
    '''
    Parameters
    ----------
    A : FLOAT
        ANTOINE EQUATION PARAMETER A.
    B : FLOAT
        ANTOINE EQUATION PARAMETER B.
    C : FLOAT
        ANTOINE EQUATION PARAMETER C.
    T_or_P_val : FLOAT
        TEMPERATURE AT A GIVEN VAPOUR PRESSURE IN K OR VAPOUR PRESSURE IN MMHG
    Find_T : BOOLEAN, optional
        WHETHER YOU ARE SOLVING FOR TEMPERATURE OR NOT. The default value is False

    Returns
    -------
    Tvap : FLOAT 
        TEMPERATURE AT A GIVEN VAPOUR PRESSURE IN K. 
    Pvap : FLOAT 
        VAPOUR PRESSURE AT A GIVEN TEMPERATURE IN mmHg. 
    '''
    
    #Check if temperature or pressure is to be computed. If T is False, then pressure will be computed.
    if Find_T == False:
        T = T_or_P_val
        Pvap = np.exp(A-(B/(C+T))) #Computes the vapour pressure in mmHg from Antoine Equation
        return(Pvap) 
    
    #If T is True, then temperature will be computed.
    elif Find_T == True:
        P = T_or_P_val
        Tvap = (B/(A-np.log(P)))-C #Computes the corresponding temperature in K from Antoine Equation
        return(Tvap)
    
    else:
        return("Invalid Input")

#Test cases    
#print(Antoine(15.9008,2788.51,-52.34,390))
#print(Antoine(16.1156,3395.57,-59.44,361,T = True))

#Define a function to compute the relative volatility of each component
def relative_volatility(P, Pvap_k, Pvap_n):
    '''
    Parameters
    ----------
    P : FLOAT
        PRESSURE IN mmHg.
    Pvap_k : ARRAY
        ARRAY OF VAPOUR PRESSURES OF EACH COMPONENT IN mmHg.
    Pvap_n : FLOAT
        VAPOUR PRESSURE OF KEY COMPONENT IN mmHg.

    Returns
    -------
    K_n : FLOAT
        VOLATILITY OF THE KEY COMPONENT
    relative_volatility_k : ARRAY
        ARRAY OF RELATIVE VOLATILITIES OF EACH COMPONENT.
    '''
    #print(Pvap_k)
    print (Pvap_n)
    
    K_k = np.array(Pvap_k)/P #Compute K_k for each component
    K_n = Pvap_n/P #Compute K_n for the key component
    relative_volatility_k = K_k/K_n #Compute the relative volatility of each component
    
    return(K_n, relative_volatility_k)

#Test cases
#print(relative_volatility(750,[2084.9,904.1,345],904.1))

#Define a function to compute the liquid and vapor molar flowrates and composition
def average_volatility(relative_volatility, fk, eps_n):
    '''    
    Parameters
    ----------
    relative_volatility : ARRAY
        ARRAY OF RELATIVE VOLATILITIES OF EACH COMPONENT.
    fk : ARRAY
        ARRAY OF MOLAR FLOWRATES OF EACH COMPONENTS (UNITS ARE INDIFFERENT)
    eps_n : FLOAT
        OVERHEAD SPLIT FRACTION FOR KEY COMPONENT.

    Returns
    -------
    avg_volatility : FLOAT
        AVERAGE VOLATILITY OF MIXTURE
    xk : ARRAY
        ARRAY OF LIQUID MOLE FRACTIONS
    yk : ARRAY
        ARRAY OF VAPOUR MOLE FRACTIONS
    V: FLOAT
        TOTAL VAPOUR MOLAR FLOWRATE (UNITS ARE INDIFFERENT)
    L: FLOAT
        TOTAL LIQUID MOLAR FLOWRATE (UNITS ARE INDIFFERENT)

    '''
    eps_k = (np.array(relative_volatility)*eps_n)/(1+(np.array(relative_volatility)-1)*eps_n) #Compute the overhead split fraction for each component
    
    #Initialize the empty arrays for the vapour and liquid molar flowrates
    vk = np.zeros(np.shape(fk))
    lk = np.zeros(np.shape(fk))
    
    for i in range(len(fk)):
        vk[i] = eps_k[i]*fk[i] #Calculate the vapour molar flowrate for each component
        lk [i] = (1-eps_k[i])*fk[i] #Calculate the liquid molar flowrate for each component
        
    V = np.sum(vk) #Compute the sum of the component vapour flowrates
    L = np.sum(lk) #Compute the sum of the component liquid flowrates
    F = V + L #Calculate the total feed flowrate
    
    #Check that the mass balance is consistent with the total feed flowrate
    #print(np.sum(fk))
    #print(F)
    
    if V == 0:
        xk = lk/L #Compute the liquid mole fractions of each component
        yk = np.zeros(np.shape(fk))
    
    elif L == 0:
        yk = vk/V #Compute the vapour mole fractions of each component
        xk = np.zeros(np.shape(fk))
    
    else:
        xk = lk/L #Compute the liquid mole fractions of each component
        yk = vk/V #Compute the vapour mole fractions of each component
        
    volatility = np.zeros(np.shape(xk)) #Initialize an empty array for the product of the liquid mole fraction and relative volatility
    
    #Compute the product of the liquid mole fraction and relative volatility for each component
    for j in range(len(xk)):
        volatility[j] = xk[j]*relative_volatility[j]
        
    avg_volatility = np.sum(volatility) #Compute the average volatility by summing up each product
    
    return(avg_volatility, xk, yk, V, L)

#Test cases
#print(average_volatility([2.30605022,1,0.38159496],[30,50,40],0.8))

#Define a function to compute the bubble point temperature or pressure
def bubble_point(T_or_P_val, T_or_P_val_init, fk, A, B, C, n, Find_T = False):
    '''
    
    Parameters
    ----------
    T_or_P_val : FLOAT
        TEMPERATURE IN K OR PRESSURE IN MMHG.
    T_or_P_val_init : FLOAT
        INITIAL GUESS FOR TEMPERATURE IN K OR PRESSURE IN MMHG
    fk : ARRAY
        ARRAY OF MOLAR FLOWRATES OF EACH COMPONENTS (UNITS ARE INDIFFERENT).
    A : ARRAY
        ANTOINE EQUATION PARAMETER A.
    B : ARRAY
        ANTOINE EQUATION PARAMETER B.
    C : ARRAY
        ANTOINE EQUATION PARAMETER C.
    n : INTEGER
        INDEX OF KEY COMPONENT. INDICES START FROM 0.
    Find_T : BOOLEAN, optional
        WHETHER YOU ARE SOLVING FOR TEMPERATURE OR NOT. The default is False.

    Returns
    -------
    T_new : FLOAT
        FINAL BUBBLE POINT TEMPERATURE IN K
    P_new : FLOAT
        FINAL BUBBLE POINT PRESSURE IN mmHg

    '''

    tol = 1e-7 #Set a tolerance for the final solution
    eps_n = 0 #For bubble point, the eps values are equal to zero
    iteration = 0 #Initialize a counter
    
    #If Find_T = False, we are given temperature and need to guess an initial pressure
    if Find_T == False:
        T = T_or_P_val #Define the temperature from the inputs
        P0 = T_or_P_val_init #Define the initial pressure guess from the inputs
        
        error_P = 1 #Initialize error for pressure value
        
        #Initialize empty vapour pressure array
        Pvap = np.zeros(np.shape(A))
        #Calculate the vapour pressure for each component
        for k in range(len(Pvap)):
            Pvap[k] = Antoine(A[k], B[k], C[k], T)
        
        #While the absolute error between the pressures is greater than the tolerance, continue to iterate
        while error_P > tol:
            
            iteration += 1
            rv = relative_volatility(P0, Pvap, Pvap[n])
            av = average_volatility(rv[1], fk, eps_n)
            
            avg_alpha = av[0] #Compute the average relative volatility of the mixture
            xk = av[1] #Compute the liquid mole fractions of each component
            alpha = rv[1] #Compute the relative volatility of the components
            
            x_ind = np.argmax(xk) #Determine component with maximum feed composition
            
            P_new = (avg_alpha*Pvap[x_ind])/alpha[x_ind] #Compute the bubble point pressure
            
            error_P = abs(P0-P_new) #Calculate the absolute error between the initial guess and calculated pressure
            
            P0 = P_new #Set the initial guess equal to the calculated pressure
            
        print("After",iteration,"iterations, the pressure has converged to", P_new, "mmHg.")   
        return(P_new)
    
    #If Find_T = True, we are given pressure and need to guess an initial temperature
    elif Find_T == True:
        P = T_or_P_val #Define the pressure from the inputs
        T0 = T_or_P_val_init #Define the initial temperature guess from the inputs
        
        error_T = 1 #Initialize error for temperature value
                
        #While the absolute error between the temperatures is greater than the tolerance, continue to iterate
        while error_T > tol:
            
            iteration += 1
            
            #Initialize empty vapour pressure array
            Pvap = np.zeros(np.shape(A))
            #Calculate the vapour pressure for each component
            for k in range(len(Pvap)):
                Pvap[k] = Antoine(A[k], B[k], C[k], T0)
            
            rv = relative_volatility(P, Pvap, Pvap[n])
            av = average_volatility(rv[1], fk, eps_n)
                        
            avg_alpha = av[0] #Compute the average relative volatility of the mixture
            xk = av[1] #Compute the liquid mole fractions of each component
            alpha = rv[1] #Compute the relative volatility of the components
            x_ind = np.argmax(xk) #Determine component with maximum feed composition
            
            Pvap_n = P/avg_alpha #Compute the bubble point vapor pressure
            
            T_new = Antoine(A[x_ind], B[x_ind], C[x_ind], Pvap_n, Find_T = True) #Compute the bubble point temperature
            
            error_T = abs(T0-T_new) #Calculate the absolute error between the initial guess and calculated temperature
            
            T0 = T_new #Set the initial guess equal to the calculated temperature
            
        print("After",iteration,"iterations, the temperature has converged to", T_new, "K.")   
        return(T_new)
    
    else:
        print("Invalid Input")
        
#Test Cases
#print(bubble_point(391, 390, [30,50,40],[15.9008,16.0137,16.1156],[2788.51,3096.52,3395.57],[-52.34,-53.67,-59.44],1))
#print(bubble_point(750, 390, [30,50,40],[15.9008,16.0137,16.1156],[2788.51,3096.52,3395.57],[-52.34,-53.67,-59.44],1, Find_T=True))

#Define a function to compute the bubble point temperature or pressure
def dew_point(T_or_P_val, T_or_P_val_init, fk, A, B, C, n, Find_T = False):
    '''
    
    Parameters
    ----------
    T_or_P_val : FLOAT
        TEMPERATURE IN K OR PRESSURE IN MMHG.
    T_or_P_val_init : FLOAT
        INITIAL GUESS FOR TEMPERATURE IN K OR PRESSURE IN MMHG
    fk : ARRAY
        ARRAY OF MOLAR FLOWRATES OF EACH COMPONENTS (UNITS ARE INDIFFERENT).
    A : ARRAY
        ANTOINE EQUATION PARAMETER A.
    B : ARRAY
        ANTOINE EQUATION PARAMETER B.
    C : ARRAY
        ANTOINE EQUATION PARAMETER C.
    n : INTEGER
        INDEX OF KEY COMPONENT. INDICES START FROM 0.
    Find_T : BOOLEAN, optional
        WHETHER YOU ARE SOLVING FOR TEMPERATURE OR NOT. The default is False.

    Returns
    -------
    T_new : FLOAT
        FINAL DEW POINT TEMPERATURE IN K
    P_new : FLOAT
        FINAL DEW POINT PRESSURE IN mmHg

    '''

    tol = 1e-7 #Set a tolerance for the final solution
    eps_n = 1 #For dew point, the eps values are equal to zero
    iteration = 0 #Initialize a counter
    
    #If Find_T = False, we are given temperature and need to guess an initial pressure
    if Find_T == False:
        T = T_or_P_val #Define the temperature from the inputs
        P0 = T_or_P_val_init #Define the initial pressure guess from the inputs
        
        error_P = 1 #Initialize error for pressure value
        
        #Initialize empty vapour pressure array
        Pvap = np.zeros(np.shape(A))
        #Calculate the vapour pressure for each component
        for k in range(len(Pvap)):
            Pvap[k] = Antoine(A[k], B[k], C[k], T)
        
        #While the absolute error betweent the pressures is greater than the tolerance, continue to iterate
        while error_P > tol :
            
            iteration += 1
            rv = relative_volatility(P0, Pvap, Pvap[n])
            av = average_volatility(rv[1], fk, eps_n)
            
            alpha = rv[1] #Compute the relative volatility of the mixture
            yk = av[2] #Compute the liquid mole fractions of each component
            
            y_ind = np.argmax(yk) #Determine component with maximum feed composition
            
            #Computes the summation term to calculate the new pressure
            summation = 0
            for i in range(len(yk)):
                summation += yk[i]/alpha[i]
            
            P_new = Pvap[y_ind]/(summation) #Compute the dew point pressure
            
            error_P = abs(P0 - P_new) #Calculate the absolute error between the initial guess and calculated pressure
            
            P0 = P_new #Set the initial guess equal to the calculated pressure
            
        print("After",iteration,"iterations, the pressure has converged to", P_new, "mmHg.")   
        return(P_new)
    
    #If Find_T = True, we are given pressure and need to guess an initial temperature
    elif Find_T == True:
        P = T_or_P_val #Define the pressure from the inputs
        T0 = T_or_P_val_init #Define the initial temperature guess from the inputs
        
        error_T = 1 #Initialize error for temperature value
                
        #While the absolute error between the temperatures is greater than the tolerance, continue to iterate
        while error_T > tol :
            
            iteration += 1
            
            #Initialize empty vapour pressure array
            Pvap = np.zeros(np.shape(A))
            #Calculate the vapour pressure for each component
            for k in range(len(Pvap)):
                Pvap[k] = Antoine(A[k], B[k], C[k], T0)
            
            rv = relative_volatility(P, Pvap, Pvap[n])
            av = average_volatility(rv[1], fk, eps_n)
            
            alpha = rv[1] #Compute the average relative volatility of the mixture
            yk = av[2] #Compute the liquid mole fractions of each component
            
            y_ind = np.argmax(yk) #Determine component with maximum feed composition
            
            #Computes the summation term to calculate the new vapor pressure
            summation = 0
            for i in range(len(yk)):
                summation += yk[i]/alpha[i]
            
            Pvap_n = P*(summation) #Compute the dew point vapor pressure
            
            T_new = Antoine(A[y_ind], B[y_ind], C[y_ind], Pvap_n, Find_T = True) #Compute the dew point temperature
            
            error_T = abs(T0-T_new) #Calculate the absolute error between the initial guess and calculated temperature
            
            T0 = T_new #Set the initial guess equal to the calculated temperature
            
        print("After",iteration,"iterations, the temperature has converged to", T_new, "K.")   
        return(T_new)
    
    else:
        print("Invalid Input")
        
#Test Cases
#print(dew_point(391, 390, [30,50,40],[15.9008,16.0137,16.1156],[2788.51,3096.52,3395.57],[-52.34,-53.67,-59.44],1))
#print(dew_point(750, 200, [30,50,40],[15.9008,16.0137,16.1156],[2788.51,3096.52,3395.57],[-52.34,-53.67,-59.44],1, Find_T=True))

#Define a function to perform case 1 flash calculations
def case1_flash(eps_n, T_or_P_val, T_or_P_val_init, fk, A, B, C, n, Find_T = True):
    '''
    
    Parameters
    ----------
    eps_n : FLOAT
        OVERHEAD SPLIT FRACTION FOR KEY COMPONENT.
    T_or_P_val : FLOAT
        TEMPERATURE IN K OR PRESSURE IN MMHG.
    T_or_P_val_init : FLOAT
        INITIAL GUESS FOR TEMPERATURE IN K OR PRESSURE IN MMHG
    fk : ARRAY
        ARRAY OF MOLAR FLOWRATES OF EACH COMPONENTS (UNITS ARE INDIFFERENT).
    A : ARRAY
        ANTOINE EQUATION PARAMETER A.
    B : ARRAY
        ANTOINE EQUATION PARAMETER B.
    C : ARRAY
        ANTOINE EQUATION PARAMETER C.
    n : INTEGER
        INDEX OF KEY COMPONENT. INDICES START FROM 0.
    Find_T : BOOLEAN, optional
        WHETHER YOU ARE SOLVING FOR TEMPERATURE OR NOT. The default is True.

    Returns
    -------
    T_new : FLOAT
        FINAL TEMPERATURE IN K
    P_new : FLOAT
        FINAL PRESSURE IN mmHg
    xk : ARRAY
        ARRAY OF LIQUID MOLE FRACTIONS
    yk : ARRAY
        ARRAY OF VAPOUR MOLE FRACTIONS
    V: FLOAT
        TOTAL VAPOUR MOLAR FLOWRATE (UNITS ARE INDIFFERENT)
    L: FLOAT
        TOTAL LIQUID MOLAR FLOWRATE (UNITS ARE INDIFFERENT)

    '''
    tol = 1e-7 #Set a tolerance for the final solution
    iteration = 0 #Initialize a counter
    
    #If Find_T = True, we are given pressure and need to guess an initial temperature
    if Find_T == True:
        
        P = T_or_P_val #Define the pressure from the inputs
        T0 = T_or_P_val_init #Define the initial temperature guess from the inputs
        #Initialize errors for average volatility and temperature
        error_T = 1
        error_a = 1
        #While the errors for alpha and T are larger than the tolerance, continue to iterate
        while (error_T > tol) and (error_a >tol):
            iteration += 1
            #Initialize empty vapour pressure array
            Pvap = np.zeros(np.shape(A))
            #Calculate the vapour pressure for each component
            for k in range(len(Pvap)):
                Pvap[k] = Antoine(A[k], B[k], C[k], T0)
            
            rv = relative_volatility(P,Pvap,Pvap[n])
            Kn = rv[0] #Compute the K value for the key component 
            alpha = rv[1] #Compute the relative volatility array for all components
            
            av = average_volatility(alpha,fk,eps_n) 
            avg_alpha_act = av[0] #Compute average volatility
            xk = av[1] #Compute liquid mole fractions for each component
            yk = av[2] #Compute vapour mole fractions for each component
            V = av[3] #Compute total vapour molar flowrate
            L = av[4] #Compute total liquid molar flowrate
            avg_alpha_theor = 1/Kn #Calculate expected average volatility
            
            x_ind = np.argmax(xk) #Determine component with maximum liquid composition
            
            Pvl = (alpha[x_ind]/avg_alpha_act)*P #Calculate new vapour pressure of said component
            T_new = (B[x_ind]/(A[x_ind]-np.log(Pvl)))-C[x_ind] #Calculate new temperature for the specified vapour pressure
                    
            error_T = abs(T0 - T_new) #Calculate the absolute error between the temperatures
            
            error_a = abs(avg_alpha_act - avg_alpha_theor) #Calculate the absolute error between the average volatilities
            
            T0 = T_new #Set the initial guess equal to the calculated temperature
            
        print("After",iteration,"iterations, the temperature has converged to", T_new, "K.")    
        return(T_new, xk, yk, V, L)
    
    #If Find_T = False, we are given temperature and need to guess an initial pressure
    elif Find_T == False:
        
        T = T_or_P_val #Define the temperature from the inputs
        P0 = T_or_P_val_init #Define the initial pressure guess from the inputs
        
        #Initialize errors for average volatility and pressure
        error_P = 1
        error_a = 1
        
        #Initialize empty vapour pressure array
        Pvap = np.zeros(np.shape(A))
        #Calculate the vapour pressure for each component
        for k in range(len(Pvap)):
            Pvap[k] = Antoine(A[k], B[k], C[k], T)
            
        #While the errors for alpha and P are larger than the tolerance, continue to iterate
        while (error_P > tol) and (error_a > tol):
            iteration += 1
            
            rv = relative_volatility(P0,Pvap,Pvap[n])
            Kn = rv[0] #Compute the K value for the key component 
            alpha = rv[1] #Compute the relative volatility array for all components
            
            av = average_volatility(alpha,fk,eps_n) 
            avg_alpha_act = av[0] #Compute average volatility
            xk = av[1] #Compute liquid mole fractions for each component
            yk = av[2] #Compute vapour mole fractions for each component
            V = av[3] #Compute total vapour molar flowrate
            L = av[4] #Compute total liquid molar flowrate
            avg_alpha_theor = 1/Kn #Calculate expected average volatility
            
            x_ind = np.argmax(xk) #Determine component with maximum liquid composition
            
            P_new = (avg_alpha_act/alpha[x_ind])*Pvap[x_ind] #Calculate new  pressure of said component
                    
            error_P = abs(P0 - P_new) #Calculate the absolute error between the pressures
            
            error_a = abs(avg_alpha_act - avg_alpha_theor)
            
            P0 = P_new #Set the initial guess equal to the calculated pressure
            
        print("After",iteration,"iterations, the pressure has converged to", P_new, "mmHg.")    
        return(P_new, xk, yk, V, L)
    else:
        return("Invalid input")

#Test Case
#print(case1_flash(0.8,750,390,[30,50,40],[15.9008,16.0137,16.1156],[2788.51,3096.52,3395.57],[-52.34,-53.67,-59.44],1))
#print(case1_flash(0.8,390,760,[30,50,40],[15.9008,16.0137,16.1156],[2788.51,3096.52,3395.57],[-52.34,-53.67,-59.44],1, Find_T = False))

#Define a function to perform case 2 flash calculations
def case2_flash(T, P, fk, A, B, C, n):
    '''
    
    Parameters
    ----------
    T : FLOAT
        TEMPERATURE IN K.
    P : FLOAT
        PRESSURE IN MMHG
    fk : ARRAY
        ARRAY OF MOLAR FLOWRATES OF EACH COMPONENTS (UNITS ARE INDIFFERENT).
    A : ARRAY
        ANTOINE EQUATION PARAMETER A.
    B : ARRAY
        ANTOINE EQUATION PARAMETER B.
    C : ARRAY
        ANTOINE EQUATION PARAMETER C.
    n : INTEGER
        INDEX OF KEY COMPONENT. INDICES START FROM 0.
    
    Returns
    -------
    eps_n : FLOAT
        VAPOR SPLIT FRACTION
    xk : ARRAY
        ARRAY OF LIQUID MOLE FRACTIONS
    yk : ARRAY
        ARRAY OF VAPOUR MOLE FRACTIONS
    V: FLOAT
        TOTAL VAPOUR MOLAR FLOWRATE (UNITS ARE INDIFFERENT)
    L: FLOAT
        TOTAL LIQUID MOLAR FLOWRATE (UNITS ARE INDIFFERENT)

    '''
    tol = 1e-7 #Set a tolerance for the final solution
    iteration = 0 #Initialize a counter
    eps0 = 0.5 #Initial guess for split fraction
    error_eps = 1 #Initialize error for split fraction
    
    #Initialize empty vapour pressure array
    Pvap = np.zeros(np.shape(A))
    #Calculate the vapour pressure for each component
    for k in range(len(Pvap)):
        Pvap[k] = Antoine(A[k], B[k], C[k], T)
    
    rv = relative_volatility(P,Pvap,Pvap[n])
    Kn = rv[0] #Compute the K value for the key component 
    alpha = rv[1] #Compute the relative volatility array for all components
    
    #While the error for the split fraction is larger than the tolerance, continue to iterate
    while (error_eps > tol):
        iteration += 1
            
        av = average_volatility(alpha,fk,eps0) 
        xk = av[1] #Compute liquid mole fractions for each component
        yk = av[2] #Compute vapour mole fractions for each component
        V = av[3] #Compute total vapour molar flowrate
        L = av[4] #Compute total liquid molar flowrate
             
        phi = V/(V+L) #Calculate the vaporization fraction
        theta = (Kn*phi)/(1-phi) 
        eps_new = theta/(1+theta) #Calculate the new split fraction value
        
        error_eps = abs(eps0 - eps_new) #Calculate the absolute error between the split fractions
        
        eps0 = eps_new #Set the initial guess equal to the calculated split fraction
            
    print("After",iteration,"iterations, the split fraction has converged to", eps_new, ".")    
    return(eps_new, xk, yk, V, L)

#Test Case
#print(case2_flash(385, 750,[30,50,40],[15.9008,16.0137,16.1156],[2788.51,3096.52,3395.57],[-52.34,-53.67,-59.44],1))

#Define a function to perform case 3 flash calculations
def case3_flash(phi, T_or_P_val, T_or_P_val_init, fk, A, B, C, n, Find_T = True):
    '''
    
    Parameters
    ----------
    phi : FLOAT
        VAPOR SPLIT FRACTION FOR KEY COMPONENT.
    T_or_P_val : FLOAT
        TEMPERATURE IN K OR PRESSURE IN MMHG.
    T_or_P_val_init : FLOAT
        INITIAL GUESS FOR TEMPERATURE IN K OR PRESSURE IN MMHG
    fk : ARRAY
        ARRAY OF MOLAR FLOWRATES OF EACH COMPONENTS (UNITS ARE INDIFFERENT).
    A : ARRAY
        ANTOINE EQUATION PARAMETER A.
    B : ARRAY
        ANTOINE EQUATION PARAMETER B.
    C : ARRAY
        ANTOINE EQUATION PARAMETER C.
    n : INTEGER
        INDEX OF KEY COMPONENT. INDICES START FROM 0.
    Find_T : BOOLEAN, optional
        WHETHER YOU ARE SOLVING FOR TEMPERATURE OR NOT. The default is True.

    Returns
    -------
    T_new : FLOAT
        FINAL TEMPERATURE IN K
    P_new : FLOAT
        FINAL PRESSURE IN mmHg
    eps_n : FLOAT
        OVERHEAD SPLIT FRACTION FOR KEY COMPONENT
    xk : ARRAY
        ARRAY OF LIQUID MOLE FRACTIONS
    yk : ARRAY
        ARRAY OF VAPOUR MOLE FRACTIONS
    V: FLOAT
        TOTAL VAPOUR MOLAR FLOWRATE (UNITS ARE INDIFFERENT)
    L: FLOAT
        TOTAL LIQUID MOLAR FLOWRATE (UNITS ARE INDIFFERENT)

    '''
    tol = 1e-7 #Set a tolerance for the final solution
    iteration = 0 #Initialize a counter
        
    #If Find_T = True, we are given pressure and need to guess an initial temperature
    if Find_T == True:
        
        P = T_or_P_val #Define the pressure from the inputs
        T0 = T_or_P_val_init #Define the initial temperature guess from the inputs
        #Initialize errors for average volatility and temperature
        error_T = 1
        error_a = 1
        
        #While the errors for alpha and T are larger than the tolerance, continue to iterate
        while (error_T > tol) and (error_a >tol):
            iteration += 1
            #Initialize empty vapour pressure array
            Pvap = np.zeros(np.shape(A))
            #Calculate the vapour pressure for each component
            for k in range(len(Pvap)):
                Pvap[k] = Antoine(A[k], B[k], C[k], T0)
            
            rv = relative_volatility(P,Pvap,Pvap[n])
            Kn = rv[0] #Compute the K value for the key component 
            alpha = rv[1] #Compute the relative volatility array for all components
            
            theta = (Kn*phi)/(1-phi) 
            eps_n = theta/(1+theta) #Calculate the new split fraction value
            
            av = average_volatility(alpha,fk,eps_n) 
            avg_alpha_act = av[0] #Compute average volatility
            xk = av[1] #Compute liquid mole fractions for each component
            yk = av[2] #Compute vapour mole fractions for each component
            V = av[3] #Compute total vapour molar flowrate
            L = av[4] #Compute total liquid molar flowrate
            avg_alpha_theor = 1/Kn #Calculate expected average volatility
            
            x_ind = np.argmax(xk) #Determine component with maximum liquid composition
            
            Pvl = (alpha[x_ind]/avg_alpha_act)*P #Calculate new vapour pressure of said component
            T_new = (B[x_ind]/(A[x_ind]-np.log(Pvl)))-C[x_ind] #Calculate new temperature for the specified vapour pressure
                    
            error_T = abs(T0 - T_new) #Calculate the absolute error between the temperatures
            
            error_a = abs(avg_alpha_act - avg_alpha_theor) #Calculate the absolute error between the average volatilities
            
            T0 = T_new #Set the initial guess equal to the calculated temperature
            
        print("After",iteration,"iterations, the temperature has converged to", T_new, "K.")    
        return(T_new, xk, yk, V, L, eps_n)
    
    #If Find_T = False, we are given temperature and need to guess an initial pressure
    elif Find_T == False:
        
        T = T_or_P_val #Define the temperature from the inputs
        P0 = T_or_P_val_init #Define the initial pressure guess from the inputs
        
        #Initialize errors for average volatility and pressure
        error_P = 1
        error_a = 1
        
        #Initialize empty vapour pressure array
        Pvap = np.zeros(np.shape(A))
        #Calculate the vapour pressure for each component
        for k in range(len(Pvap)):
            Pvap[k] = Antoine(A[k], B[k], C[k], T)
            
        #While the errors for alpha and P are larger than the tolerance, continue to iterate
        while (error_P > tol) and (error_a > tol):
            iteration += 1
            
            rv = relative_volatility(P0,Pvap,Pvap[n])
            Kn = rv[0] #Compute the K value for the key component 
            alpha = rv[1] #Compute the relative volatility array for all components
            
            theta = (Kn*phi)/(1-phi) 
            eps_n = theta/(1+theta) #Calculate the new split fraction value
            
            av = average_volatility(alpha,fk,eps_n) 
            avg_alpha_act = av[0] #Compute average volatility
            xk = av[1] #Compute liquid mole fractions for each component
            yk = av[2] #Compute vapour mole fractions for each component
            V = av[3] #Compute total vapour molar flowrate
            L = av[4] #Compute total liquid molar flowrate
            avg_alpha_theor = 1/Kn #Calculate expected average volatility
            
            x_ind = np.argmax(xk) #Determine component with maximum liquid composition
            
            P_new = (avg_alpha_act/alpha[x_ind])*Pvap[x_ind] #Calculate new  pressure of said component
                    
            error_P = abs(P0 - P_new) #Calculate the absolute error between the pressures
            
            error_a = abs(avg_alpha_act - avg_alpha_theor)
            
            P0 = P_new #Set the initial guess equal to the calculated pressure
            
        print("After",iteration,"iterations, the pressure has converged to", P_new, "mmHg.")    
        return(P_new, xk, yk, V, L, eps_n)
    else:
        return("Invalid input")
    
#Test Case
#print(case3_flash(0.5, 750, 390, [30,50,40],[15.9008,16.0137,16.1156],[2788.51,3096.52,3395.57],[-52.34,-53.67,-59.44],1))
#print(case3_flash(0.5, 387.36, 700, [30,50,40],[15.9008,16.0137,16.1156],[2788.51,3096.52,3395.57],[-52.34,-53.67,-59.44],1, Find_T = False))