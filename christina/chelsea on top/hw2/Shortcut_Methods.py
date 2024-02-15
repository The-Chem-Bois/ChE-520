#Shortcut Flash - do not know feed compositions

import numpy as np
from scipy import *

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

def bubble_point(T_or_P_val, T_or_P_val_init, fk, A, B, C, Find_T = False):
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
        
        #Initialize empty vapour pressure and feed mole fraction array
        Pvap = np.zeros(np.shape(A))
        z = np.zeros(np.shape(A))
        #Calculate the vapour pressure for each component
        for k in range(len(Pvap)):
            Pvap[k] = Antoine(A[k], B[k], C[k], T)
            z[k] = fk[k]/np.sum(fk)
            
        n = np.argmax(fk) #Determine most abundant component in the feed
        
        #While the absolute error between the pressures is greater than the tolerance, continue to iterate
        while error_P > tol:
            
            K = Pvap/P0 #Compute volatility of each component
                        
            alpha = K/K[n] #Compute relative volatility of each component
            
            avg_alpha = sum(z*alpha) #Compute average relative volatility
            
            P_new = (avg_alpha*Pvap[n])/alpha[n] #Compute the bubble point pressure
            
            error_P = abs(P0 - P_new) #Calculate the absolute error between the initial guess and calculated pressure
            
            P0 = P_new #Set the initial guess equal to the calculated pressure
            
        print("After",iteration,"iterations, the bubble point pressure has converged to", P_new, "mmHg.")   
        return(P_new)
    
    #If Find_T = True, we are given pressure and need to guess an initial temperature
    elif Find_T == True:
        P = T_or_P_val #Define the pressure from the inputs
        T0 = T_or_P_val_init #Define the initial temperature guess from the inputs
        
        error_T = 1 #Initialize error for temperature value
                
        #While the absolute error between the temperatures is greater than the tolerance, continue to iterate
        while error_T > tol:
            
            iteration += 1
            
            #Initialize empty vapour pressure and feed mole fraction array
            Pvap = np.zeros(np.shape(A))
            z = np.zeros(np.shape(A))
            #Calculate the vapour pressure for each component
            for k in range(len(Pvap)):
                Pvap[k] = Antoine(A[k], B[k], C[k], T0)
                z[k] = fk[k]/np.sum(fk)
                
            n = np.argmax(fk) #Determine most abundant component in the feed
            
            K = Pvap/P #Compute relative volatility of each component
            
            alpha = K/K[n] #Compute relative volatility of each component
            
            avg_alpha = sum(z*alpha) #Compute average relative volatility
            
            Pvap_n = P/avg_alpha #Compute the bubble point vapor pressure
            
            T_new = Antoine(A[n], B[n], C[n], Pvap_n, Find_T = True) #Compute the bubble point temperature
            
            error_T = abs(T0 - T_new) #Calculate the absolute error between the initial guess and calculated temperature
            
            T0 = T_new #Set the initial guess equal to the calculated temperature
            
        print("After",iteration,"iterations, the bubble point temperature has converged to", T_new, "K.")   
        return(T_new)
    
    else:
        print("Invalid Input")

#Define a function to compute the bubble point temperature or pressure
def dew_point(T_or_P_val, T_or_P_val_init, fk, A, B, C, Find_T = False):
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
        
        #Initialize empty vapour pressure and feed mole fraction array
        Pvap = np.zeros(np.shape(A))
        z = np.zeros(np.shape(A))
        #Calculate the vapour pressure for each component
        for k in range(len(Pvap)):
            Pvap[k] = Antoine(A[k], B[k], C[k], T)
            z[k] = fk[k]/np.sum(fk)
            
        n = np.argmax(fk) #Determine most abundant component in the feed
        
        #While the absolute error betweent the pressures is greater than the tolerance, continue to iterate
        while error_P > tol :
            
            iteration += 1
            
            K = Pvap/P0 #Compute relative volatility of each component
            
            alpha = K/K[n] #Compute relative volatility of each component
                      
            #Computes the summation term to calculate the new pressure
            summation = sum(z/alpha)
            
            P_new = Pvap[n]/(summation) #Compute the dew point pressure
            
            error_P = abs(P0 - P_new) #Calculate the absolute error between the initial guess and calculated pressure
            
            P0 = P_new #Set the initial guess equal to the calculated pressure
            
        print("After",iteration,"iterations, the dew point pressure has converged to", P_new, "mmHg.")   
        return(P_new)
    
    #If Find_T = True, we are given pressure and need to guess an initial temperature
    elif Find_T == True:
        P = T_or_P_val #Define the pressure from the inputs
        T0 = T_or_P_val_init #Define the initial temperature guess from the inputs
        
        #Given denominator term, create function for fsolve that outputs the difference in pressures        
        def func(T):
        #Initialize empty vapour pressure and feed mole fraction array
            Pvap = np.zeros(np.shape(A))
            z = np.zeros(np.shape(A))
            #Calculate the vapour pressure for each component
            for k in range(len(Pvap)):
               Pvap[k] = np.exp(A[k]-(B[k]/(C[k]+T)))
               z[k] = fk[k]/np.sum(fk)
                            
            Pcalc = 1/sum(z/Pvap) #Compute the dew point vapor pressure
            
            return(Pcalc-P) 
            
        T_new = fsolve(func,T0) #Use fsolve to iterate and find the temperature
              
        print("The dew point temperature has converged to", float(T_new), "K.")   
        return(T_new)
    
    else:
        print("Invalid Input")


#Define a function to compute the split fraction of a flash unit without knowing the feed composition
def shortcut_flash(P,T,A,B,C,epsn,n):
    '''
    

    Parameters
    ----------
    P : FLOAT
        PRESSURE IN MMHG
    T : FLOAT
        TEMPERATURE IN K
    A : ARRAY
        ANTOINE EQUATION PARAMETER A.
    B : ARRAY
        ANTOINE EQUATION PARAMETER B.
    C : ARRAY
        ANTOINE EQUATION PARAMETER C.
    epsn : FLOAT
        SPLIT FRACTION OF KEY COMPONENT.
    n : INTEGER
        INDEX OF KEY COMPONENT.

    Returns
    -------
    eps : ARRAY
        SPLIT FRACTIONS OF EACH COMPONENT.

    '''
    #Initialize empty vapour pressure and feed mole fraction array
    Pvap = np.zeros(np.shape(A))
    #Calculate the vapour pressure for each component
    for k in range(len(Pvap)):
        Pvap[k] = Antoine(A[k], B[k], C[k], T)
    
    K = Pvap/P 
    
    alpha = K/K[n] #Calculate relative volatility of each component
    
    eps = (alpha*epsn)/(1+(alpha-1)*epsn) #Calculate split fraction of each component
    
    return(eps)

#Define a function to compute the split fraction of a absorber unit without knowing the feed composition
def shortcut_absorber(n,P,T,A,B,C):
    '''
    
    Parameters
    ----------
    n : INTEGER
        INDEX OF KEY COMPONENT.
    P : FLOAT
        PRESSURE IN MMHG
    T : FLOAT
        TEMPERATURE IN K.
    A : ARRAY
        ANTOINE EQUATION PARAMETER A.
    B : ARRAY
        ANTOINE EQUATION PARAMETER B.
    C : ARRAY
        ANTOINE EQUATION PARAMETER C.

    Returns
    -------
    epsV : ARRAY
        SPLIT FRACTIONS OF EACH COMPONENT IN VAPOR PHASE.
    epsL : ARRAY
        SPLIT FRACTIONS OF EACH COMPONENT IN LIQUID PHASE.

    '''
    r = 0.99 #Fixed recovery of typically 0.99
    Ae = 10 #Absorber factor
    
    #Initialize an empty array for vapour pressures
    Pvap = np.zeros(np.shape(A))
        
    #Calculate the vapour pressure for each component
    for k in range(len(Pvap)):
        Pvap[k] = np.exp(A[k]-(B[k]/(C[k]+T)))
       
    K = Pvap/P #Compute relative volatility of each component
        
    alpha = K/K[n] #Compute relative volatility of each component
    
    #Calculate number of stages from Kremser equation
    N = np.log((r-Ae)/(Ae*(r-1)))/np.log(Ae)
    
    Ak = Ae/alpha
    Bk = (1-Ak**(N+1))/(1-Ak)
    
    #Compute split fractions in each phase
    epsV = 1/Bk
    epsL = 1 - epsV
    
    return(epsV,epsL)
