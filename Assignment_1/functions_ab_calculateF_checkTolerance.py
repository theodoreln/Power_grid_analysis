# -*- coding: utf-8 -*-
"""
Created on Sun Feb 11 19:04:56 2024

@author: Zoe Wieser
"""

import numpy as np

# 2. the calculate_F() function
# Calculates the mismatch between specified values of Psp and Qsp, and the calculated values of P and Q
def calculate_F(Ybus,Sbus,V,pv_index,pq_index):
    # =============================================================================
    # Ybus: bus admittance matrix (N×N) 
    # Sbus: specified apparent power injection vector (N×1)
    # V: complex bus voltage vector (N×1) 
    # pv_index: indice for PV-busses
    # pq_index: indice for PQ-busses
    # =============================================================================
    
    Delta_S = Sbus - V * (Ybus.dot(V)).conj() 
    
    # real_index = np.isreal(Delta_S)
    # print(real_index)
    # imag_index = np.iscomplex(Delta_S)
    # print(imag_index)    
    
    Delta_P = Delta_S.real
    Delta_Q = Delta_S.imag
    
    # print(Delta_P,Delta_Q)
    # print('/n''to compare:''/n')
    # print(Delta_S[real_index],Delta_S[imag_index])
    
    F= np.concatenate((Delta_P[pv_index],Delta_P[pq_index],Delta_Q[pq_index]),axis=0)
    
    return F


# 3. the CheckTolerance() function
# Checks whether the greatest mismatch in the mismatch vector is smaller than the specified tolerance (return 1 if that's the case, 0 otherwise). 
# Displays as well the absolute value of the greatest mismatch as well as the iteration number.
def CheckTolerance(F,n,err_tol):
    # =============================================================================
    # F: mismatch vector
    # n: iteration counter
    # err_tol: specified error tolerance
    # =============================================================================
    
    normF = np.linalg.norm(F, np.inf) # returns the infinite norm of the vector F: max(abs(F))
    
    success = (normF <= err_tol ) # returns 1 if normF <= specified tolerance ; 0 otherwise
    
    return success