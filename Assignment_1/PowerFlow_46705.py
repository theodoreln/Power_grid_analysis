"""
46705 - Power Grid Analysis
This file contains the definitions of the functions needed to
carry out Power Flow calculations in python.

How to carry out Power Flow in a new *.py file? 
See the example in table 1 in the assignment text
"""

import numpy as np


# 1. the PowerFlowNewton() function
def PowerFlowNewton(Ybus,Sbus,V0,pv_index,pq_index,max_iter,err_tol):
    ''' String here with purpose '''
    # implement your code here

    return V,success,n


# 2. the calculate_F() function
def calculate_F(Ybus,Sbus,V,pv_index,pq_index):

    return F


# 3. the CheckTolerance() function
def CheckTolerance(F,n,err_tol):

    return success

# 4. the generate_Derivatives() function
def generate_Derivatives(Ybus,V):

    return J_ds_dVm,J_dS_dTheta


# 5. the generate_Jacobian() function
def generate_Jacobian(J_dS_dVm,J_dS_dTheta,pv_index,pq_index):

    return J


# 6. the Update_Voltages() function
def Update_Voltages(dx,V,pv_index,pq_index):

    return V



####################################################
#  Displaying the results in the terminal window   #
####################################################
def DisplayResults(V,lnd):

    return
