# -*- coding: utf-8 -*-
"""
Created on Tue Feb 13 09:32:14 2024

@author: group12_46705
"""

import PowerFlow_46705 as pf # importPower Flowfunctions
import LoadNetworkData as lnd #load the network data to global variables

max_iter = 30 #Iteration settings
err_tol = 10**(-4)

# Load the Network data...
# filename = 'TestSystem.txt'
filename = 'Kundur_two_area_system.txt'
lnd.LoadNetworkData(filename) # makes Ybus available as lnd.Ybus and etc.

# Carry out the power flow analysis...
V,success,n = pf.PowerFlowNewton(lnd.Ybus,lnd.Sbus,lnd.V0,lnd.pv_index,lnd.pq_index,max_iter,err_tol)

# Display results if the power flow analysis converged
if success:
    print('Success')
    pf.DisplayResults(V,lnd) #Now we are just passing the lnd module as input