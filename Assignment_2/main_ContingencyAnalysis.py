# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 10:05:37 2024

@author: ZoéWieser
"""

import numpy as np
import PowerFlow_46705 as pf # import PowerFlow functions
import LoadNetworkData as lnd # load the network data to global variables
max_iter = 30 # Iteration settings
err_tol = 1e−4

# Load the Network data...
filename = "./TestSystem4SA.txt"
lnd.LoadNetworkData(filename) # makes Ybus available as lnd.Ybus etc.

#%%
###################################################################
# Part I: Study the base case and display results (with % loading)#
###################################################################

V,success,n = pf.PowerFlowNewton(lnd.Ybus,lnd.Sbus,lnd.V0,lnd.pv_index,lnd.pq_index,max_iter,err_tol)

if success: # Display results if the power flow analysis converged
pf.DisplayResults_and_loading(V,lnd)

#%%
######################################################################
# Part II: Simplifiedcontingencyanalysis(onlybranchoutages) #
######################################################################

print('*'*50)
print('* Contingency Analysis *')
print('*'*50)

for i in range(len(lnd.br_f)): # sweep over branches
    fr_ind = lnd.br_f[i]
    to_ind = lnd.br_t[i]
    br_ind = i
    Ybr_mat = lnd.br_Ymat[i]
    Ybus_mod, Yfr_mod, Yto_mod = pf.apply_contingency_to_Y_matrices(lnd.Ybus,lnd.Y_fr,lnd.Y_to,\
                                                                 fr_ind,to_ind,br_ind,Ybr_mat)
    
    str_status = '-'*63 + '\nTripping of branch {:} (bus {:} - bus {:})'.format(i+1, lnd.ind_to_bus[fr_ind],\
                                                                                lnd.ind_to_bus[to_ind])
    
    try: # try the load flow, if it fails, display message
        V,success,n = pf.PowerFlowNewton(Ybus_mod,lnd.Sbus,lnd.V0,lnd.pv_index,lnd.pq_index,\
                                         max_iter,err_tol,print_progress=False)
        
    except:
        str_ = '--> Load Flow error (Jacobian) when branch {:} (bus {:} - bus {:}) is tripped'
        str_ = str_.format(i+1, lnd.ind_to_bus[fr_ind], lnd.ind_to_bus[to_ind])
        print(str_status + ' [CONVERGENCE ISSUES!]')
        print(str_)
    
    else:
        if success: # Display results if the power flow analysis converged
            violations = pf.System_violations(V,Ybus_mod,Yfr_mod,Yto_mod,lnd)
        
            if not violations: # no violations, printstatusandmoveon
                print(str_status + '[OK!]')
            
            else: # ifviolation, displaythem
                print(str_status + '[Violations!]')
                for str_ in violations:
                    print(str_)
                    
        else: # no convergence...
            str_ = '--> No load-flow convergence when branch {:} (bus {:} - bus {:}) is tripped'.
            str_ = str_.format(i+1, lnd.ind_to_bus[fr_ind], lnd.ind_to_bus[to_ind])
            print(str_status + ' [CONVERGENCE ISSUES!]')
            print(str_)