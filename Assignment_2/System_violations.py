# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 10:49:58 2024

@author: ZoéWieser
"""

def System_violations(V,Ybus,Y_from,Y_to,lnd):
    ###
    # Inputs:
    # V = results from the load flow
    # Ybus = the bus admittance matrix used in the load flow
    # Y_from,Y_to = tha admittance matrices(modified) used to determine the branch flows
    # lnd = the LoadNetworkData object for easy access to other model data
    ###
    
    # Store variables as more convenient names
    br_f=lnd.br_f; br_t=lnd.br_t;   # from and to branch indices
    ind_to_bus=lnd.ind_to_bus;      # the ind_to_bus mapping object
    bus_to_ind=lnd.bus_to_ind;      # the bus_to_ind mapping object
    br_MVA = lnd.br_MVA             # object containing the MVA ratings of the branches
    br_id = lnd.br_id               # branch id?
    gen_MVA = lnd.Gen_MVA           # object containing the MVA ratings of the generators
    
    
    # Line flows and generators injection....
    ## S_k = V_k . (I_k)* apparent power
    ## SLD = (P + jQ)/mva_base with P, Q load powers in MW
    S_to = np.round(np.absolute(V[br_t]*(Y_to.dot(V)).conj()),3)        # the flow in the to end.. >modified
    S_from = np.round(np.absolute(V[br_f]*(Y_from.dot(V)).conj()),3)     # the flow in the from end >modified
    S_inj = V*(Ybus.dot(V)).conj()              # the injected power in the nodes >modified
    SLD=lnd.S_LD                                # the defined loads on the PQ busses (apparent power of loads [pu])
    S_gen = np.round(np.absolute(S_inj + SLD),3)                        # the generator arrays = injection (modified) + loads (unchanged) >updated
    
    
    violations = []     # empty list that will store strings describing each violation
    
       
    # 1. Check flow in all branches (both ends) and report if limits are violated
    for i in range(len(br_f)):
        if S_from[i] > br_MVA[i]:
            str_ = 'Branch flow limit (from) violated: FROM bus {0} to bus {1}'.format(ind_to_bus[br_f[i]], ind_to_bus[br_t[i]]) 
            violations.append(str_)
            
        if S_to[i] > br_MVA[i]:
            str_ = 'Branch flow limit (to) violated: from bus {0} to bus {1}'.format(ind_to_bus[br_f[i]], ind_to_bus[br_t[i]]) 
            violations.append(str_)
    
    
    # 2. Check output of all generators and see if limits are exceeded
    for i in range(len(gen_MVA)):
        
        if S_gen[i] > gen_MVA[i] and gen_MVA[i] != 0 :
            str_ = 'Generation limit violated: generator at bus {}'.format(ind_to_bus[i])
            violations.append(str_)
        
        
    # 3. Check voltages on all busses and see if it remains between 0.9 and 1.1 pu
    for i in range(len(bus_to_ind)):
        Vm = abs(V[i])     # voltage magnitude at the bus indexed i
        
        if Vm < 0.9:
            str_ = 'Bus voltage limit violated: voltage too low at bus {}'.format(ind_to_bus[i])
            violations.append(str_)
        
        if Vm > 1.1:
            str_ = 'Bus voltage limit violated: voltage too high at bus {}'.format(ind_to_bus[i])
            violations.append(str_)  
    
    
    return violations   # return the list with description of all of the violations

