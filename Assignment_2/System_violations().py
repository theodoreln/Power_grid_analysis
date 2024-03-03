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
    # Y_from,Y_to = tha admittance matrices used to determine the branch flows
    # lnd = the LoadNetworkData object for easy access to other model data
    ###
    
    # Store variables as more convenient names
    br_f=lnd.br_f; br_t=lnd.br_t;   # from and to branch indices
    ind_to_bus=lnd.ind_to_bus;      # the ind_to_bus mapping object
    bus_to_ind=lnd.bus_to_ind;      # the bus_to_ind mapping object
    br_MVA = lnd.br_MVA             # object containing the MVA ratings of the branches
    br_id = lnd.br_id               # (you have to update LoadNetworkData for this)
    
    
    # Line flows and generators injection....
    S_to = V[br_t]*(Y_to.dot(V)).conj()         # the flow in the to end..
    S_from = V[br_f]*(Y_from.dot(V)).conj()     # the flow in the from end
    S_inj = V*(Ybus.dot(V)).conj()              # the injected power
    SLD=lnd.S_LD                                # the defined loads on the PQ busses
    S_gen = S_inj + SLD                         # the generator arrays
    
    
    violations = []     # empty list that will store strings describing each violation
    
    
    ##################################################################################
    #
    # YOUR CODE COMES HERE:
    # 1. Check flow in all branches (both ends) and report if limits are violated
    # 2. Check output of all generators and see if limits are exceeded
    # 3. Check voltages on all busses and see if it remains between 0.9 and 1.1 pu
    #
    ####################################################################################
    
    return violations   # return the list with description of all of the violations