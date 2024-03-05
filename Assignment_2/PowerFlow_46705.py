# -*- coding: utf-8 -*-
"""
46705 - Power Grid Analysis - Assignment 2
This file contains the definitions of the functions needed to
carry out Contingency Analysis calculations in python.
"""

import numpy as np
from tabulate import tabulate


# Function that displays the results of the Power Flow, 
def DisplayResults_and_loading(V,lnd):
    
    Ybus=lnd.Ybus ; Y_from=lnd.Y_fr ; Y_to=lnd.Y_to ; br_f=lnd.br_f ; br_t=lnd.br_t ;
    buscode=lnd.buscode; SLD=lnd.S_LD ; ind_to_bus=lnd.ind_to_bus;
    bus_to_ind=lnd.bus_to_ind ; MVA_base=lnd.MVA_base ; bus_labels=lnd.bus_labels ;
    ref=lnd.ref; Gen_MVA = lnd.Gen_MVA; br_MVA = lnd.br_MVA

    # Busses
    N = len(V) # Number of busses
    bus_no = [ind_to_bus[i] for i in range(N)] # Bus numbers
    bus_no[ref[0]] = "*"+str(bus_no[ref[0]])+"*"   # Highliting the reference bus

    Vm = np.round(np.absolute(V),3) # Voltage magnitude
    Theta = np.round(np.angle(V)*180/np.pi,2) #Voltage angle

    S_inj = V*(Ybus.dot(V)).conj() # Generation and load power
    Sm = np.round(np.absolute(S_inj),2)

    P_inj_g = np.array([]) #Generation active Power
    P_inj_l = np.array([]) #Load active Power
    Q_inj_g = np.array([]) #Generation reactive Power
    Q_inj_l = np.array([]) #Load reactive Power
    loading_inj = np.array([]) #Loading of generator at bus
    S_gen = S_inj-SLD
    print(S_inj, S_gen, SLD)

    for index,s in enumerate(S_gen):
        if np.real(s)>0: 
            P_inj_g = np.concatenate((P_inj_g, [str(np.round(np.real(s),3))]))
            Q_inj_g = np.concatenate((Q_inj_g, [str(np.round(np.imag(s),3))]))
            loading_inj = np.concatenate((loading_inj,["{0}%".format(str(np.round(Sm[index]*100/Gen_MVA[index],2)))]))
        else:
            P_inj_g = np.concatenate((P_inj_g, ["-"]))
            Q_inj_g = np.concatenate((Q_inj_g, ["-"]))
            loading_inj = np.concatenate((loading_inj,["-"]))

        if SLD[index] != 0:
            P_inj_l = np.concatenate((P_inj_l, [str(np.round(np.real(SLD[index]),3))]))
            Q_inj_l = np.concatenate((Q_inj_l, [str(np.round(np.imag(SLD[index]),3))]))
        else:
            P_inj_l = np.concatenate((P_inj_l, ["-"]))
            Q_inj_l = np.concatenate((Q_inj_l, ["-"]))
        

    # Branches
    n_br = len(br_f) # Number of Branches
    branch_no = np.arange(1,(n_br+1)) # Branch numbers

    S_to = V[br_t]*Y_to.dot(V).conj() # Apparent Power flowing into the receiving end busses
    loading_to= np.round(S_to*100/br_MVA,3) # Loading of line at to branch
    S_from = V[br_f]*Y_from.dot(V).conj() # Apparent Power flowing from the receiving end busses
    loading_from= np.round(S_to*100/br_MVA,3) # Loading of line at from branch
    # Having number of branches and not there index
    mapping_function = np.vectorize(lambda x: ind_to_bus[x])
    br_from = mapping_function(br_f)
    br_to = mapping_function(br_t)


    bus_data = np.vstack((bus_no, bus_labels, Vm, Theta, P_inj_g, Q_inj_g, loading_inj, P_inj_l, Q_inj_l))
    bus_data = bus_data.transpose()


    branch_data = np.vstack((branch_no.astype(str), br_from, br_to, np.round(np.real(S_from),3), np.round(np.imag(S_from),3),loading_from, np.round(np.real(S_to),3), np.round(np.imag(S_to),3), loading_to))
    branch_data = branch_data.transpose()


    # Display Bus results
    print("================================================================================")
    print("                           | Bus results |                             ")
    print("================================================================================")
    print("Bus    Bus            Voltage                Generation                Load   ")

    print(tabulate(bus_data, headers=["#", "Label", "Mag(pu)", "Ang(deg)","P(pu)","Q(pu)","loading","P(pu)","Q(pu)"],stralign="center"))
    print("\n")
    # Display Branch flow
    print("=============================================================")
    print("                      | Branch flow |                        ")
    print("=============================================================")
    print("Branch  From    To       From Bus Inject.      To Bus Inject.  ")
    print(tabulate(branch_data, headers=[" # ", "Bus", "Bus", "P(pu)", "Q(pu)", "loading", "P(pu)", "Q(pu)", "loading"],stralign="center"))



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
    S_to = V[br_t]*(Y_to.dot(V)).conj()         # the flow in the to end.. >modified
    S_from = V[br_f]*(Y_from.dot(V)).conj()     # the flow in the from end >modified
    S_inj = V*(Ybus.dot(V)).conj()              # the injected power in the nodes >modified
    SLD=lnd.S_LD                                # the defined loads on the PQ busses (apparent power of loads [pu])
    S_gen = S_inj + SLD                         # the generator arrays = injection (modified) + loads (unchanged) >updated
    
    
    violations = []     # empty list that will store strings describing each violation
    
       
    # 1. Check flow in all branches (both ends) and report if limits are violated
    for i in range(len(br_f)):
        if S_from > br_MVA[i]:
            str_ = 'branch flow limit (from) violated: FROM bus {0} to bus {1}'.format(ind_to_bus[br_f[i]], ind_to_bus[br_t[i]]) 
            violations += str_
            
        if S_to > br_MVA[i]:
            str_ = 'branch flow limit (to) violated: from bus {0} TO bus {1}'.format(ind_to_bus[br_f[i]], ind_to_bus[br_t[i]]) 
            violations += str_
    
    
    # 2. Check output of all generators and see if limits are exceeded
    for i in range(len(gen_MVA)):
        P = S_gen.real[i]   # active power of the generator at bus indexed i
        Q = S_gen.imag[i]   # reactive power of the generator at bus indexed i
        
        if P > gen_MVA[i]:
            str_ = 'generation limit violated: active power at bus {}'.format(ind_to_bus[i])
            violations += str_
        
        if Q > gen_MVA[i]:
            str_ = 'generation limit violated: reactive power at bus {}'.format(ind_to_bus[i])  
            violations += str_
        
        
    # 3. Check voltages on all busses and see if it remains between 0.9 and 1.1 pu
    for i in range(len(bus_to_ind)):
        Vm = abs(V[i])     # voltage magnitude at the bus indexed i
        
        if Vm < 0.9:
            str_ = 'bus voltage limit violated: voltage too low at bus {}'.format(ind_to_bus[i])
            violations += str_
        
        if Vm > 1.1:
            str_ = 'bus voltage limit violated: voltage too high at bus {}'.format(ind_to_bus[i])
            violations += str_    
    
    
    return violations   # return the list with description of all of the violations



def apply_contingency_to_Y_matrices(Ybus,Yfr,Yto,ind_fr,ind_to,br_ind,Ybr_mat):
    # input:
    # The original admittance matrices: Ybus,Yfr,Yto
    # The from and to end indices for the branch (fr_ind, to_ind)
    # The indice for where the branch is in the branch list (br_ind)
    # The 2x2 admittance matrix for the branch Ybr_mat
    ##########################################################
    # cpoying the original matrices
    Ybus_mod = Ybus.copy() 
    Yfr_mod = Yfr.copy()
    Yto_mod = Yto.copy() 
    #deleting the entry in the Ybus
    Ybus_mod[ind_fr, ind_fr] -= Ybr_mat[0,0]
    Ybus_mod[ind_fr, ind_to] -= Ybr_mat[0,1]
    Ybus_mod[ind_to, ind_fr] -= Ybr_mat[1,0]
    Ybus_mod[ind_to, ind_to] -= Ybr_mat[1,1]
    #deleting the row fr_ind of the branch
    Yfr_mod[br_ind, ind_to] = 0
    Yfr_mod[br_ind, ind_fr] = 0
    #deleting the row to_ind of the branch
    Yto_mod[br_ind, ind_to] = 0
    Yto_mod[br_ind, ind_fr] = 0
   
    return Ybus_mod,Yfr_mod,Yto_mod




"""
46705 - Power Grid Analysis - Assignment 1
This following part of the file contains the definitions of the functions 
needed to carry out Power Flow calculations in python.
"""

# 1. the PowerFlowNewton() function
def PowerFlowNewton(Ybus,Sbus,V0,pv_index,pq_index,max_iter,err_tol):
    success = 0 #Initialization of status flag and iteration counter
    n = 0
    V = V0
    print(' iteration maximum P & Q mismatch (pu)')
    print(' --------- ---------------------------')
    
    # Determine mismatch between initial guess and and specified value for P and Q
    F = calculate_F(Ybus, Sbus, V, pv_index, pq_index)
    
    # Check if the desired tolerance is reached
    success = CheckTolerance(F ,n ,err_tol)
    
    # Start the Newton iteration loop
    while (not success) and (n < max_iter):
        n += 1 # Update counter
        
        # Compute derivatives and generate the Jacobian matrix
        J_dS_dVm, J_dS_dTheta = generate_Derivatives (Ybus, V)
        J = generate_Jacobian (J_dS_dVm ,J_dS_dTheta ,pv_index ,pq_index)
        
        # Compute the update step
        dx = np.linalg.solve(J ,F)
        
        # Update voltages and check if tolerance is now reached
        V = Update_Voltages(dx ,V, pv_index , pq_index)
        F = calculate_F( Ybus , Sbus , V , pv_index , pq_index)
        success = CheckTolerance( F , n , err_tol)
    
    if success : #print out message concerning wether the power flow converged or not
        print('The Newton Rapson Power Flow Converged in %d iterations!' % (n , ) )
    else :
        print('No Convergence !!!\n Stopped after %d iterations without solution...' % (n , ) )

    return V,success,n


# 2. the calculate_F() function
# Calculates the mismatch between specified values of Psp and Qsp, and the calculated values of P and Q
def calculate_F(Ybus,Sbus,V,pv_index,pq_index):
    Delta_S = Sbus - V * (Ybus.dot(V)).conj() 

    Delta_P = Delta_S.real
    Delta_Q = Delta_S.imag
    
    F= np.concatenate((Delta_P[pv_index],Delta_P[pq_index],Delta_Q[pq_index]),axis=0)
    
    return F


# 3. the CheckTolerance() function
# Checks whether the greatest mismatch in the mismatch vector is smaller than the specified tolerance (return 1 if that's the case, 0 otherwise). 
# Displays as well the absolute value of the greatest mismatch as well as the iteration number.
def CheckTolerance(F,n,err_tol):
    normF = np.linalg.norm(F, np.inf) # returns the infinite norm of the vector F: max(abs(F))
    
    success = (normF <= err_tol ) # returns 1 if normF <= specified tolerance ; 0 otherwise
    
    return success


# 4. the generate_Derivatives() function
def generate_Derivatives(Ybus,V):
    J_ds_dVm = np.diag(V/np.absolute(V)).dot(np.diag((Ybus.dot(V)).conj())) + np.diag(V).dot(Ybus.dot(np.diag(V/np.absolute(V))).conj())
                
    J_dS_dTheta = 1j*np.diag(V).dot((np.diag(Ybus.dot(V))-Ybus.dot(np.diag(V))).conj())

    return J_ds_dVm,J_dS_dTheta


# 5. the generate_Jacobian() function
def generate_Jacobian(J_dS_dVm,J_dS_dTheta,pv_index,pq_index):
    pvpq_ind = np.append(pv_index, pq_index)
    J_11 = np.real(J_dS_dTheta[np.ix_(pvpq_ind, pvpq_ind)])
    J_12 = np.real(J_dS_dVm[np.ix_(pvpq_ind, pq_index)])
    J_21 = np.imag(J_dS_dTheta[np.ix_(pq_index, pvpq_ind)])
    J_22 = np.imag(J_dS_dVm[np.ix_(pq_index, pq_index)])
    J = np.block([[J_11,J_12],[J_21,J_22]])

    return J


# 6. the Update_Voltages() function
def Update_Voltages(dx,V,pv_index,pq_index):
    N1 = 0
    N2 = len(pv_index) 
    N3 = N2
    N4 = N3 + len(pq_index) 
    N5 = N4
    N6 = N5 + len(pq_index)
    
    Theta = np.angle(V) 
    Vm = np.absolute(V)    
    
    if len(pv_index)>0:
        Theta[pv_index] += dx[N1:N2]   
            
    if len(pq_index)>0:
        Theta[pq_index] += dx[N3:N4]
        Vm[pq_index] += dx[N5:N6]

    V = Vm * np.exp(1j*Theta)   
    return V
