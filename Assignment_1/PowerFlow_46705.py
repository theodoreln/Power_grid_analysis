"""
46705 - Power Grid Analysis
This file contains the definitions of the functions needed to
carry out Power Flow calculations in python.

How to carry out Power Flow in a new *.py file? 
See the example in table 1 in the assignment text
"""

import numpy as np
from tabulate import tabulate

# 1. the PowerFlowNewton() function
def PowerFlowNewton(Ybus,Sbus,V0,pv_index,pq_index,max_iter,err_tol):
    success = 0 #Initialization of status flag and iteration counter
    n = 0
    V = V0
    print(' iteration maximum P & Q mismatch (pu)')
    print(' --------- ---------------------------')
    # Determine mismatch between initial guess and and specified value for P and Q
    F = calculate_F( Ybus , Sbus , V , pv_index , pq_index)
    # Check if the desired tolerance is reached
    success = CheckTolerance ( F , n , err_tol)
    # Start the Newton iteration loop
    while (not success) and (n < max_iter) :
        n += 1 # Update counter
        # Compute derivatives and generate the Jacobian matrix
        J_dS_dVm , J_dS_dTheta = generate_Derivatives ( Ybus , V)
        J = generate_Jacobian ( J_dS_dVm , J_dS_dTheta , pv_index , pq_index)
        # Compute the update step
        dx = np . linalg . solve(J , F)
        # Update voltages and check if tolerance is now reached
        V = Update_Voltages ( dx , V , pv_index , pq_index)
        F = calculate_F( Ybus , Sbus , V , pv_index , pq_index)
        success = CheckTolerance ( F , n , err_tol)
    
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



####################################################
#  Displaying the results in the terminal window   #
####################################################

''' I'm not 100% sure if the calculations are correct. 
 But when we have extracted the data, we can compare the results.

 And the assignment doesn't say that bus_label is an input. 
 I wouldn't know how to get the information otherwise, though. So I added it.'''

# Function that displays the results of the Power Flow, 
def DisplayResults(V,lnd):
    
    Ybus=lnd.Ybus ; Y_from=lnd.Y_fr ; Y_to=lnd.Y_to ; br_f=lnd.br_f ; br_t=lnd.br_t ;
    buscode=lnd.buscode; SLD=lnd.S_LD ; ind_to_bus=lnd.ind_to_bus;
    bus_to_ind=lnd.bus_to_ind ; MVA_base=lnd.MVA_base ; bus_labels=lnd.bus_labels ;
    ref=lnd.ref

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

    for s in S_inj:
        if np.real(s)>=0: 
            P_inj_g = np.concatenate((P_inj_g, [str(np.round(np.real(s),3))]))
            P_inj_l = np.concatenate((P_inj_l, ["-"]))
            Q_inj_g = np.concatenate((Q_inj_g, [str(np.round(np.imag(s),3))]))
            Q_inj_l = np.concatenate((Q_inj_l, ["-"]))
        else:
            P_inj_l = np.concatenate((P_inj_l, [str(np.round(-np.real(s),3))]))
            P_inj_g = np.concatenate((P_inj_g, ["-"]))
            Q_inj_l = np.concatenate((Q_inj_l, [str(np.round(-np.imag(s),3))]))
            Q_inj_g = np.concatenate((Q_inj_g, ["-"]))


    # Branches
    n_br = len(br_f) # Number of Branches
    branch_no = np.arange(1,(n_br+1)) # Branch numbers

    S_to = V[br_t]*Y_to.dot(V).conj() # Apparent Power flowing into the receiving end busses
    S_from = V[br_f]*Y_from.dot(V).conj() # Apparent Power flowing from the receiving end busses
    
    # Having number of branches and not there index
    mapping_function = np.vectorize(lambda x: ind_to_bus[x])
    br_from = mapping_function(br_f)
    br_to = mapping_function(br_t)


    bus_data = np.vstack((bus_no, bus_labels, Vm, Theta, Sm, P_inj_g, Q_inj_g, P_inj_l, Q_inj_l))
    bus_data = bus_data.transpose()


    branch_data = np.vstack((branch_no.astype(str), br_from, br_to, np.round(np.real(S_from),3), np.round(np.imag(S_from),3), np.round(np.real(S_to),3), np.round(np.imag(S_to),3)))
    branch_data = branch_data.transpose()


    # Display Bus results
    print("================================================================================")
    print("                           | Bus results |                             ")
    print("================================================================================")
    print("Bus    Bus            Voltage          Power     Generation          Load   ")

    print(tabulate(bus_data, headers=["#", "Label", "Mag(pu)", "Ang(deg)", "S(pu)","P(pu)","Q(pu)","P(pu)","Q(pu)"],stralign="center"))
    print("\n")
    # Display Branch flow
    print("=======================================================")
    print("                   | Branch flow |                     ")
    print("=======================================================")
    print("Branch  From    To    From Bus Inject.   To Bus Inject.  ")
    print(tabulate(branch_data, headers=[" # ", "Bus", "Bus", "P(pu)", "Q(pu)", "P(pu)", "Q(pu)"],stralign="center"))





""" 
To try Display Results:
V= np.array([1,1,1])
Ybus=np.array([[1+1j, -2+2j, 1-1j],
               [-2+2j, 2+2j, 2+2j],
               [1-1j, 2+2j, -1-2j]])
Y_from= np.array([1+1j,2-1j,2+1j])
Y_to= np.array([-1+1j,2+1j,-2-1j])
br_f = np.array([1,1,2])-1
br_t = np.array([2,3,3])-1
buscode = np.array([3,1,2])
bus_label = np.array(["BUS1HV","BUS2HV","BUS3HV"])

DisplayResults(V,Ybus,Y_from,Y_to,br_f,br_t,buscode, bus_label)"""
