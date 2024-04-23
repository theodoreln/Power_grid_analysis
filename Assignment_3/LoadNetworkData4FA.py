# -*- coding: utf-8 -*-
import numpy as np
import ReadNetworkData4FA as rd4fa

# filename = "./TestSystem4FA.txt"
# bus_data,load_data,gen_data,line_data,tran_data,mva_base,bus_to_ind,ind_to_bus = rd4fa.read_network_data_4fa_from_file(filename)

def LoadNetworkData4FA(filename):
    global Ybus,Sbus,V0,buscode,pq_index,pv_index,Y_fr,Y_to,br_f,br_t,br_Y,S_LD, \
           ind_to_bus,bus_to_ind,MVA_base,bus_labels,Ybus0,Ybus2,Zbus0,Zbus1,Zbus2          
    # read in the data from the file...
    bus_data,load_data,gen_data,line_data,tran_data,mva_base,bus_to_ind,ind_to_bus = \
    rd4fa.read_network_data_4fa_from_file(filename)

    ############################################################################################## 
    # Construct the Ybus (positive-sequence), Ybus0 (zero-sequence), and Ybus2 (negative-sequence)
    # matrices from elements in the line_data and trans_data
    # Keep/modify code from the Python power flow program as needed
    ##########################################################################  
    N = len(bus_data) # Number of buses
    Ybus0 = np.zeros((N,N),dtype=complex)
    Ybus1 = np.zeros((N,N),dtype=complex)
    Ybus2 = np.zeros((N,N),dtype=complex)
    
    # Information about the line
    for line in line_data:
        # Taking the data
        bus_from, bus_to, id_, R, X, B, X2, X0 = line
        id_from, id_to = bus_to_ind[bus_from], bus_to_ind[bus_to]
        # Filling the Ybus matrix
        Y0 = 1/(1j*X0)
        Ybus0[id_from,id_from] += Y0
        Ybus0[id_to, id_to] += Y0
        Ybus0[id_from, id_to] += -Y0
        Ybus0[id_to, id_from] += -Y0
        Y1 = 1/(1j*X)
        Ybus1[id_from,id_from] += Y1
        Ybus1[id_to, id_to] += Y1
        Ybus1[id_from, id_to] += -Y1
        Ybus1[id_to, id_from] += -Y1
        Y2 = 1/(1j*X2)
        Ybus2[id_from,id_from] += Y2
        Ybus2[id_to, id_to] += Y2
        Ybus2[id_from, id_to] += -Y2
        Ybus2[id_to, id_from] += -Y2
        
    # Information about the generator
    for gen in gen_data:
        # Taking the data
        bus_nr, MVA_size, P_gen, X, X2, X0, Xn, GRND = gen
        id_bus = bus_to_ind[bus_nr]
        # Filling the Ybus matrix
        if GRND == 1:
            Y0 = 1/(1j*(X0+3*Xn))
        elif GRND == 0 :
            Y0 = 1/(1j*X0)
        Ybus0[id_bus,id_bus] += Y0
        Y1 = 1/(1j*X)
        Ybus1[id_bus,id_bus] += Y1
        Y2 = 1/(1j*X2)
        Ybus2[id_bus,id_bus] += Y2
        
    # Information about the transformer
    for tran in tran_data:
        # Taking the data
        bus_from, bus_to, id_, R, X, n, ang_deg, fr_co, to_co, X2, X0 = tran
        id_from, id_to = bus_to_ind[bus_from], bus_to_ind[bus_to]
        # Filling the Ybus matrix
        if fr_co == 2 and to_co == 2 :
            Y0 = 1/(1j*X0)
            Ybus0[id_from,id_from] += Y0
            Ybus0[id_to, id_to] += Y0
            Ybus0[id_from, id_to] += -Y0
            Ybus0[id_to, id_from] += -Y0
        if fr_co == 3 and to_co == 2 :
            Y0 = 1/(1j*X0)
            Ybus0[id_to, id_to] += Y0
        if fr_co == 2 and to_co == 3 :
            Y0 = 1/(1j*X0)
            Ybus0[id_from, id_from] += Y0
        Y1 = 1/(1j*X)
        Ybus1[id_from,id_from] += Y1
        Ybus1[id_to, id_to] += Y1
        Ybus1[id_from, id_to] += -Y1
        Ybus1[id_to, id_from] += -Y1
        Y2 = 1/(1j*X2)
        Ybus2[id_from,id_from] += Y2
        Ybus2[id_to, id_to] += Y2
        Ybus2[id_from, id_to] += -Y2
        Ybus2[id_to, id_from] += -Y2
    
    # End your code by deriving the Zbus matrices (Zbus0, Zbus1, Zbus2)  
    Zbus0 = np.linalg.inv(Ybus0)
    Zbus1 = np.linalg.inv(Ybus1)
    Zbus2 = np.linalg.inv(Ybus2)
    
    return   