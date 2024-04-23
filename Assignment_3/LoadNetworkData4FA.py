# -*- coding: utf-8 -*-
import numpy as np
import ReadNetworkData4FA as rd4fa

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
    Ybus = np.zeros((N,N),dtype=complex)
    Ybus0 = np.zeros((N,N),dtype=complex)
    Ybus2 = np.zeros((N,N),dtype=complex)
    # Continue with your code here...
    # ...
    # End your code by deriving the Zbus matrices (Zbus0, Zbus1, Zbus2)  
    Zbus0 = np.linalg.inv(Ybus0)
    Zbus1 = np.linalg.inv(Ybus)
    Zbus2 = np.linalg.inv(Ybus2)
    
    return   