import numpy as np 
import cmath
import math
import ReadNetworkData as rd

# bus_data , load_data , gen_data , line_data , tran_data , mva_base , bus_to_ind , ind_to_bus = rd.read_network_data_from_file('TestSystem.txt')

def LoadNetworkData(filename):
    # Create global variables
    global Ybus , Sbus , V0 , buscode , bus_kv , pq_index , pv_index , ref , Y_fr , Y_to , br_f , br_t , br_Y , S_LD , ind_to_bus , bus_to_ind , MVA_base , bus_labels
    # Read in the data from the file
    bus_data , load_data , gen_data , line_data , tran_data , mva_base , bus_to_ind , ind_to_bus = \
    rd . read_network_data_from_file(filename)
    
    # Informations about the data
    N = len(bus_data) #Number of buses
    M_lines = len(line_data) #Number of lines
    M_trans = len(tran_data) #Number of transformers
    M_branches = M_lines + M_trans 
    MVA_base = mva_base
    
    # Get the index of the buses 
    bus_labels = np.array([bus_data[i][1] for i in range(N)])
    bus_kv = np.array([bus_data[i][2] for i in range(N)])
    buscode = np.array([bus_data[i][3] for i in range(N)])
    pq_index = np.where(buscode == 1)[0] # Find indices for all PQ-busses
    pv_index = np.where(buscode == 2)[0] # Find indices for all PV-busses
    ref = np.where(buscode == 3)[0] # Find index for ref bus
    
    # Ybus computation
    Ybus = np.zeros((N,N),dtype=complex)
    Y_fr = np.zeros((M_branches,N), dtype=complex)
    Y_to = np.zeros((M_branches,N), dtype=complex)
    br_f = np.zeros(M_branches, dtype=int)
    br_t = np.zeros(M_branches, dtype=int)
    for i in range(M_lines):
        # Taking the data
        id_from, id_to, R, X, B = bus_to_ind[line_data[i][0]], bus_to_ind[line_data[i][1]], line_data[i][3], line_data[i][4], line_data[i][5]
        Z = R + 1j*X
        Y = 1/Z
        b = 1j*B/2
        # Filling the Ybus matrix
        Ybus[id_from,id_from] += Y+b
        Ybus[id_to, id_to] += Y+b
        Ybus[id_from, id_to] += -Y
        Ybus[id_to, id_from] += -Y
        # Filling the Y_fr and Y_to matrix
        Y_fr[i, id_from] += Y+b
        Y_fr[i, id_to] += -Y
        Y_to[i, id_from] += -Y
        Y_to[i, id_to] += Y+b
        # Indices of the from and to of lines
        br_f[i] = id_from
        br_t[i] = id_to

    for i in range(M_trans) :
        # Taking the data
        id_from, id_to, R, X, n, ang_deg = bus_to_ind[tran_data[i][0]], bus_to_ind[tran_data[i][1]], tran_data[i][3], tran_data[i][4], tran_data[i][5], tran_data[i][6]
        Z = R + 1j*X
        Y = 1/Z
        a,b = cmath.rect(n, ang_deg*math.pi/180).real,  cmath.rect(n, ang_deg*math.pi/180).imag
        # Filling the Ybus matrix
        Ybus[id_from, id_from] += Y/(a*a+b*b)
        Ybus[id_from, id_to] += -Y/(a-1j*b)
        Ybus[id_to, id_from] += -Y/(a+1j*b)
        Ybus[id_to, id_to] += Y 
        # Filling the Y_fr and Y_to matrix
        Y_fr[M_lines+i, id_from] += Y/(a*a+b*b)
        Y_fr[M_lines+i, id_to] += -Y/(a-1j*b)
        Y_to[M_lines+i, id_from] += -Y/(a+1j*b)
        Y_to[M_lines+i, id_to] += Y
        # Indices of the from and to of trans
        br_f[M_lines+i] = id_from
        br_t[M_lines+i] = id_to

    
    # Apparent power computation 
    S_LD = np.array([0]*N, dtype=complex)
    Sbus = np.array([0]*N, dtype=complex)
    for i in range(len(gen_data)):
        id_bus, mva_size, p_gen = bus_to_ind[gen_data[i][0]], gen_data[i][1], gen_data[i][2]
        if id_bus not in ref :
            Sbus[id_bus] = p_gen/MVA_base
    for i in range(len(load_data)):
        id_bus, p, q = bus_to_ind[load_data[i][0]], load_data[i][1], load_data[i][2]
        S_LD[id_bus] = (p+1j*q)/MVA_base
        Sbus[id_bus] = -(p+1j*q)/MVA_base
    
    # Initial guess for V vector
    V0 = np.ones(N, dtype=complex)

    return

# LoadNetworkData('TestSystem.txt')
