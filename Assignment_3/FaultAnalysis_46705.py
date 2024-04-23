"""
46705 - Power Grid Analysis
This file contains the definitions of the functions needed to
carry out Fault Analysis calculations in python.
"""

import numpy as np

# 1. the FaultAnalysis() function
def FaultAnalysis(Zbus0,Zbus1,Zbus2,bus_to_ind,fault_bus,fault_type,Zf,Vf):
    ''' String here with purpose '''
    # calculate sequence fault currents
    Iseq = Calculate_Sequence_Fault_Currents(Zbus0,Zbus1,Zbus2,bus_to_ind,fault_bus,fault_type,Zf,Vf)
    # calculate sequence fault voltages
    Vseq_mat = Calculate_Sequence_Fault_Voltages(Zbus0,Zbus1,Zbus2,bus_to_ind,fault_bus,Vf,Iseq)
    # convert sequence currents to phase (fault) currents
    Iph = Convert_Sequence2Phase_Currents(Iseq)
    # convert sequence voltages to phase line-to-ground (fault) voltages
    Vph_mat = Convert_Sequence2Phase_Voltages(Vseq_mat)    
    return Iph, Vph_mat, Iseq, Vseq_mat

# 1.1. the Calculate_Sequence_Fault_Currents() function
def Calculate_Sequence_Fault_Currents(Zbus0,Zbus1,Zbus2,bus_to_ind,fault_bus,fault_type,Zf,Vf):
#  fault_type: 0 = 3-phase balanced fault; 1 = Single Line-to-Ground fault;
#              2 = Line-to-Line fault;     3 = Double Line-to-Ground fault.    
    # Iseq current array: 
    # Iseq[0] = zero-sequence; Iseq[1] = positive-sequence; Iseq[2] = negative-sequence
    Iseq = np.zeros(3,dtype=complex)
    fb = bus_to_ind[fault_bus]
    
    # Find the Thevenin equivalent 
    Zth0 = Zbus0[fb,fb]
    Zth1 = Zbus1[fb,fb]
    Zth2 = Zbus2[fb,fb]
    
    if fault_type == 0:
        Iseq[0] = 0
        Iseq[1] = Vf / Zth1
        Iseq[2] = 0
    elif fault_type == 1:
        Iseq[0] = Vf / (Zth0 + Zth1 + Zth2 + 3*Zf)
        Iseq[1] = Vf / (Zth0 + Zth1 + Zth2 + 3*Zf)
        Iseq[2] = Vf / (Zth0 + Zth1 + Zth2 + 3*Zf)
    elif fault_type == 2:
        Iseq[0] = 0
        Iseq[1] = Vf / (Zth1 + Zth2 + Zf)
        Iseq[2] = -Vf / (Zth1 + Zth2 + Zf)
    elif fault_type == 3:
        Iseq[1] = Vf / (Zth1 + (1 / (1/Zth2 + 1/(Zth0 + 3*Zf))))
        Iseq[0] = -Iseq[1] * Zth2 / (Zth0 + 3*Zf + Zth2)
        Iseq[2] = -Iseq[1] * (Zth0 + 3*Zf) / (Zth0 + 3*Zf + Zth2)
    else:
        print('Unknown Fault Type')
    return Iseq

# 1.2 the Calculate_Sequence_Fault_Voltages() function
def Calculate_Sequence_Fault_Voltages(Zbus0,Zbus1,Zbus2,bus_to_ind,fault_bus,Vf,Iseq):
    N = len(bus_to_ind)
    fb = bus_to_ind[fault_bus]
    Vseq_mat = np.zeros((N,3),dtype=complex)
    
    for i in range(N) :
        Vseq_mat[i,0] = - Zbus0[i,fb]*Iseq[0]
        Vseq_mat[i,1] = Vf - Zbus1[i,fb]*Iseq[1]
        Vseq_mat[i,2] = - Zbus2[i,fb]*Iseq[2]
    return Vseq_mat

# 1.3. the Convert_Sequence2Phase_Currents() function
def Convert_Sequence2Phase_Currents(Iseq):
    a = np.exp(1j*120*(np.pi/180))
    a_sq = np.exp(1j*240*(np.pi/180))
    Iph = np.zeros((1,3), dtype=complex)
    a_matr = np.array([[1,1,1],[1, a_sq, a],[1,a,a_sq]])
    Iph = np.dot(a_matr,Iseq)  
    return Iph

# 1.4 the Convert_Sequence2Phase_Voltages() function
def Convert_Sequence2Phase_Voltages(Vseq_mat):
    m = Vseq_mat.shape[0]
    a = np.exp(1j*120*(np.pi/180))
    a_sq = np.exp(1j*240*(np.pi/180))
    Vph_mat = np.zeros((m,3), dtype=complex)
    a_matr = np.array([[1,1,1],[1, a_sq, a],[1,a,a_sq]])
    for row in range(m):
        Vph_mat[row,:] = np.dot(a_matr,Vseq_mat[row,:])  
    return Vph_mat

# ####################################################
# #  Displaying the results in the terminal window   #
# ####################################################
# 2. the DisplayFaultAnalysisResults() function
from tabulate import tabulate
def DisplayFaultAnalysisResults(Iph,Vph_mat,fault_bus,fault_type,Zf,Vf):
    # Fault type
    if fault_type == 0:
        type = "3-Phase Balanced fault"
        phase = 'a,b and c'
    elif fault_type == 1:
        type = "Single Line-to-Ground fault"
        phase = 'a'
    elif fault_type == 2:
        type = "Line-to-Line fault"
        phase = 'b and c'
    elif fault_type == 3:
        type = "Double Line-to-Ground fault"
        phase = 'b and c'
    else:
        print('Unknown Fault Type')


    # Phase current
    mag_and_ang= []
    Iph = np.round(Iph, 3)
    for i in Iph:
        mag_and_ang.append(np.round(np.absolute(i),3))
        mag_and_ang.append(np.round(np.angle(i)*180/np.pi,3))
                    

    # Phase voltages
    mag_a = []  
    mag_b = []
    mag_c = []
    ang_a = []
    ang_b = []
    ang_c = []

    for row in range(Vph_mat.shape[0]):
        mag_a.append(np.round(np.absolute(Vph_mat[row,0]),3))
        mag_b.append(np.round(np.absolute(Vph_mat[row,1]),3))                  
        mag_c.append(np.round(np.absolute(Vph_mat[row,2]),3)) 
        ang_a.append(np.round(np.angle(Vph_mat[row,0])*180/np.pi,3))
        ang_b.append(np.round(np.angle(Vph_mat[row,1])*180/np.pi,3))
        ang_c.append(np.round(np.angle(Vph_mat[row,2])*180/np.pi,3))
    voltage_data = np.vstack((range(1,Vph_mat.shape[0]+1),mag_a,ang_a,mag_b,ang_b,mag_c,ang_c))
    voltage_data = voltage_data.transpose()


    print('==============================================================')
    print('|                  Fault Analysis Results                    |')
    print('==============================================================')
    print('| ', type, ' at Bus ', fault_bus, 'phase ',phase, '.|' )
    print('| Prefault Voltage: Vf = ',Vf, ' (pu)                    |' )
    print('| Fault Impedance: Zf = ',Zf, '  (pu)                    |' )
    print('==============================================================')
    print('|                      Phase Currents                        |')
    print('|                    ------------------                      |')
    print('|------Phase a------|------Phase b------|------Phase c-------|')
    print("|Mag(pu)", "Ang(deg)|", "Mag(pu)", "Ang(deg)|", "Mag(pu)", "Ang(deg)|")
    print(mag_and_ang[0],"  ",mag_and_ang[1],"  ",mag_and_ang[2],"  ",mag_and_ang[3],"  ",mag_and_ang[4],"  ",mag_and_ang[5])
    print("\n")
    print('==============================================================')
    print('|               Phase Line-to-Ground Voltages                 |')
    print('|              ------------------------------                |')
    print('|------Phase a------|------Phase b------|------Phase c-------|')
    print(tabulate(voltage_data, headers=["|Bus|","|Mag(pu)", "Ang(deg)|", "|Mag(pu)", "Ang(deg)|", "|Mag(pu)", "Ang(deg)|"],stralign="center"))
    print("\n")
    print('==============================================================')  
    return