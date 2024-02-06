"""
The ReadNetworkData library contains a function that reads in network data
stored in a external *.txt file. Used for assignemnt work in the course
46705 - Power Grid Analysis
"""

import numpy as np
import re
import csv

def read_network_data_from_file(filename):
    bus_data = []
    load_data = []
    gen_data = []
    line_data = []
    tran_data = []
    mva_base = 0.0
    data_input_type = '' #label keeping track of the type of model data being read
    #Open the network file and process the data
    with open(filename,newline='')as f:
        csv_reader = csv.reader(f)  # read in all of the data
        for line in csv_reader:
            #check if the line is a header line:
            if "//BEGIN " in line[0]: #New section of data starts with //BEGIN
                if 'BEGIN BUS DATA' in line[0]:
                    data_input_type = 'BUS'
                elif 'BEGIN LOAD DATA' in line[0]:
                    data_input_type = 'LOAD'
                elif 'BEGIN GENERATOR DATA' in line[0]:
                    data_input_type = 'GEN'
                elif 'BEGIN LINE DATA' in line[0]:
                    data_input_type = 'LINE'
                elif 'BEGIN TRANSFORMER DATA' in line[0]:
                    data_input_type = 'TRAN'
                elif 'BEGIN MVA SYSTEM' in line[0]:
                    data_input_type = 'MVA'                   
                else:
                    data_input_type = 'UNSPECIFIED'
                
                continue #read next line which contains data

            
            #Check if there is a comment at the end of the line // is in the data (that is comment) and remove that part
            dummy = []
            for k in line:
                if '//' in k: 
                    dummy.append(k.split('//')[0])
                    break
                else:
                    dummy.append(k)
            line = dummy

            if line[0] is '': # If the line was commented out, skip it..             
                continue
            if line[0].isspace(): # If the line has only whitespace, skip it..             
                continue
            
            
            if data_input_type is 'MVA':
                mva_base = parse_mva_data(line)
            elif data_input_type is 'BUS':
                a = parse_bus_data(line)
                bus_data.append(a)
            elif data_input_type is 'LOAD':
                b = parse_load_data(line)
                load_data.append(b)                
            elif data_input_type is 'GEN':
                d = parse_gen_data(line)
                gen_data.append(d)                
                continue
            elif data_input_type is 'LINE':
               c =  parse_transmission_line_data(line)
               line_data.append(c)
            elif data_input_type is 'TRAN':
               f =  parse_transformer_data(line)
               tran_data.append(f)
            else:
                print('DATA TYPE is unspecified!!!!')
                print('Input not treaded:', line)
            
    

    #create mappping objects from bus_nr to matrix indices and vice verse
    bus_to_ind = {} # empty dictionary
    ind_to_bus = {}
    for ind, line in zip(range(len(bus_data)),bus_data): 
        bus_nr = line[0]
        bus_to_ind[bus_nr] = ind
        ind_to_bus[ind] = bus_nr
            
    return bus_data,load_data,gen_data,line_data,tran_data,mva_base, bus_to_ind, ind_to_bus



#test data for branch:

'''parse_line_data(row_)
parses a string that contains the transmission line data
The input contains following FROM_BUS, TO_BUS, ID, R, X, B_HALF
'''
def parse_transmission_line_data(row_):
    # unpack the values:    
    fr,to,br_id,R,X,B_half = row_[0:6]
    fr = int(fr)    #convert the srting to int
    to = int(to)    #convert the string to int
    br_id = re.findall("'([^']*)'",br_id)[0]
#    br_id = re.findall(r"'.*'",br_id)[0] #regular expression to find the id str
#    br_id = br_id[1:-1]     #get the id out of the string
    R = float(R) #convert the R str to float
    X = float(X) #convert the X str to float
    B_half = float(B_half) #convert the B_half str to float
    return [fr,to,br_id,R,X,B_half] 


def parse_mva_data(row_):
    mva_size = row_[0]
    mva_size = float(mva_size)        #convert string to float
    return mva_size 


def parse_gen_data(row_):
    # unpack the values:    
    bus_nr, mva_size, p_gen = row_[0:3]
    bus_nr = int(bus_nr)     #convert the srting to int
    mva_size = float(mva_size)        #convert string to float
    p_gen = float(p_gen)        #convert string to float
    return [bus_nr, mva_size, p_gen] 



#//BEGIN BUS DATA,(BUS_NR, LABEL, KV_BASE, BUSCODE)
def parse_bus_data(row_):
    # unpack the values:    
    bus_nr, label, kv_base, buscode = row_[0:4]
    bus_nr = int(bus_nr)    #convert the srting to int
    label = re.findall(r'\b.*\b',label)[0] #regular expression to get the label
    kv_base = float(kv_base)        #convert string to float
    buscode = int(buscode)          #convert the buscode str to int
    return [bus_nr, label, kv_base, buscode] 



#//BEGIN LOAD DATA (BUS_NR, P_load MW, Q_load MVAR)  
def parse_load_data(row_):
    # unpack the values:    
    bus_nr, p_ld, q_ld = row_[0:3]
    bus_nr = int(bus_nr)     #convert the srting to int
    p_ld = float(p_ld)        #convert string to float
    q_ld = float(q_ld)        #convert string to float
    return [bus_nr, p_ld, q_ld] 


def parse_transformer_data(row_):
    # unpack the values:    
    fr,to,br_id,R,X,n,ang1 = row_[0:7]
    fr = int(fr)    #convert the srting to int
    to = int(to)    #convert the string to int
    br_id = re.findall("'([^']*)'",br_id)[0]
#    br_id = re.findall(r"'.*'",br_id)[0] #regular expression to find the id str
#    br_id = br_id[1:-1]     #get the id out of the string
    R = float(R) #convert the R str to float
    X = float(X) #convert the X str to float
    n = float(n) #convert the n str to float
    ang1 = float(ang1) #convert the ang1 str to float
    return [fr,to,br_id,R,X,n,ang1] 

    

# This code is executed only if the file is executed as a stand-alone. The below will not
# be executed if the file is imported in other files. 
if __name__=="__main__":    
    file_name = "TestSystem.txt"
    bus_data,load_data,gen_data,line_data, tran_data,mva_base, bus_to_ind, ind_to_bus = \
    read_network_data_from_file(file_name)

    MVA_base = mva_base   #OBS...... should be a part of the network data....
    N = len(bus_data)
    Ybus = np.zeros((N,N),dtype=complex)
    
    for line in line_data:
        bus_fr, bus_to, id_, R,X,B_2 = line #unpack
        ind_fr = bus_to_ind[bus_fr]    
        ind_to = bus_to_ind[bus_to] 
        Z_se = R + 1j*X; Y_se = 1/Z_se
        Y_sh_2 = 1j*B_2
        
        #Update the matrix:
        Ybus[ind_fr,ind_fr]+= Y_se + Y_sh_2
        Ybus[ind_to,ind_to]+= Y_se + Y_sh_2
        Ybus[ind_fr,ind_to]+= -Y_se
        Ybus[ind_to,ind_fr]+= -Y_se
        
        
    #Get transformer data as well...
    bus_kv = []
    buscode = []
    bus_labels = []
    for line in bus_data:
        b_nr, label, kv, code = line
        buscode.append(code)
        bus_labels.append(label)
        bus_kv.append(kv)
        
    buscode = np.array(buscode)
    bus_kv = np.array(bus_kv)

    # Make the bus injection array and the inital guess for the voltages
    Sbus = np.zeros(N,dtype=complex)
    V0 = np.ones(N,dtype=complex)
    
    for line in load_data:
        bus_nr, PLD, QLD = line
        ind_nr = bus_to_ind[bus_nr]
        SLD =(PLD+1j*QLD)/MVA_base
        Sbus[ind_nr] += -SLD # load is a negative injection...
        
    for line in gen_data:
        bus_nr, MVA_size, p_gen = line
        ind_nr = bus_to_ind[bus_nr]
        SLD =(p_gen)/MVA_base
        Sbus[ind_nr] += SLD # load is a negative injection...
        
     
    pq_index = np.where(buscode== 1)[0]# Find indices for all PQ-busses 
    pv_index = np.where(buscode== 2)[0]# Find indices for all PV-busses
    ref = np.where(buscode== 3)[0] # Find index for ref bus
    

    #bus-branch matrices
    N_branches = len(line_data) + len(tran_data)
    br_f = -np.ones(N_branches,dtype=int)
    br_t = -np.ones(N_branches,dtype=int)
    
    Y_fr = np.zeros((N_branches,N),dtype=complex)
    Y_to = np.zeros((N_branches,N),dtype=complex)
        
    for line,i in zip(line_data,range(len(line_data))):
        bus_fr, bus_to, id_, R,X,B_2 = line #unpack
        ind_fr = bus_to_ind[bus_fr]    
        ind_to = bus_to_ind[bus_to] 
        Z_se = R + 1j*X; Y_se = 1/Z_se
        Y_sh_2 = 1j*B_2
        # update the entries
        Y_fr[i,ind_fr] =  Y_se + Y_sh_2       
        Y_fr[i,ind_to] = -Y_se
        Y_to[i,ind_to] =  Y_se + Y_sh_2       
        Y_to[i,ind_fr] = -Y_se
        br_f[i] = ind_fr
        br_t[i] = ind_to
    
    for line,i in zip(tran_data,range(len(line_data),N_branches)):
        bus_fr, bus_to, id_, R,X,n,ang1 = line #unpack
        ind_fr = bus_to_ind[bus_fr]    
        ind_to = bus_to_ind[bus_to] 
        Zeq = R+1j*X; Yeq = 1/Zeq
        c = n*np.exp(1j*ang1/180*np.pi)
        ### adminttance matrix
        Yps_mat = np.zeros((2,2),dtype=complex)
        Yps_mat[0,0] = Yeq/np.abs(c)**2
        Yps_mat[0,1] = -Yeq/c.conj()
        Yps_mat[1,0] = -Yeq/c
        Yps_mat[1,1] = Yeq
        # indices
        ind_ = np.array([ind_fr,ind_to])
        #update
        Ybus[np.ix_(ind_,ind_)] += Yps_mat
        br_f[i] = ind_fr
        br_t[i] = ind_to
        # update the entries
        Y_fr[i,ind_fr] =  Yps_mat[0,0]      
        Y_fr[i,ind_to] =  Yps_mat[0,1]
        Y_to[i,ind_to] =  Yps_mat[1,1]       
        Y_to[i,ind_fr] =  Yps_mat[1,0]
        
    


    
    
