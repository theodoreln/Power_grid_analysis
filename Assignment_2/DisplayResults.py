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

    S_inj = np.round(V*(Ybus.dot(V)).conj(),3) # Generation and load power
    # Sm = np.round(np.absolute(S_inj),2)

    P_inj_g = np.array([]) #Generation active Power
    P_inj_l = np.array([]) #Load active Power
    Q_inj_g = np.array([]) #Generation reactive Power
    Q_inj_l = np.array([]) #Load reactive Power
    loading_inj = np.array([]) #Loading of generator at bus
    S_gen = S_inj+SLD
    S_gen_m = np.round(np.absolute(S_gen),2)
    print(S_inj, S_gen, SLD)

    for index,s in enumerate(S_gen):
        if np.real(s)>0: 
            P_inj_g = np.concatenate((P_inj_g, [str(np.round(np.real(s),3))]))
            Q_inj_g = np.concatenate((Q_inj_g, [str(np.round(np.imag(s),3))]))
            loading_inj = np.concatenate((loading_inj,["{0}%".format(str(np.round(S_gen_m[index]*100/Gen_MVA[index],2)))]))
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

    loading_from = np.array([])
    loading_to = np.array([])
    S_to = V[br_t]*Y_to.dot(V).conj() # Apparent Power flowing into the receiving end busses
    S_to_m = np.absolute(S_to)
    for index,s in enumerate(S_to) :
        loading_to = np.concatenate((loading_to,["{0}%".format(str(np.round(S_to_m[index]*100/br_MVA[index],2)))])) # Loading of line at to branch
    S_from = V[br_f]*Y_from.dot(V).conj() # Apparent Power flowing from the receiving end busses
    S_from_m = np.absolute(S_from)
    for index,s in enumerate(S_from) :
        loading_from = np.concatenate((loading_from,["{0}%".format(str(np.round(S_from_m[index]*100/br_MVA[index],2)))])) # Loading of line at from branch
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
    print("Branch  From    To         From Bus Inject.         To Bus Inject.  ")
    print(tabulate(branch_data, headers=[" # ", "Bus", "Bus", "P(pu)", "Q(pu)", "loading", "P(pu)", "Q(pu)", "loading"],stralign="center"))

