# -*- coding: utf-8 -*-
"""
Created on Sun Mar  3 17:10:02 2024

@author: julia
"""

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
