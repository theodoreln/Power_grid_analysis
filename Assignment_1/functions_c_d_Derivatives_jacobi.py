# -*- coding: utf-8 -*-
"""
Created on Tue Feb  6 10:38:55 2024

@author: julia
"""
import numpy as np

#J_dS_dVm derivatives with respect to voltage magnitude
#J_dS_dTheta derivatives with respect to voltage angles

def generate_Derivatives(Ybus, V):
    J_ds_dVm = np.diag(V/np.absolute(V)).dot(np.diag((Ybus.dot(V)).conj())) 
    + np.diag(V).dot(Ybus.dot(np.diag(V/np.absolute(V))).conj())
                
    J_dS_dTheta = 1j*np.diag(V).dot((np.diag(Ybus.dot(V))-Ybus.dot(np.diag(V))).conj())
    
    return J_ds_dVm, J_dS_dTheta
                                                               
J_dS_dVm, J_dS_dTheta = generate_Derivatives(Ybus,V)

def generate_Jacobian(J_dS_dVm, J_dS_dTheta, pv_index, pq_index):
    pvpq_ind = np.append(pv_index, pq_index)
    J_11 = np.real(J_dS_dTheta[np.ix_(pvpq_ind, pvpq_ind)])
    J_12 = np.real(J_dS_dVm[np.ix_(pvpq_ind, pq_index)])
    J_21 = np.imag(J_dS_dTheta[np.ix_(pq_index, pvpq_ind)])
    J_22 = np.imag(J_dS_dVm[np.ix_(pq_index, pq_index)])
    J = np.block([[J_11,J_12],[J_21,J_22]])
    
    return J

J = generate_Jacobian(J_dS_dVm,J_dS_dTheta,pv_index, pq_index)