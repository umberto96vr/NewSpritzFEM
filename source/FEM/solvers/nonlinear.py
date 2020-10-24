#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 24 11:09:56 2020

@author: Bettinardi Umberto
@mail:   umberto.bettinardi@studenti.unipd.it
         bettinardi96@gmail.com
"""

import numpy as np
import time
from FEM import FEM_engine
from scipy import linalg

class SolverError(Exception):
    def __init__(self, text):
        self.text = text
        Exception.__init__(self, text)

# -----------------------------------------------------------------------------

# def nonlinear(mesh, BCs, MaterialSets, parameters):
    
#     raise SolverError("Nonlinear solver not yet implemented!")
    
    
def nonlinear(mesh, BCs, MaterialSets, Procedures):
    """
    Linear solver.

    """
    S = 9.6296*1e-13
    L = 10000
    
    systemDofs = mesh.dofsNode*len(mesh.points)
    M = np.zeros(shape = (systemDofs,systemDofs))
    K = np.zeros(shape = (systemDofs,systemDofs))
    F = np.ones((systemDofs,1))*S*L/(mesh.Elements)
    F[0] = S*L/(mesh.Elements)*0.5
    F[-1] = S*L/(mesh.Elements)*0.5
    
    print("Assemblying global stiffness matrix...\n")
    
    start_a = time.time()
    
    for e in range(mesh.Elements):
    
        k = FEM_engine.stiffness_matrix(e, mesh, MaterialSets)
        m = FEM_engine.mass_matrix(e, mesh, MaterialSets)
                
        # Get global dof associate with element e.
        dof = FEM_engine.DofMap(e,mesh)
        
        # Assemble the e-th local matrix into the global one.
        #print("dof : {}\n".format(dof))
        K = FEM_engine.assemble(K,k,dof)
        M = FEM_engine.assemble(M,m,dof)
        #print(M)
        
    
    end_a = time.time()
    
    dt          = Procedures["solver"]['dt']
    time_tot    = Procedures["solver"]['Tot time']
    U_old       = BCs.init_cond.reshape(systemDofs,1)
    
    U = U_old
    A       = M/dt
    K_tilde = A + K    
        
    n_step = int(time_tot/dt)
    
    for t in range(n_step): #CONTINUARE DA QUI
        
        b = F + np.matmul(A,U_old) 
        Kr,b = BCs.apply(K_tilde, b, mesh)
        LU = linalg.lu_factor(Kr)
        U_new   = linalg.lu_solve(LU,b)        
        U = np.append(U, U_new, axis = 1)    
        U_old = U_new   
        
    R = np.matmul(K_tilde, U_new)
    
    return U, R, K