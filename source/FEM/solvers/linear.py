# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 19:59:59 2020

@author: Bettinardi Umberto
@mail:   umberto.bettinardi@studenti.unipd.it
         bettinardi96@gmail.com
"""
import numpy as np
from scipy import linalg
import time

import FEM.FEM_engine as FEM_engine

class SolverError(Exception):
    def __init__(self, text):
        self.text = text
        Exception.__init__(self, text)


def linear(mesh, BCs, MaterialSets, Procedures):
    """
    Linear solver.

    """
            
    # Initializing arrays
    systemDofs = mesh.dofsNode*len(mesh.points)
    K = np.zeros(shape = (systemDofs,systemDofs), dtype = np.float32)
    F = np.zeros(shape = (systemDofs,1), dtype = np.float32)
    
    print("Assemblying global stiffness matrix...")
    
    start_a = time.time()
    
    for e in range(len(mesh.elements)):
    
        k = FEM_engine.stiffness_matrix(e, mesh, MaterialSets)
                
        # Get global dof associate with element e.
        dof = FEM_engine.DofMap(e,mesh)
        
        # Assemble the e-th local matrix into the global one.
        K = FEM_engine.assemble(K,k,dof)
    
        
    end_a = time.time()
    
    print("Global stiffness matrix assembled in {}s".format(end_a - start_a))
    
    Kr = K.copy()
    
    (Kr, F) = BCs.apply(Kr, F, mesh)
    print("Solving F = Ku...")
    
    start_s = time.time()
    
    factorization = Procedures["solver"]["factorization"]
    
    if factorization == "LU":
        LU  = linalg.lu_factor(Kr)
        del Kr
        U   = linalg.lu_solve(LU,F)
        del LU
        end_s = time.time()
        print("\nLU solver: {}s".format(end_s - start_s))
    elif factorization == "cholesky":
        CHO  = linalg.cho_factor(Kr)
        del Kr
        U   = linalg.cho_solve(CHO,F)
        del CHO
        end_s = time.time()
        print("\nCholesky solver: {}s".format(end_s - start_s))
    else:
        LU  = linalg.lu_factor(Kr)
        del Kr
        U   = linalg.lu_solve(LU,F)
        del LU
        end_s = time.time()
        print("\nsolver: {}s".format(end_s - start_s))      
    

    R = np.matmul(K,U)
    
    end_s = time.time()
    
    print("Linear system solved in {}s".format(end_s - start_s))
    
    return U,R,K
        
    