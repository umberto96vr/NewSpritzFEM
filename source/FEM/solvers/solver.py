# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 19:51:19 2020

@author: Umberto Bettinardi
         bettinardi96@gmail.com
"""
class SolverError(Exception):
    def __init__(self, text):
        self.text = text
        Exception.__init__(self, text)


def run(mesh, BCs, MaterialSets, Procedures):
    """
    Solve the algebric equation K U = f

    Parameters
    ----------
    mesh : object containing mesh informations.\n
    BCs : object containing boundary conditions informations.\n
    material : object containing material informations.\n
    parameters : dictionary containing solver settings.\n

    Returns
    -------
    U : Displacements vector.\n
    R : Reaction forces vector.\n
    K : Stiffness matrix.\n

    """
        
    # Defining possible solver options
    possible_solvers = ['linear','linear_disp', 'nonlinear']

    solverType = Procedures["solver"]["type"]
    
    if solverType in possible_solvers:
    
        exec("from FEM.solvers."+solverType+" import "+solverType)
                
        (U, R, K) = eval(solverType+"(mesh, BCs, MaterialSets, Procedures)")
        
        return U, R, K
    
    else:
        raise SolverError("Invalid solver type!")

