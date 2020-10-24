# -*- coding: utf-8 -*-
"""
Created on Fri Mar 20 14:50:27 2020

@author: Umberto Bettinardi
         bettinardi96@gmail.com
"""
import sys
import numpy as np

class BoundaryConditions():
    
    """ Class storing Von Neuman and Dirichlet boundary conditions
        
        Methods:
            self.apply(K, F, mesh):
                Set the generalized force to the "node" numbered node to 0.
        Attributes:
        	self.data:
        		contains BCs infos with the following syntax:
        		node number | BC type | local dof | value

        		BC type: "Dirichlet" --> Dirichlet BC, "Neumann" ---> Neumann BC. 
    
    """
    
    
    def __init__(self):
        
        self.data  = None

# ------------------------------------------------------------------------------

    def apply(self, K, F, mesh):
        
        
        for bc in self.data: # this select a row at a time in data

            dof = bc[0]*mesh.dofsNode + bc[2] # dof_global = node number * number of dofs per node + dof_local

            if bc[1] == "Neumann":
                F[dof] += bc[3]
            elif bc[1] == "Dirichlet":
                F -= K[:, [dof]] * bc[3] # modify RHS to take into account prescribed value
                K[:,dof] = np.zeros(len(K)) # zero-out column
                K[dof,:] = np.zeros(len(K)) # zero-out row
                K[dof,dof] = 1 # set diagonal to 1
                F[dof] = bc[3] # enforce value
            else:
                print("BC not recognized")
                sys.exit(0)
        return K, F

    def set_bc(self, BCtype, label, mesh, dof, value):
        """
        Generates a list containing BCs infos

        Parameters
        ----------
        BCtype : Type of BC --> Neumann / Dirichlet.\n
        label : label used in gmsh file to identify the physical region associated\n
                with the given BC.
        mesh : mesh object.
        dof : list containing local dof to wich the BC has to be applied.
        value : BC value.

        Returns
        -------
        listBCs : list containing BCs info. After summing up all the BCs lists the \n
                resulting one has to be transformed into a numpy array and the to be \n
                passed to assigned to the BoundaryConditions.data

        """
        
        tag = mesh.field_data[label][0]
        dim = mesh.field_data[label][1]

        if dim == 0:
            # array containing indices of elements on the boundary
            on_boundary = np.nonzero(mesh.cell_data_dict["gmsh:physical"]["vertex"] == tag)[0]
            # array containing indices of nodes on the boundary
            nodes = np.unique(mesh.cells_dict["vertex"][on_boundary])
            
        elif dim == 1:
            try:
                on_boundary = np.nonzero(mesh.cell_data_dict["gmsh:physical"]["line"] == tag)[0]
                nodes = np.unique(mesh.cells_dict["line"][on_boundary])
            except KeyError:
                pass
            try:
                on_boundary = np.nonzero(mesh.cell_data_dict["gmsh:physical"]["line3"] == tag)[0]
                nodes = np.unique(mesh.cells_dict["line3"][on_boundary])
            except KeyError:
                pass
        else:
            print("Error in BC_engine.set(): dimension not coded")
            
        listBCs = []
        
        for n in nodes:
            for d in range(len(dof)):
                listBCs.append([n,BCtype,dof[d],value[d]])

        return listBCs
    
# -----------------------------------------------------------------------------

    def findDuplicates(self):
        """
        This function finds duplicates in the list of constrained DOFS.
        
        The check is done on the node number and the direction.

        """

        import collections
 
        # extract first and third columns 
        list = [0,2]
        constrainedNodes = [[l[i] for i in list] for l in self.data]
        
        # from list of lists into list of tuples
        nested_lst_of_tuples = [tuple(l) for l in constrainedNodes]
        
        counted=collections.Counter(nested_lst_of_tuples)

        # select only those DOFS that appear more than once in the list
        # https://stackoverflow.com/questions/18807079/selecting-elements-of-a-python-dictionary-greater-than-a-certain-value/18807120
        list = dict((k, v) for k, v in counted.items() if v >1)
    
        if not list:
            pass # list is empty -> each DOF is constrained only once
        else:
            print("ERROR: Some DOFs have been constrained more than once\n", list)
            print("\nThe dictionary contains tuples with the node number and the direction")
            print("and the number of times the DOF appears in the list of constrained DOFs")
            sys.exit()