# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 18:30:32 2020

@author: Bettinardi Umberto
@mail:   umberto.bettinardi@studenti.unipd.it
         bettinardi96@gmail.com
"""

# This reimplements gmsh/demos/boolean/boolean.geo in Python.
import os
import sys
import time
sys.path.append("../../source")

import numpy as np

from FEM import mesh_engine 
from FEM import BC_engine
import FEM.solvers.solver as solver
import FEM.post.postprocess as post
import meshio

np.set_printoptions(precision=4)

# Parameters
print("\n--> Pre-processing...\n")
print("Initializing model parameters...")

Procedures = {"solver": {"type": "linear"}, 
              
              "postprocessing":{"undeformed mesh": {"output": "screen"},
                                "deformed mesh": {"magnification factor": 50, "output": "screen"},
                                "axial stress": {"output": "screen"}
                                }
              }

#---------------------------------------------------------------

# Material info
print("Defining material sets...")
MaterialSets = {

    '112': {'element': 'triangle',
          'plane deformation'    : 'plane stress',
          'material behavior'    : 'isotropic linear elastic',
          'elastic properties'   : {"Young modulus" : 3*1e7, "poisson ratio" : 0.25},
          'geometric properties' : {'thikness': 0.5},
          'stiffness matrix'     : {'evaluation': 'closed form'}
          
          },
    '2': {'element': 'triangle',
          'plane deformation'    : 'plane stress',
          'material behavior'    : 'isotropic linear elastic',
          'elastic properties'   : {"Young modulus" : 2.06*1e6, "poisson ratio" : 0.3},
          'geometric properties' : {'thikness': 0.5},
          'stiffness matrix'     : {'evaluation': 'closed form'}
          
          }
    }

#---------------------------------------------------------------
print("Meshing...")
mesh_file = "bimaterial2"
os.system("gmsh -2 " + mesh_file + ".geo" + " -o " + mesh_file + ".msh")
meshIO = meshio.read(mesh_file + ".msh")

mesh = mesh_engine.Mesh()
mesh.elements = meshIO.cells_dict["triangle"]
mesh.elementMaterialTag = meshIO.cell_data_dict['gmsh:physical']['triangle']
mesh.elementType = np.array([["triangle"] for i in range(len(mesh.elements))])

mesh.points = meshIO.points

mesh.NodesElement = 3
mesh.dofsNode     = 2
mesh.d            = 2
mesh.Nodes        = len(mesh.points)


# Load info
print("Applying BCs...\n")

BCs = BC_engine.BoundaryConditions()

P = 1000

Dirichlet = 1
Neuman = 0
BCs.data = np.array([[0, Dirichlet, 0,  0],
                     [0, Dirichlet, 1,  0],
                     [1, Dirichlet, 1,  0],
                     [3, Dirichlet, 0,  0],
                     [3, Dirichlet, 1,  0],
                     [2, Neuman   , 1, -P]])#,
                     #[4, Neuman   , 1, q1],
                     #[0, Neuman   , 1, q1]])

#---------------------------------------------------------------

print("--> Solving...\n")

# Global solver time assesment
start = time.time()

U, R, K = solver.run(mesh, BCs, MaterialSets, Procedures)

end = time.time()


print("Total time: {}s".format(end-start))

print("\n-- Post-processing...\n")

U_print = U.copy()

print("   . solution U:\n{}\n".format(U_print.reshape(mesh.Nodes,2)))
print("   . reaction forces R:\n{}\n".format(R.reshape(mesh.Nodes,2)))

meshio.write(mesh_file+".vtk", mesh)