# -*- coding: utf-8 -*-
"""
Created on Sat May 16 14:52:12 2020

@author: Umberto Bettinardi
         bettinardi96@gmail.com
"""
import sys
import time
# import warnings
# warnings.filterwarnings("ignore")
sys.path.append("../../source")

import numpy as np
np.set_printoptions(threshold=6)
import meshio
from FEM import mesh_engine 
from FEM import BC_engine
import FEM.solvers.solver as solver
from assignment_support_functions import u_x, u_y
from assignment_support_functions import change_mesh_size, print_results, plot_results, post_proc
import matplotlib.pyplot as plt
import datetime


data_file = datetime.datetime.now().strftime("%Y-%m-%d_%H_%M")
figure_folder = "./figures/"
mesh_file = "assignment2"

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

#---------------------------------------------------------------

# Problem data
sigma_applied = 100
Young = 2.1*1e7
Poisson = 0.29
radius = 2.5
thickness = 1

#---------------------------------------------------------------
# MODEL PARAMETERS
#---------------------------------------------------------------

print("\n--> Pre-processing...\n")
print("Initializing model parameters...")

Procedures = {"solver": {"type": "linear",
                         "factorization": "cholesky"}, 
              
              "postprocessing":{"undeformed mesh": {"output": "screen"},
                                "deformed mesh": {"magnification factor": 50, "output": "screen"},
                                "axial stress": {"output": "screen"}
                                }
              }
#---------------------------------------------------------------

# Material info
print("Defining material sets...")

MaterialSets = {
    'piastra': {'element': 'triangle',
          'plane deformation'    : 'plane strain',
          'material behavior': 'isotropic linear elastic',
          'elastic properties'   : {"Young modulus" : Young, "poisson ratio" : Poisson},
          'geometric properties' : {'thickness': thickness},
          'stiffness matrix'     : {'evaluation': 'numerical integration',
                                            'domain'    : 'triangular',
                                            'rule'      : 'Gauss-Legendre',
                                            'points'    : 1}
          
          }
    }
#---------------------------------------------------------------

H = [1.7, 0.71, 0.35, 0.17]#, 0.17, 0.1]#, 0.13]#, 0.12, 0.11, 0.1]
ERRORE  = []
ELEMEN  = []
DOFS    = []
REL_ERR = []

for h in H:

    print("\nMesh size: {}\n".format(h))
        
    # Mesh info
    print("Meshing...")
    
    # Upgrading mesh size in .geo file
    change_mesh_size(mesh_file, h)   # see assignment_support_function.change_mesh_size()
    mesh = mesh_engine.GMSH(mesh_file)  
    #--------------------------------------------------------------------------------------------------
    
    # BCs enforcement    
    print("Applying BCs...\n")
    BCs = BC_engine.BoundaryConditions()
    
    # Non-homogeneous Dirichlet BCs enforcement on the boundary
    name = "EnforcedDisplacement"
    tag = mesh.field_data[name][0]    
    on_boundary = np.nonzero(mesh.cell_data_dict["gmsh:physical"]["line"] == tag)[0]
    nodes_enf_disp = np.unique(mesh.cells_dict["line"][on_boundary])
    
    BCdata = []
 
    for i in nodes_enf_disp:
        
        x = mesh.points[i,0]
        y = mesh.points[i,1]
        
        ux = u_x(x, y, radius, Young, Poisson, sigma_applied)
        uy = u_y(x, y, radius, Young, Poisson, sigma_applied)        
                
        BCdata.append([int(i), "Dirichlet", 0, ux ])
        BCdata.append([int(i), "Dirichlet", 1, uy ])
    del x,y, ux, uy
    
    # In case symmetries are exploited corresponding BCs are enforced
    if mesh_file == 'assignment2':
        BCdata += BCs.set_bc("Dirichlet", 'x-symmetry', mesh, [1], [0.])
        BCdata += BCs.set_bc("Dirichlet", 'y-symmetry', mesh, [0], [0.])
        
    BCs.data = BCdata
    
    #---------------------------------------------------------------
    # SOLVER
    
    print("--> Solving...\n")

    start = time.time()    
    U, R, K = solver.run(mesh, BCs, MaterialSets, Procedures)    
    end = time.time()    
    
    print("Total time: {}s".format(end-start))    
    print("\n-- Post-processing...\n")
        
    U_print = U.copy().reshape(len(mesh.points),mesh.d)
    print("   . solution U:\n{}\n".format(U_print))
    
    # -----------------------------------------------------------------------------
    # Stress computing 
    
    # Retrieving true displacements
    U_true = U_print.copy()
    for i in range(len(mesh.points)):
        
        # Points coordinates
        x_point = mesh.points[i,0]
        y_point = mesh.points[i,1]
        
        # Points displacements
        uu_x = u_x(x_point, y_point, radius, Young, Poisson, sigma_applied)
        uu_y = u_y(x_point, y_point, radius, Young, Poisson, sigma_applied)
        U_true[i,0] = uu_x
        U_true[i,1] = uu_y
        
    del x_point, y_point, uu_x, uu_y
    
    # Retrieving stress and error norms
    stress, stress_true, von_mises, err_H0, err_S = post_proc(mesh, MaterialSets, U_print,
                                                              radius, Young, Poisson, sigma_applied)
    
    # Storing values in arrays    
    REL_ERR.append(err_H0/err_S)
    ERRORE.append(err_H0)
    ELEMEN.append(len(mesh.elements))
    DOFS.append(len(K))
     
    # Output to paraview
    U_print = np.append(U_print, np.zeros((len(mesh.points),1)),axis = 1)
    U_true  = np.append( U_true, np.zeros((len(mesh.points),1)),axis = 1)
    
    mesh.point_data = {'Displacements_FEM': U_print, 
                       'Displacements_true': U_true}  
    
    mesh.cell_data = {"stress_true": stress_true,"stress_fem": stress, "von_mises": von_mises}      
    cells = {'triangle': mesh.elements}
    
    # Writing the vtk file
    meshio.write_points_cells("./fileVtk/{}/".format(
        MaterialSets['piastra']['plane deformation'].replace(" ","_"))+ mesh_file +"_e{}.vtk".format(
            len(mesh.elements)), mesh.points, cells, mesh.point_data, mesh.cell_data, binary=True)
         
    
# Plotting results
plot_results(mesh_file, figure_folder, 600, ELEMEN, ERRORE, DOFS, H)

# Printing results to file
print_results(mesh_file, ELEMEN, ERRORE, DOFS, H, REL_ERR)

# ...see assignment_support_function for further details

#del I, Is, err_H0, err_S, x_c, y_c, v, vs