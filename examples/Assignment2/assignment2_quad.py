# -*- coding: utf-8 -*-
"""
Created on Wed May 27 10:26:58 2020

@author: Umberto Bettinardi
         bettinardi96@gmail.com

Solve exercise 4 of assignment 2 with Q4 elements
"""


import sys
import time
sys.path.append("../../source")

import numpy as np
import meshio
from FEM import mesh_engine 
from FEM import BC_engine
import FEM.solvers.solver as solver
from FEM import FEM_utilities
from assignment_support_functions import u_x, u_y, s_x, s_y, s_xy, change_mesh_size  
import matplotlib.pyplot as plt
import datetime

data_file = datetime.datetime.now().strftime("%Y-%m-%d_%H_%M")
figure_folder = "./figures/"
mesh_file = "assignment2_quad"

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

#---------------------------------------------------------------

# Problem data
sigma_applied = 100
Young = 2.1*1e7
Poisson = 0.29
radius = 2.5
thikness = 1

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
    'piastra': {'element': 'quad',
          'plane deformation'    : 'plane strain',
          'material behavior': 'isotropic linear elastic',
          'elastic properties'   : {"Young modulus" : Young, "poisson ratio" : Poisson},
          'geometric properties' : {'thikness': thikness},
          'stiffness matrix'     : {'evaluation': 'numerical integration',
                                            'domain'    : 'rectangular',
                                            'rule'      : 'Gauss-Legendre',
                                            'points'    : 2}
          
          }
    }
#---------------------------------------------------------------
# 100, 500, 2000, 8000 elements
H = [ 1.1, 0.49, 0.24, 0.12]#, 0.1, 0.09]#, 0.085, 0.078, 0.07]
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
        
    # In case symmetries are exploited corresponding BCs are enforced
    if mesh_file != 'assignment2_full':
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
    
    # Retrieving element number
    ElementNumber = len(mesh.elements) 
        
    # Retrieving strain size
    if mesh.d == 1:
        strainComponents = 1
    elif mesh.d == 2:
        strainComponents = 3
    elif mesh.d == 3:
        strainComponents = 6
        
    # Initializin results arrays
    stress      = np.zeros([ElementNumber, strainComponents])
    stress_true = np.zeros([ElementNumber, strainComponents])
    von_mises   = np.zeros(len(stress))
    elem_area   = np.zeros(len(stress))
    err_H0 = 0
    err_S  = 0
    
    for e in range(len(mesh.elements)):
        
        # Loading element type and material set number
        elType    = mesh.elementType[e]
        elMatSet  = mesh.elementMaterialTag[e]
        
        key = MaterialSets[str(elMatSet)]   
        
        # Loading stiffness matrix evaluation infos
        if key['stiffness matrix']["evaluation"] == 'numerical integration':
            evaluation          = "numerical integration"
            domain              = key['stiffness matrix']["domain"]
            rule                = key['stiffness matrix']["rule"]
            gauss_points_number = 3#key['stiffness matrix']["points"]
            
        elif key['stiffness matrix']["evaluation"] == 'closed form':
            evaluation          = "closed form"
            domain              = None
            rule                = None
            gauss_points_number = None
        else:
            raise TypeError("Keyword not recognized for stiffness matrix evaluation:\n", key['stiffness matrix']["evaluation"])
        
        # Initializing variables
        elementPoints  = mesh.points[mesh.elements[e]][:,:mesh.d]     
        nodesInElement = len(elementPoints)
        elementDofs = nodesInElement * mesh.dofsNode
        
        # Retrieving strain size
        strainComponents = 3
        
        D = np.zeros([strainComponents, strainComponents])
        
        # Loading element type and material set number
        elMatSet  = mesh.elementMaterialTag[e]
        key = MaterialSets[str(elMatSet)]
        
        # Loading thikness t, Young's modulus E and Poisson ration v.
        E = key['elastic properties'  ]['Young modulus']
        v = key['elastic properties'  ]['poisson ratio']
        t = key['geometric properties'  ]['thikness']        
                
        # Stiffness tensor evaluation
        if key['plane deformation'] == 'plane stress':
            
            D = (E/(1 - v**2))*np.array([[1, v,         0],
                                          [v, 1,         0],
                                          [0, 0, 0.5*(1-v)]])
        
        elif key['plane deformation'] == 'plane strain':
            
            D = (E/((1 - 2*v)*(1+v)))*np.array([[(1-v),     v,       0],
                                                [    v, (1-v),       0],
                                                [    0,     0, (1-2*v)/2.]])
            
        U_ = U_print[mesh.elements[e]].reshape([8,1])
        
        # Defining volume factor
        #Vf = t          
        
        # Loading quadrature scheme weights and evaluation_points
        integration = FEM_utilities.Integration()
        (w, int_p) = integration.quadrature_rule(rule, domain, gauss_points_number)        
            
        # Initializing local stiffness matrix
        sigma = np.zeros([1,strainComponents])
                
        I = 0
        Is = 0
        for i in range(len(int_p)):
            for j in range(len(int_p)):
                
                # Compute integration points location
                int_p_vect = [int_p[i], int_p[j]]
            
                # Computing shape functions derivatives vector                     
                (detJ, dN, N) = FEM_utilities.shape_functions(elementPoints, int_p_vect, 
                                                              nodesInElement, mesh.d, elType)
                
                x_c = np.dot(N[0,:],elementPoints[:,0])
                y_c = np.dot(N[1,:],elementPoints[:,1])
                 
                stress_true[e,0] = s_x(x_c, y_c, radius, sigma_applied)
                stress_true[e,1] = s_y(x_c, y_c, radius, sigma_applied)
                stress_true[e,2] = s_xy(x_c, y_c, radius, sigma_applied)
                
                B = np.array([[dN[0,0],       0, dN[0,1],       0, dN[0,2],       0, dN[0,3],       0], 
                              [      0, dN[1,0],       0, dN[1,1],       0, dN[1,2],       0, dN[1,3]],
                              [dN[1,0], dN[0,0], dN[1,1], dN[0,1], dN[1,2], dN[0,2], dN[1,3], dN[0,3]]])
                
                # Taking the average value
                sigma += (1/(len(int_p)**2)) * np.linalg.multi_dot([D,B,U_]).reshape(1,3)
                
                v = (stress_true[e,:] - np.linalg.multi_dot([D,B,U_]).reshape(1,3)).reshape(3,1)
                vs = stress_true[e,:].reshape(3,1)
                
                I += detJ*w[i]*w[j]*t* float(v.T @ v)       
                Is += detJ*w[i]*w[j]*t* float(vs.T @ vs) 
            stress[e,:] = sigma.reshape(1,3)
            von_mises[e] = np.sqrt(sigma[0,0]**2 + sigma[0,1]**2 - sigma[0,0]*sigma[0,1] + 3*sigma[0,2]**2)
        I *= t
        Is *= t
        err_H0 += I
        err_S += Is
        
    err_H0 = np.sqrt(err_H0)
    err_S = np.sqrt(err_S)
    REL_ERR.append(err_H0/err_S)
    ERRORE.append(err_H0)
    ELEMEN.append(len(mesh.elements))
    DOFS.append(len(K))
    #--------------------------------------------------------------------------------------------------

    # PARAVIEW .vtk file generation
    
    # Generate additional data for meshio for .vtk file generation
    U_vtk = np.append(U_print, np.zeros([len(U_print[:,0]),1]), axis = 1)
    
    mesh.point_data = {'Displacements_FEM': U_vtk}
    #mesh.cell_data = {"stress_true": stress_true,"stress_fem": stress, "von_mises": von_mises}      
    cells = {'quad': mesh.elements}
    mesh.cell_data = {"stress_fem": stress, "von_mises": von_mises}
    
    # Write the vtk file
    meshio.write_points_cells("./fileVtk/{}/".format(
        MaterialSets['piastra']['plane deformation'].replace(" ","_"))+ mesh_file +"_e{}.vtk".format(
            len(mesh.elements)), mesh.points, cells, mesh.point_data, mesh.cell_data, binary=True)
    
    #-------------------------------------------------------------------------------------------------- 

if len(ELEMEN)>=2:
    #--------------------------------------------------------------------------------------------------        
    conv_rate = -np.log(ERRORE[-1]/ERRORE[-2])/np.log(DOFS[-1]/DOFS[-2])
    #--------------------------------------------------------------------------------------------------    
    
    # Writing results to file
    result_file = open('./results/'+mesh_file+'_results_'+data_file+'.txt', 'w')
    
    result_file.write("{}.geo - results\n\n".format(mesh_file))
    result_file.write("@author: Umberto Bettinardi\n"+
                      "         bettinardi96@gmail.com\n\n")
    result_file.write("Analysis date: {}\n\n".format(datetime.datetime.now().strftime("%Y-%m-%d %H:%M")))
    result_file.write("\t{:.8s}\t  {:.8s}\t    {:.8s}\t\t{:.8s}\t {:.9s}\n".format("h-size","Elements", "Dofs", "err_H0", "rel_err"))
    for i in range(len(ELEMEN)):
        result_file.write("{:8.3f}\t{:8d}\t{:8d}\t{:8.6f}\t{:9.6f}\n".format(H[i],ELEMEN[i],DOFS[i],ERRORE[i], REL_ERR[i]))
    result_file.write("\nConvergence rate:\t{}".format(conv_rate))    
    result_file.close()    
    #-------------------------------------------------------------------------------------------------- 
    
    # Plotting convergence diagramm
    fig, ax = plt.subplots()
    ax.grid(which = 'both')
    ax.set_title("Convergence analysis: Q4 elements.")
    ax.set_ylabel(r"$ ||\sigma_{e} - \sigma_{h}||_{H^0}$")
    ax.set_xlabel("Number of dofs")
    ax.loglog(DOFS,ERRORE,'-*b')   
    textstr = "Convergence rate: {:4.3f}".format(conv_rate)
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    ax.text(0.5, 0.9, textstr, transform=ax.transAxes, fontsize=14,
            verticalalignment='top', bbox=props)
    fig.savefig(figure_folder+"convergence_{}.eps".format(mesh_file), dpi = 600)
    fig.savefig(figure_folder+"convergence_{}.png".format(mesh_file), dpi = 600)
    fig.savefig(figure_folder+"convergence_{}.pdf".format(mesh_file))
    
    conv_rate = np.log(ERRORE[-1]/ERRORE[-2])/np.log(H[-1]/H[-2])
    fig2, ax2 = plt.subplots()
    ax2.grid(which = 'both')
    ax2.set_title("Convergence analysis: Q4 elements.")
    ax2.set_ylabel(r"$ ||\sigma_{e} - \sigma_{h}||_{H^0}$")
    ax2.set_xlabel("Element size")
    ax2.loglog(H,ERRORE,'-*b')
    textstr = "Convergence rate: {:4.3f}".format(conv_rate)
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    ax2.text(0.05, 0.9, textstr, transform=ax2.transAxes, fontsize=14,
            verticalalignment='top', bbox=props)
    fig2.savefig(figure_folder+"convergence_{}_H.eps".format(mesh_file), dpi = 600)
    fig2.savefig(figure_folder+"convergence_{}_H.png".format(mesh_file), dpi = 600)
    fig2.savefig(figure_folder+"convergence_{}_H.pdf".format(mesh_file))
    
    
    plt.show()



    