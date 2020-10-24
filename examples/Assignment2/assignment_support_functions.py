# -*- coding: utf-8 -*-
"""
Created on Sun May 24 15:34:13 2020

@author: Bettinardi Umberto
@mail:   umberto.bettinardi@studenti.unipd.it
         bettinardi96@gmail.com
"""
import sys
sys.path.append("../../source")

from FEM.FEM_utilities import shape_functions, Integration

import numpy as np
import matplotlib.pyplot as plt
import datetime

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

# ---------------------------------------------------------------------------

"""
Infinite plate with circular hole uniformily loaded: Kirsch solution.
Below the functions giving the displacements and the stresses
for the solution are coded.

Parameters:
    R:   circular hole radius.
    E:   Young's modulus.
    v:   Poisson's ratio.
    sig: nominal stress applied at the boundaries in the x direction.
    x,y: coordinates with respect to the hole center.
"""

def u_x(x, y, R, E, v, sig):
    """
    Returns the displacement for the given point along 
    the x-axis
    """
    
    E_ = E/(1-v**2)
    v_ = v/(1-v)
    if x != 0:       
        t = np.arctan2(y,x) #(y/x)
    else:
        t = np.pi/2
    r = np.sqrt(x**2 + y**2)
    
    U_x = ( (1 + v_)/E_ )* sig * ( 
                (1/(1+v_))*r*np.cos(t) +
                (2/(1+v_))*(R**2 / r)*np.cos(t) +
                0.5*(R**2 / r)*np.cos(3*t) - 
                0.5 * (R**4 / r**3)*np.cos(3*t))
    if x == 0:
        U_x = 0
    
    return U_x

def u_y(x, y, R, E, v, sig):
    """
    Returns the displacement for the given point along 
    the y-axis
    """
    
    E_ = E/(1-v**2)
    v_ = v/(1-v)    
    if x != 0:       
        t = np.arctan2(y,x) #(y/x)
    else:
        t = np.pi/2
    r = np.sqrt(x**2 + y**2)

    U_y = ( (1 + v_)/E_ )* sig * ( 
                (-v_ /(1+v_))*r*np.sin(t) -
                ((1 - v_)/(1+v_))*(R**2 / r)*np.sin(t) +
                0.5*(R**2 / r)*np.sin(3*t) - 
                0.5 * (R**4 / r**3)*np.sin(3*t))
    
    if y == 0:
        U_y = 0
    
    return U_y

def s_x(x, y, R, sig): 
    """
    Returns the stress for the given point along 
    the x-axis
    """
    if x != 0:       
        t = np.arctan2(y,x) #(y/x)
    else:
        t = np.pi/2
    r = np.sqrt(x**2 + y**2)
    
    S_x = sig*(1 - (R**2/r**2)*(1.5*np.cos(2*t) + np.cos(4*t)) + 
               1.5*(R**4/r**4)*np.cos(4*t))
    
    return S_x

def s_y(x, y, R, sig):
    """
    Returns the displacement for the given point along 
    the y-axis
    """
    
    if x != 0:       
        t = np.arctan2(y,x) #(y/x)
    else:
        t = np.pi/2
    r = np.sqrt(x**2 + y**2)
    
    S_y = sig*( -(R**2/r**2)*(0.5*np.cos(2*t) - np.cos(4*t)) - 
               1.5*(R**4/r**4)*np.cos(4*t))
    
    return S_y

def s_xy(x, y, R, sig):
    """
    Returns the shear stress for the given point in the x-y plane
    """
    
    
    if x != 0:       
        t = np.arctan2(y,x) #(y/x)
    else:
        t = np.pi/2
    r = np.sqrt(x**2 + y**2)
    
    S_xy = sig*( -(R**2/r**2)*(0.5*np.sin(2*t) + np.sin(4*t)) + 
               1.5*(R**4/r**4)*np.sin(4*t))
    
    return S_xy

# ---------------------------------------------------------------------------

def change_mesh_size(filename, mesh_size):
    """
    This function is specifically designed for the a conventinal
    .geo file structure adopted in this assignment.
    """ 
    with open(filename +'.geo', 'r') as file:
        # read a list of lines into data
        data = file.readlines()
        
    # changing mesh size
    data[1] = 'h = {};\n'.format(mesh_size)
    
    # and write everything back
    with open(filename +'.geo', 'w') as file:
        file.writelines( data )
        
def plot_results(mesh_file, figure_folder, dpi, ELEMEN, ERRORE, DOFS, H):
    """
    Plot the convergence diagramm and store them in a given folder.
    All the figure formatting is predefined

    Parameters
    ----------
    mesh_file :         mesh_file name.
    figure_folder :     figure folder's name.
    dpi:                figure resolution.
    ELEMEN :            array containing number of elements.
    ERRORE :            array storing L2 error norm values.
    DOFS :              array containing number of dof.
    H :                 array containing elements sizes.
    """
    
    if len(ELEMEN)>=2:
        
        # Connvergence rate as defined in the theory of FEM convergence analysis
        conv_rate = float(-np.log(ERRORE[-1]/ERRORE[-2])/np.log(DOFS[-1]/DOFS[-2]))
        
        # Connvergence rate as defined in the assignments 2 text
        conv_rate_H = float(np.log(ERRORE[-1]/ERRORE[-2])/np.log(H[-1]/H[-2]))
    
        
        # Plotting convergence diagramm...
        
        #... with respect to the number of dofs
        fig, ax = plt.subplots()
        ax.grid(which = 'both')
        ax.set_title("Convergence analysis: T3 elements.")
        ax.set_ylabel(r"$ ||\sigma_{e} - \sigma_{h}||_{H^0}$")
        ax.set_xlabel("Number of dofs")
        ax.loglog(DOFS,ERRORE,'-*b')   
        textstr = "Convergence rate: {:4.3f}".format(conv_rate)
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        ax.text(0.5, 0.9, textstr, transform=ax.transAxes, fontsize=14,
                verticalalignment='top', bbox=props)
        
        #... with respect to the number of elements
        fig2, ax2 = plt.subplots()
        ax2.grid(which = 'both')
        ax2.set_title("Convergence analysis: T3 elements.")
        ax2.set_ylabel(r"$ ||\sigma_{e} - \sigma_{h}||_{H^0}$")
        ax2.set_xlabel("Element size")
        ax2.loglog(H,ERRORE,'-*b')
        textstr = "Convergence rate: {:4.3f}".format(conv_rate_H)
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        ax2.text(0.05, 0.9, textstr, transform=ax2.transAxes, fontsize=14,
                verticalalignment='top', bbox=props)
        
        
        # Saving figures in different formats
        fig.savefig(figure_folder+"convergence_{}.eps".format(mesh_file), dpi = 600)
        fig.savefig(figure_folder+"convergence_{}.png".format(mesh_file), dpi = 600)
        fig.savefig(figure_folder+"convergence_{}.pdf".format(mesh_file))
        
        fig2.savefig(figure_folder+"convergence_{}_H.eps".format(mesh_file), dpi = 600)
        fig2.savefig(figure_folder+"convergence_{}_H.png".format(mesh_file), dpi = 600)
        fig2.savefig(figure_folder+"convergence_{}_H.pdf".format(mesh_file))
                
        plt.show()
        
def print_results(mesh_file, ELEMEN, ERRORE, DOFS, H, REL_ERR):    
    """
    Print the results in a .txt file whose name will be in the form:
    mesh_file_result_%Y-%m-%d_%H_%M
    
    Parameters
    ----------
    mesh_file :         mesh_file name.
    ELEMEN :            array containing number of elements.
    ERRORE :            array storing L2 error norm values.
    DOFS :              array containing number of dof.
    H :                 array containing elements sizes.
    REL_ERR :           array storing relative error values.
    """    
    data_file = datetime.datetime.now().strftime("%Y-%m-%d_%H_%M")
    
    if len(ELEMEN)>=2:
        
        # Connvergence rate as defined in the theory of FEM convergence analysis-        
        conv_rate_H = float(np.log(ERRORE[-1]/ERRORE[-2])/np.log(H[-1]/H[-2]))
        
        # Writing results to file
        with open('./results/'+mesh_file+'_results_'+data_file+'.txt', 'w') as result_file:
        
            result_file.write("{}.geo - results\n\n".format(mesh_file))
            result_file.write("@author: Umberto Bettinardi\n"+
                            "         bettinardi96@gmail.com\n\n")
            result_file.write("Analysis date: {}\n\n".format(datetime.datetime.now().strftime("%Y-%m-%d %H:%M")))
            result_file.write("{:<8s}\t{:<8s}\t{:<8s}\t{:<8s}\t {:<9s}\n".format("h-size","Elements", "Dofs", "err_H0", "rel_err"))
            for i in range(len(ELEMEN)):
                result_file.write("{:^8.3f}\t{:^8d}\t{:^8d}\t{:^8.6f}\t{:^9.6f}\n".format(H[i],ELEMEN[i],DOFS[i],ERRORE[i], REL_ERR[i]))
            result_file.write("\nConvergence rate:{}".format(conv_rate_H))    
            
def post_proc(mesh, MaterialSets, U_print, R, E, v, sig):
    
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
    err_H0 = 0
    err_S  = 0
    
    for e in range(ElementNumber):
        
        # Loading element type and material set number
        elMatSet  = mesh.elementMaterialTag[e]
        key = MaterialSets[str(elMatSet)]
        
        # Loading thickness t, Young's modulus E and Poisson ration v.
        E = key['elastic properties'  ]['Young modulus']
        v = key['elastic properties'  ]['poisson ratio']
        t = key['geometric properties'  ]['thickness']
        
        # Loading numerical integration parameters
        evaluation          = key['stiffness matrix']["evaluation"]
        
        if evaluation != 'numerical integration':
            sys.exit()
            
        domain                  = key['stiffness matrix']["domain"]
        rule                    = key['stiffness matrix']["rule"]
        gauss_points_number_str = key['stiffness matrix']["points"] - 1
        
        if gauss_points_number_str < 1:
            gauss_points_number_str = 1
            
        
        # Loading quadrature scheme weights and evaluation_points
        integration = Integration()
        (w, int_p_array) = integration.quadrature_rule(rule, domain, gauss_points_number_str)
        
        # Initializing variables
        elementPoints = mesh.points[mesh.elements[e]][:,:mesh.d]  
        elementNodes  = len(elementPoints[:,0])  
        
        sigma = 0
        sigma_true = np.zeros([1,3])
        
        for i in range(gauss_points_number_str):
            int_p = int_p_array[:,i]
        
            # Computing shape functions: for T3 elements stresses are constant.
            detJ, dN, N = shape_functions(elementPoints, int_p, elementNodes, mesh.d, key['element'])
            
            B = np.array([[dN[0,0],       0, dN[0,1],       0, dN[0,2],       0], 
                          [      0, dN[1,0],       0, dN[1,1],       0, dN[1,2]],
                          [dN[1,0], dN[0,0], dN[1,1], dN[0,1], dN[1,2], dN[0,2]]])
            
            # Stiffness tensor evaluation
            if key['plane deformation'] == 'plane stress':
                
                D = (E/(1 - v**2))*np.array([[1, v,         0],
                                             [v, 1,         0],
                                             [0, 0, 0.5*(1-v)]])
            
            elif key['plane deformation'] == 'plane strain':
                
                D = (E/((1 - 2*v)*(1+v)))*np.array([[(1-v),     v,       0],
                                                    [    v, (1-v),       0],
                                                    [    0,     0, (1-2*v)/2.]])
                
            U_elem = U_print[mesh.elements[e]].reshape([6,1])
            
            sigma_i = np.linalg.multi_dot([D, B, U_elem]).T 
            sigma += sigma_i
        
            # Computing gauss point position
            x_c = N[0,:] @ elementPoints[:,0]#mesh.points[e,0].mean()
            y_c = N[1,:] @ elementPoints[:,1]#mesh.points[e,1].mean()
            
            # Computing "real" stresses
            sigma_true_i = np.zeros([1,3])
            sigma_true_i[0,0] = s_x( x_c, y_c, R, sig)
            sigma_true_i[0,1] = s_y( x_c, y_c, R, sig)
            sigma_true_i[0,2] = s_xy(x_c, y_c, R, sig)
            sigma_true += sigma_true_i
            
            # Computing error
                
            # Integrand function                
            v  = (sigma_true_i[0,:] - sigma_i[0,:]).reshape(3,1) #L2 contribution
            vs = (sigma_true_i[0,:]).reshape(3,1)              # stress RMS
            
            # Adding error contribution        
            err_H0 += 0.5*w[i]* float(( v.T @  v))*detJ*t
            err_S  += 0.5*w[i]* float((vs.T @ vs))*detJ*t
            
        sigma_true /= gauss_points_number_str
        sigma /= gauss_points_number_str
        
        stress[e] = sigma.reshape(1,3)        
        stress_true[e] = sigma_true
        sigma = sigma.T
        von_mises[e] = np.sqrt(sigma[0,0]**2 + sigma[1,0]**2 - sigma[0,0]*sigma[1,0] + 3*sigma[2,0]**2)
        
    err_H0 = np.sqrt(err_H0)
    err_S = np.sqrt(err_S)
    return stress, stress_true, von_mises, err_H0, err_S