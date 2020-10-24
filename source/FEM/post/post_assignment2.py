# -*- coding: utf-8 -*-
"""
Created on Sat Jun 27 10:25:13 2020

@author: Bettinardi Umberto
@mail:   umberto.bettinardi@studenti.unipd.it
         bettinardi96@gmail.com
"""
import numpy as np

def post_proc(mesh, MaterialSets, U):
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
        
        #dofs = int(mesh.dofsNode*len(mesh.elements[e]))
        
        # Loading element type and material set number
        elMatSet  = mesh.elementMaterialTag[e]
        key = MaterialSets[str(elMatSet)]
        
        # Loading thikness t, Young's modulus E and Poisson ration v.
        E = key['elastic properties'  ]['Young modulus']
        v = key['elastic properties'  ]['poisson ratio']
        t = key['geometric properties'  ]['thikness']
        
        
        n1 = mesh.elements[e][0]
        n2 = mesh.elements[e][1]
        n3 = mesh.elements[e][2]
        
        x1 = mesh.points[n1,0]
        y1 = mesh.points[n1,1]
        x2 = mesh.points[n2,0]
        y2 = mesh.points[n2,1]
        x3 = mesh.points[n3,0]
        y3 = mesh.points[n3,1]
        
        y23 = y2 - y3
        y31 = y3 - y1
        y12 = y1 - y2
        x32 = x3 - x2        
        x13 = x1 - x3
        x21 = x2 - x1
        
        A = 0.5*(x13*y23 - y31*x32)
        detJ = 2*A          
        elem_area[e] = A
        
        B = (1/detJ)*np.array([[y23,   0, y31,   0, y12,   0],
                               [  0, x32,   0, x13,   0, x21],
                               [x32, y23, x13, y31, x21, y12]])
        
        # Stiffness tensor evaluation
        if key['plane deformation'] == 'plane stress':
            
            D = (E/(1 - v**2))*np.array([[1, v,         0],
                                         [v, 1,         0],
                                         [0, 0, 0.5*(1-v)]])
        
        elif key['plane deformation'] == 'plane strain':
            
            D = (E/((1 - 2*v)*(1+v)))*np.array([[(1-v),     v,       0],
                                                [    v, (1-v),       0],
                                                [    0,     0, (1-2*v)/2.]])
            
        U_ = U_print[mesh.elements[e]].reshape([6,1])
        #U_ = U_true[mesh.elements[e]].reshape([6,1])
        
        sigma = np.linalg.multi_dot([D, B, U_])
        
        stress[e] = sigma.reshape(1,3)
        
        von_mises[e] = np.sqrt(sigma[0,0]**2 + sigma[1,0]**2 - sigma[0,0]*sigma[1,0] + 3*sigma[2,0]**2)
        
        # Computing element centroid
        x_c = (x1 + x2 + x3)/3
        y_c = (y1 + y2 + y3)/3
        
        # Computing "real" stresses
        stress_true[e,0] = s_x(x_c, y_c, radius, sigma_applied)
        stress_true[e,1] = s_y(x_c, y_c, radius, sigma_applied)
        stress_true[e,2] = s_xy(x_c, y_c, radius, sigma_applied)
        
        # Computin error
        n = 1
        def xx(xi,eta,X):
                
                return X[0] + (X[1]-X[0])*xi + (X[2]-X[0])*eta       
        
        def yy(xi,eta,Y):
            
            return Y[0] + (Y[1]-Y[0])*xi + (Y[2]-Y[0])*eta
        
        quad_scheme_GL = {
                            1: np.array([
                                        [1],
                                        [0.3333333333333333],
                                        [0.3333333333333333]]),
                            3: np.array([
                                         [0.3333333333333333,  
                                          0.3333333333333333,  
                                          0.3333333333333333],
                                         [0, 0.5, 0.5],
                                         [0.5, 0, 0.5]])}
        
        I = 0   # Initialize integral value
        Is = 0
        
        X = mesh.points[mesh.elements[e],0]
        Y = mesh.points[mesh.elements[e],1]
            
        # Computation loop
        
        for i in range(n):
            
            xi_i = quad_scheme_GL[n][1][i]
            eta_i = quad_scheme_GL[n][2][i]
            w_i = quad_scheme_GL[n][0][i]
            
            x_c = xx(xi_i, eta_i, X)
            y_c = yy(xi_i, eta_i, Y)
            
            stress_true[e,0] = s_x(x_c, y_c, radius, sigma_applied)
            stress_true[e,1] = s_y(x_c, y_c, radius, sigma_applied)
            stress_true[e,2] = s_xy(x_c, y_c, radius, sigma_applied)
            
            v = (stress_true[e,:] - stress[e,:]).reshape(3,1)
            vs = (stress_true[e,:]).reshape(3,1)
            
            I +=  0.5*w_i* float((v.T @ v))
            Is +=  0.5*w_i* float((vs.T @ vs))
        
        I *= detJ*t
        Is *= detJ*t
        err_H0 += I
        err_S  += Is