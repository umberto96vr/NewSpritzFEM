# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 19:51:19 2020

@author: Umberto Bettinardi
         bettinardi96@gmail.com
"""
import sys

import numpy as np

import FEM.FEM_utilities as FEM_utilities
#import FEM.material_data as material_data

# import time

# Additional utilities
integration = FEM_utilities.Integration()

#----------------------------------------------------------

class FemEngineError(Exception):
    def __init__(self, text):
        self.text = text
        Exception.__init__(self, text)

#----------------------------------------------------------

def stiffness_matrix(e, mesh, MaterialSets):
    """
    Return the local stiffness matrix.
    
    Parameters
    ----------
    e:           element number
    mesh:        object containing mesh info
    material:    object containing material info
    utilities:   object containing integration quadrature scheme info
                 and shape functions derivatives
    
    Returns
    ---------
    k : local stiffness matrix
    
    """
    # Loading element type and material set number
    elType    = mesh.elementType[e]
    elMatSet  = mesh.elementMaterialTag[e]
    
    key = MaterialSets[str(elMatSet)]   
    
    # Loading stiffness matrix evaluation infos
    if key['stiffness matrix']["evaluation"] == 'numerical integration':
        evaluation          = "numerical integration"
        domain              = key['stiffness matrix']["domain"]
        rule                = key['stiffness matrix']["rule"]
        gauss_points_number = key['stiffness matrix']["points"]
        
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
    if mesh.d == 1:
        strainComponents = 1
    elif mesh.d == 2:
        strainComponents = 3
    elif mesh.d == 3:
        strainComponents = 6
    else:
        raise FemEngineError("Error, spatial dimension different from 1, 2, 3: mesh.d = {}\n".format(mesh.d))
        
    D = np.zeros([strainComponents, strainComponents])
    
    
    if evaluation == 'closed form':  
                
        if elType == 'spring':
            
            k_s = key['elastic properties'  ]['Spring stiffness']
    
            return k_s*np.array([(1, -1), (-1, 1)])
        
        elif elType == 'bar':
            
            # Loading cross section A and Young's modulus E
            A = key["geometric properties"]["area"]
            E = key["elastic properties"  ]["Young modulus"]
            
            # Loading points coordinates
            n1 = mesh.elements[e][0]
            n2 = mesh.elements[e][1]
                        
            # Computing element length                    
            L = mesh.points[n2,0] - mesh.points[n1,0]
            
            if L < 0:
                
                raise FemEngineError("Error: bar lenght is negative!")
                                
            else:
                
                return ((E*A)/L)*np.array([(1, -1), (-1, 1)]) 
        
        elif elType == 'euler_bernoulli':
            
            # Loading cross section A and Young's modulus E
            A = key["geometric properties"]["area"]
            J = key["geometric properties"]["inertia"]
            E = key["elastic properties"  ]["Young modulus"]
            
            # Loading points coordinates
            n1 = mesh.elements[e][0]
            n2 = mesh.elements[e][1]                  
            x1 = mesh.points[n1,0]
            x2 = mesh.points[n2,0]            
            y1 = mesh.points[n1,1]
            y2 = mesh.points[n2,1]
                        
            # Computing element length
            L = np.sqrt( pow(x2-x1,2)+ pow(y2-y1,2) )
            
            # Computing cosine
            
            #TODO 
            # Check that ly and my are evaluated correctly
            lx = (x2-x1)/L
            mx = (y2-y1)/L
            ly = -mx
            my = lx #np.sqrt(1 - mx**2)
            
            N = (E*A)/L
            M = (12*E*J)/(L**3)
            T = (6*E*J)/(L**2) 
            S =  (2*E*J)/(L)
            
            K_0 = np.array([[ N, 0,    0, -N,  0,   0],
                            [ 0, M,    T,  0, -M,   T],
                            [ 0, T,  2*S,  0, -T,   S],
                            [-N, 0,    0,  N,  0,   0],
                            [ 0,-M,   -T,  0,  M,  -T],
                            [ 0, T,    S,  0, -T, 2*S]]) 
            
            T_mat = np.array([[lx, mx, 0,  0,  0, 0],
                              [ly, my, 0,  0,  0, 0],
                              [ 0,  0, 1,  0,  0, 0],
                              [ 0,  0, 0, lx, mx, 0],
                              [ 0,  0, 0, ly, my, 0],
                              [ 0,  0, 0,  0,  0, 1]])
            
            return T_mat.T @ K_0 @ T_mat
        
        elif elType == 'euler_bernoulli3D':
            
            raise FemEngineError("Error: element currently work in progress!")
            
            # Currently WORK IN PROGRESS
            
            # Loading cross section A and Young's modulus E
            A = key["geometric properties"]["area"]
            Jp = key["geometric properties"]["polar inertia"]
            I1 = key["geometric properties"]["Inertia 1"]
            I2 = key["geometric properties"]["Inertia 2"]
            E = key["elastic properties"  ]["Young modulus"]
            G = key["elastic properties"  ]["Shear modulus"]
            
            # Loading points coordinates
            n1 = mesh.elements[e][0]
            n2 = mesh.elements[e][1]                
            x1 = mesh.points[n1,0]
            x2 = mesh.points[n2,0]   
            y1 = mesh.points[n1,1]
            y2 = mesh.points[n2,1]
            z1 = mesh.points[n1,2]
            z2 = mesh.points[n2,2]
            
                        
            # Computing element length
            L = np.sqrt( pow(x2-x1,2) + pow(y2-y1,2) + pow(z2-z1,2) )
            
            # Computing cosine
            
            #TODO 
            # Cx = (x2-x1)/L
            # Cy = (y2-y1)/L
            # Cz = (z2-z1)/L
            
            
            l1 = np.sqrt(pow(x2-x1,2) + pow(z2-z1,2))/L
            m1 = np.sqrt(1 - pow(l1,2))
            l2 = (x2 - x1) / np.sqrt(pow(x2-x1,2) + pow(z2-z1,2))
            m2 = np.sqrt(1 - pow(l2,2))
            
            N   = (E*A)/L
            Nt  = (G*Jp)/L
            M1  = (12*E*I1)/(L**3)
            T1  = (6*E*I1)/(L**2) 
            S1  =  (2*E*I1)/(L)
            M2  = (12*E*I2)/(L**3)
            T2  = (6*E*I2)/(L**2) 
            S2  =  (2*E*I2)/(L)
            
            # CORRECTED VERSION
            K_0 = np.array([[  N,   0,   0,   0,    0,    0, -N,   0,   0,   0,    0,    0],
                            [  0,  M1,   0,   0,    0,   T1,  0, -M1,   0,   0,    0,   T1],
                            [  0,   0,  M2,   0,  -T2,    0,  0,   0, -M2,   0,  -T2,    0],
                            [  0,   0,   0,  Nt,    0,    0,  0,   0,   0, -Nt,    0,    0],
                            [  0,   0, -T2,   0, 2*S2,    0,  0,   0,  T2,   0,   S2,    0],
                            [  0,  T1,   0,   0,    0, 2*S1,  0, -T1,   0,   0,    0,   S1],
                            [ -N,   0,   0,   0,    0,    0,  N,   0,   0,   0,    0,    0],
                            [  0, -M1,   0,   0,    0,  -T1,  0,  M1,   0,   0,    0,  -T1],
                            [  0,   0, -M2,   0,   T2,    0,  0,   0,  M2,   0,   T2,    0],
                            [  0,   0,   0, -Nt,    0,    0,  0,   0,   0,  Nt,    0,    0],
                            [  0,   0, -T2,   0,   S2,    0,  0,   0,  T2,   0, 2*S2,    0],
                            [  0,  T1,   0,   0,    0,   S1,  0, -T1,   0,   0,    0, 2*S1]])
            
            
           
            # NEW VERSION
            T1 = np.array([[],
                           [],
                           [],])
            
            T_mat = np.array([[  l1*l2, m1,  l1*m2,      0,  0,      0,      0,  0,      0,      0,  0,      0],
                              [ -l2*m1, l1, -m1*m2,      0,  0,      0,      0,  0,      0,      0,  0,      0],
                              [    -m2,  0,     l2,      0,  0,      0,      0,  0,      0,      0,  0,      0],
                              [      0,  0,      0,  l1*l2, m1,  l1*m2,      0,  0,      0,      0,  0,      0],
                              [      0,  0,      0, -l2*m1, l1, -m1*m2,      0,  0,      0,      0,  0,      0],
                              [      0,  0,      0,    -m2,  0,     l2,      0,  0,      0,      0,  0,      0],
                              [      0,  0,      0,      0,  0,      0,  l1*l2, m1,  l1*m2,      0,  0,      0],
                              [      0,  0,      0,      0,  0,      0, -l2*m1, l1, -m1*m2,      0,  0,      0],
                              [      0,  0,      0,      0,  0,      0,    -m2,  0,     l2,      0,  0,      0],
                              [      0,  0,      0,      0,  0,      0,      0,  0,      0,  l1*l2, m1,  l1*m2],
                              [      0,  0,      0,      0,  0,      0,      0,  0,      0, -l2*m1, l1, -m1*m2],
                              [      0,  0,      0,      0,  0,      0,      0,  0,      0,    -m2,  0,     l2]])
            
            
            
            return T_mat.T @ K_0 @ T_mat
        
        elif elType == 'bar2d':
            
            # Loading cross section A and Young's modulus E
            A = key["geometric properties"]["area"]
            E = key["elastic properties"  ]["Young modulus"]
                        
            # Loading points coordinates
            n1 = mesh.elements[e][0]
            n2 = mesh.elements[e][1]                  
            x1 = mesh.points[n1,0]
            x2 = mesh.points[n2,0]            
            y1 = mesh.points[n1,1]
            y2 = mesh.points[n2,1]
            
            # Computing element length
            L = np.sqrt( pow(x2-x1,2)+ pow(y2-y1,2) )
            
            # Computing cosine
            l = (x2-x1)/L
            m = (y2-y1)/L
                        
            return ((E*A)/L)*np.array([[ l*l,  m*l, -l*l, -m*l],
                                       [ m*l,  m*m, -m*l, -m*m],
                                       [-l*l, -m*l,  l*l,  m*l],
                                       [-m*l, -m*m,  m*l,  m*m]]) 
        elif elType == 'diff_bar':

            # Loading 'stiffnes matrix' properties
            diff =  key["elastic properties"  ]["diffusivity"] 

            # Loading points coordinates
            n1 = mesh.elements[e][0]
            n2 = mesh.elements[e][1]
                        
            # Computing element length                    
            L = mesh.points[n2,0] - mesh.points[n1,0]
            
            if L < 0:
                
                raise FemEngineError("Error: bar lenght is negative!")
                                
            else:
                
                return (diff/L)*np.array([(1, -1), (-1, 1)]) 
                                  
        elif elType == 'bar_pullout':
            
            # Element associated with the pullout problem of a bar immersed in a
            # medium of stiffness k
            
            # Loading cross section A and Young's modulus E
            A  = key['geometric properties']['area']
            E  = key['elastic properties'  ]['Young modulus']
            kk = key['elastic properties'  ]['k']
                        
            # Loading points coordinates
            n1 = mesh.elements[e,0]
            n2 = mesh.elements[e,1]
                        
            # Computing element length                    
            L = mesh.points[n2,0] - mesh.points[n1,0]
                        
            if L < 0:
                
                raise FemEngineError("Error: bar lenght is negative!")
                
                
            else:
                if nodesInElement == 2:          
                    
                    k_p = (kk*A*L)/6 * np.array([[2, 1],
                                                 [1, 2]])
                        
                    k = ((E*A)/L)*np.array([[1, -1], [-1, 1]]) + k_p#k_pullout(L,L,A,kk) - k_pullout(0,L,A,kk)
                    
                    return k
                elif nodesInElement == 3:
                    
                    k_p = (kk*A*L)/30 * np.array([[ 4,  2, -1],
                                                  [ 2, 16,  2],
                                                  [-1,  2,  4]])
                    k_s = ((E*A)/(3*L)) * np.array([[7,  -8,  1],
                                                    [-8, 16, -8],
                                                    [1,  -8,  7]])
                    return k_s + k_p
                
        elif elType == 'quad':
            
            # Loading cross section A and Young's modulus E
            t = key["geometric properties"]["thickness"]
            E = key["elastic properties"  ]["Young modulus"] 
            v = key["elastic properties"  ]["poisson ratio"] 
            
            # Loading points coordinates
            n1 = mesh.elements[e][0]
            n2 = mesh.elements[e][1]
            n3 = mesh.elements[e][2]
            n4 = mesh.elements[e][3]
            
            x1 = mesh.points[n1,0]
            x2 = mesh.points[n2,0]
            x3 = mesh.points[n3,0]
            x4 = mesh.points[n4,0]
            y1 = mesh.points[n1,1]
            y2 = mesh.points[n2,1]
            y3 = mesh.points[n3,1]
            y4 = mesh.points[n4,1]
            
            # Computing element dimension            
            a = x2 - x1
            b = y3 - y2
            gam = a/b
            p1 = (1+v)*gam
            p2 = (1-3*v)*gam
            p3 = 2 + (1-v)*(gam**2)
            p4 = 2*(gam**2) + (1-v)
            p5 = (1-v)*(gam**2) - 4
            p6 = (1-v)*(gam**2) - 1
            p7 = 4*(gam**2) - (1-v)
            p8 = gam**2 + (1-v)
            
            k_lu = np.array([[2*p3, 3*p1, 2*p5, -3*p2, -2*p3, -3*p1, -4*p6,  3*p2],
                             [0,    2*p4, 3*p2,  4*p8, -3*p1, -2*p4, -3*p2, -2*p7],
                             [0,       0, 2*p3, -3*p1, -4*p6, -3*p2, -2*p3,  3*p1],
                             [0,       0,    0,  2*p4,  3*p2, -2*p7,  3*p1, -2*p4],
                             [0,       0,    0,     0,  2*p3,  3*p1,  2*p5, -3*p2],
                             [0,       0,    0,     0,     0,  2*p4,  3*p2,  4*p8],
                             [0,       0,    0,     0,     0,     0,  2*p3, -3*p1],
                             [0,       0,    0,     0,     0,     0,     0,  2*p4]])
            k = k_lu+k_lu.T
            
            k *= (E*t) / (24 * gam *(1 - v**2) ) 
            
            return k
        
        elif elType == 'triangle':
            
            # Material set loading
            key = MaterialSets[str(elMatSet)]
            
            # Loading thikness t, Young's modulus E and Poisson ration v.
            t = key['geometric properties']['thickness']
            E = key['elastic properties'  ]['Young modulus']
            v = key['elastic properties'  ]['poisson ratio']
            
            # Stiffness tensor evaluation
            if key['plane deformation'] == 'plane stress':
                
                D = (E/(1 - v**2))*np.array([[1, v,         0],
                                             [v, 1,         0],
                                             [0, 0, 0.5*(1-v)]])
            
            elif key['plane deformation'] == 'plane strain':
                
                D = (E/((1 - 2*v)*(1+v)))*np.array([[(1-v),     v,           0],
                                                    [    v, (1-v),           0],
                                                    [    0,     0, 0.5*(1-2*v)]])
                
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
            x32 = x3 - x2
            y31 = y3 - y1
            x13 = x1 - x3
            y12 = y1 - y2
            x21 = x2 - x1
            
            A = 0.5*(x13*y23 - y31*x32)
            
            detJ = 2*A 
            
            if detJ < 0:
                print("ERRORE!!!!!")
            
            B = (1/detJ)*np.array([[y23,   0, y31,   0, y12,   0],
                                   [  0, x32,   0, x13,   0, x21],
                                   [x32, y23, x13, y31, x21, y12]])
            
            K = t * A * (B.T @ D @ B)
            
            return K
        
        elif elType == "triangle_prof":
            
            # T6 element example coded by Prof. Angelo Simone
            
            E = key['elastic properties']["Young modulus"]
            nu = key['elastic properties']["poisson ratio"]
            thickness = key['geometric properties']["thickness"]

            # plane deformation: plane stress or plane strain?
            if key["plane deformation"] == 'plane stress':
                D = E/(1-nu**2) * np.array([[1,nu,0],[nu,1,0],[0,0,(1-nu)/2.]])
            elif key["plane deformation"] == 'plane strain':
                D = E/(1+nu)/(1-2*nu) * np.array([[1-nu,nu,0],[nu,1-nu,0],[0,0,(1-2*nu)/2.]])
            else:
                print("stiffness_matrix(): closed form: Cannot define D!")
                sys.exit()

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
            x32 = x3 - x2
            y13 = y1 - y3
            x23 = x2 - x3
            
            area = 0.5* ( (x2*y3 - x3*y2) + (x3*y1 - x1*y3) + (x1*y2 - x2*y1) )
            detJ = x13*y23 - y13*x23

            B = 1/detJ * np.array([[y23, 0,   y31, 0,   y12, 0],
                                   [0,   x32, 0,   x13, 0,   x21],
                                   [x32, y23, x13, y31, x21, y12]])

            # IFEM.Ch15.pdf, formula 15.22 
            matrix = thickness * area * np.linalg.multi_dot([B.T, D, B])
           
            return matrix           
        
        
        elif elType == "tethraedal":
            
            # Material set loading
            key = MaterialSets[str(elMatSet)]
            
            # Loading thikness t, Young's modulus E and Poisson ration v.
            E = key['elastic properties'  ]['Young modulus']
            v = key['elastic properties'  ]['poisson ratio']
            
            # Stiffness tensor evaluation
            D = E /( (1+v) * (1-2*v) )*np.array([1-v,   v,   v,     0,     0,     0],
                                                [  v, 1-v,   v,     0,     0,     0],
                                                [  v,   v, 1-v,     0,     0,     0],
                                                [  0,   0,   0, 0.5-v,     0,     0],
                                                [  0,   0,   0,     0, 0.5-v,     0],
                                                [  0,   0,   0,     0,     0, 0.5-v],)
            
            # Retrieving element coordinates
            n1 = mesh.elements[e][0]
            n2 = mesh.elements[e][1]
            n3 = mesh.elements[e][2]
            n4 = mesh.elements[e][3]
            x1 = mesh.points[n1,0]
            y1 = mesh.points[n1,1]
            z1 = mesh.points[n1,2] 
            x2 = mesh.points[n2,0]
            y2 = mesh.points[n2,1]
            z2 = mesh.points[n2,2] 
            x3 = mesh.points[n3,0]
            y3 = mesh.points[n3,1]
            z3 = mesh.points[n3,2] 
            x4 = mesh.points[n3,0]
            y4 = mesh.points[n3,1]
            z4 = mesh.points[n4,2] 
            
            # Retrieving shape function constants
            x21 = x2 - x1
            x31 = x3 - x1
            x41 = x4 - x1
            y21 = y2 - y1
            y31 = y3 - y1
            y41 = y4 - y1
            y32 = y3 - y2
            y42 = y4 - y2
            y43 = y4 - y3
            z21 = z2 - z1
            z31 = z3 - z1
            z41 = z4 - z1
            z32 = z3 - z2
            z42 = z4 - z2
            z43 = z4 - z3            
        
            a1 =  y2*z43 - y3*z42 + y4*z32
            a2 = -y1*z43 + y3*z41 - y4*z31
            a3 =  y1*z42 - y2*z41 + y4*z21
            a4 = -y1*z32 + y2*z31 - y3*z21
            
            b1 = -x2*z43 + x3*z42 - x4*z32
            b2 =  x1*z43 - x3*z41 + x4*z31
            b3 = -x1*z42 + x2*z41 - x4*z21
            b4 =  x1*z32 - x2*z31 + x3*z21
            
            c1 =  x2*y43 - x3*y42 + x4*y32
            c2 = -x1*y43 + x3*y41 - x4*y31
            c3 =  x1*y42 - x2*y41 + x4*y21
            c4 = -x1*y32 + x2*y31 - x3*y21
            
            detJ = x21*(y31*z41 - y41*z31) + y21*(x41*z31 - x31*z41) + z21*(x31*y41 - x41*y31)
            vol  = detJ/6  #Volume
            
            # B matrix
            B = (1/detJ) * np.array([a1,  0,  0, a2,  0,  0, a3,  0,  0, a4,  0,  0],
                                    [ 0, b1,  0,  0, b2,  0,  0, b3,  0,  0, b4,  0],
                                    [ 0,  0, c1,  0,  0, c2,  0,  0, c3,  0,  0, c4],
                                    [b1, a1,  0, b2, a2,  0, b3, a3,  0, b4, a4,  0],
                                    [ 0, c1, b1,  0, c2, b2,  0, c3, b3,  0, c4, b4],
                                    [c1,  0, a1, c2,  0, a2, c3,  0, a3, c4,  0, a4])
            
            K = vol * ( B.T @ D @ B )
            
            return K            
        
        else:
            raise FemEngineError("Error: element type {} not yet avaiable!".format(mesh.elementType[e]))
            
            
    elif evaluation == 'numerical integration':
        
        # These elements are based on their isoparametric formulation, their 
        # stiffness matricies are computed by means of Gauss-Legendre quadrature
        # of  \int_{V}{B^T*D*B*dV}
        
        if elType == 'bar':
            
            # Loading quadrature scheme weights and evaluation_points
            (w, int_p) = integration.quadrature_rule(rule, domain, gauss_points_number)
            
            # Loading elastic properties
            A = key['geometric properties']['area']
            E = key['elastic properties'  ]['Young modulus']
            
            # Defining volume factor
            Vf = A
            
            # Defining the stiffness tensor D
            D[0,0] = E
            
            # Loading points coordinates
            n1 = mesh.elements[e][0]
            n2 = mesh.elements[e][1]
                        
            # Computing element length                    
            L = mesh.points[n2,0] - mesh.points[n1,0]
            
            # Initializing local stiffness matrix
            elementDofs = mesh.NodesElement * mesh.dofsNode
            k = np.zeros([elementDofs, elementDofs])
            
            if L < 0:
                
                raise FemEngineError("Error: bar lenght is negative!")
                
            else:
                
                # Loop for numerical integration
                
                for i in range(len(int_p)):
                    
                    # Computing shape functions derivatives vector                     
                    (detJ, B, N) = FEM_utilities.shape_functions(mesh.elem_coordinates(e), int_p[i], 
                                                  mesh.NodesElement, mesh.d, elType)

                    # i-step numerical contribution to stiffness matrix evaluation
                    # jacobian * dot(B_transposed, B) * Young modulus * Area * weight
                    
                    # dN has to be expressed in terms of the global coordinate system (x,y,z)
                    # The integrand function is F(x,y)*|J|*d_xi * d_eta that is
                    #   B(x,y)^T * D * B(x,y) * |J| * d_xi * d_eta 
                    # = B(xi,eta)^T * (J^-1)^T* D * B(xi,eta) * J^T * |J| * d_xi * d_eta 
                    
                    k += detJ * B.T @ D @ B * Vf * w[i]
                     
            return k
        
        elif elType == 'bar_pullout':
            
            # See the description for the "closed formulation" version
            
            # Loading quadrature scheme weights and evaluation_points
            (w, int_p) = integration.quadrature_rule(rule, domain, gauss_points_number)
            
            # Loading elastic properties
            A  = key['geometric properties']['area']
            E  = key['elastic properties'  ]['Young modulus']
            kk = key['elastic properties'  ]['k']
            
            # Defining volume factor
            Vf = A
            
            # Defining the stiffness tensor D
            D[0,0] = E
            
            # Loading points coordinates
            n1 = mesh.elements[e][0]
            n2 = mesh.elements[e][1]
                        
            # Computing element length                    
            L = mesh.points[n2,0] - mesh.points[n1,0]
            
            # Initializing local stiffness matrix
            elementDofs = mesh.NodesElement * mesh.dofsNode
            k = np.zeros([elementDofs, elementDofs])
            
            if L < 0:
                
                raise FemEngineError("Error: bar lenght is negative!")
                
            else:
                
                # Loop for numerical integration
                
                for i in range(len(int_p)):
                    
                    # Computing shape functions derivatives vector                     
                    (detJ, B, N) = FEM_utilities.shape_functions(elementPoints, int_p[i], 
                                                  nodesInElement, mesh.d, elType)

                    # i-step numerical contribution to stiffness matrix evaluation
                    # jacobian * dot(B_transposed, B) * Young modulus * Area * weight
                    
                    # dN has to be expressed in terms of the global coordinate system (x,y,z)
                    # The integrand function is F(x,y)*|J|*d_xi * d_eta that is
                    #   B(x,y)^T * D * B(x,y) * |J| * d_xi * d_eta 
                    # = B(xi,eta)^T * (J^-1)^T* D * B(xi,eta) * J^T * |J| * d_xi * d_eta 
                    
                    k += detJ * (B.T @ D @ B * Vf  + kk * N.T @ N)* w[i]
            
            return k         
        
        elif elType == 'quad':
            
            # Q4 element
            
            # Loading quadrature scheme weights and evaluation_points
            (w, int_p) = integration.quadrature_rule(rule, domain, gauss_points_number)
            
            # Material set loading
            key = MaterialSets[str(elMatSet)]
            
            # Loading thikness t, Young's modulus E and Poisson ration v.
            t = key['geometric properties']['thickness']
            E = key['elastic properties'  ]['Young modulus']
            v = key['elastic properties'  ]['poisson ratio']
            
            # Stiffness tensor evaluation
            if key['plane deformation'] == 'plane stress':
                
                D = (E/(1 - v**2))*np.array([[1, v,         0],
                                              [v, 1,         0],
                                              [0, 0, 0.5*(1-v)]])
            
            elif key['plane deformation'] == 'plane strain':
                
                D = (E/((1 - 2*v)*(1+v)))*np.array([[(1-v),     v,       0],
                                                    [    v, (1-v),       0],
                                                    [    0,     0, 0.5*(1-2*v)]])
            
            # Defining volume factor
            Vf = t          
                        
            # Initializing local stiffness matrix
            k = np.zeros([elementDofs, elementDofs])
            
            
            # Loop for numerical integration
            
            for i in range(len(int_p)):
                
                for j in range(len(int_p)):
                    
                    # Compute integration points location
                    int_p_vect = [int_p[i], int_p[j]]
                
                    # Computing shape functions derivatives vector                     
                    (detJ, dN, N) = FEM_utilities.shape_functions(elementPoints, int_p_vect, 
                                                                  nodesInElement, mesh.d, elType)
    
                    B = np.array([[dN[0,0],       0, dN[0,1],       0, dN[0,2],       0, dN[0,3],       0], 
                                  [      0, dN[1,0],       0, dN[1,1],       0, dN[1,2],       0, dN[1,3]],
                                  [dN[1,0], dN[0,0], dN[1,1], dN[0,1], dN[1,2], dN[0,2], dN[1,3], dN[0,3]]])
                    
                    k += detJ * (B.T @ D @ B) * Vf * w[i]*w[j]                    
                     
            return k
        elif elType == 'triangle6':
            
            # T6 element
            
            # Loading quadrature scheme weights and evaluation_points
            (w, int_p) = integration.quadrature_rule(rule, domain, gauss_points_number)
            
            # Material set loading
            key = MaterialSets[str(elMatSet)]
            
            # Loading thikness t, Young's modulus E and Poisson ration v.
            t = key['geometric properties']['thickness']
            E = key['elastic properties'  ]['Young modulus']
            v = key['elastic properties'  ]['poisson ratio']
            
            # Stiffness tensor evaluation
            if key['plane deformation'] == 'plane stress':
                
                D = (E/(1 - v**2))*np.array([[1, v,         0],
                                              [v, 1,         0],
                                              [0, 0, 0.5*(1-v)]])
            
            elif key['plane deformation'] == 'plane strain':
                
                D = (E/((1 - 2*v)*(1+v)))*np.array([[(1-v),     v,       0],
                                                    [    v, (1-v),       0],
                                                    [    0,     0, 0.5*(1-2*v)]])
            
            # Defining volume factor
            Vf = t          
                        
            # Initializing local stiffness matrix
            k = np.zeros([elementDofs, elementDofs])
            
            
            # Loop for numerical integration
            
            for i in range(len(int_p)):
                    
                    # Compute integration points location
                    int_p_vect = int_p[i]
                
                    # Computing shape functions derivatives vector                     
                    (detJ, dN, N) = FEM_utilities.shape_functions(elementPoints, int_p_vect, 
                                                                  nodesInElement, mesh.d, elType)
    
                    B = np.array([[dN[0,0],       0, dN[0,1],       0, dN[0,2],       0, dN[0,3],       0,dN[0,4],       0, dN[0,5],       0], 
                                  [      0, dN[1,0],       0, dN[1,1],       0, dN[1,2],       0, dN[1,3],      0, dN[1,4],       0, dN[1,5]],
                                  [dN[1,0], dN[0,0], dN[1,1], dN[0,1], dN[1,2], dN[0,2], dN[1,3], dN[0,3],dN[1,4], dN[0,4], dN[1,5], dN[0,5]]])
                    
                    k += 0.5*detJ * (B.T @ D @ B) * Vf * w[i]
                     
            return k
        elif elType == 'triangle':
            
            # T3 element
            
            # Loading quadrature scheme weights and evaluation_points
            (w, int_p) = integration.quadrature_rule(rule, domain, gauss_points_number)
            
            # Material set loading
            key = MaterialSets[str(elMatSet)]
            
            # Loading thickness t, Young's modulus E and Poisson ration v.
            t = key['geometric properties']['thickness']
            E = key['elastic properties'  ]['Young modulus']
            v = key['elastic properties'  ]['poisson ratio']
            
            # Stiffness tensor evaluation
            if key['plane deformation'] == 'plane stress':
                
                D = (E/(1 - v**2))*np.array([[1, v,         0],
                                              [v, 1,         0],
                                              [0, 0, 0.5*(1-v)]])
            
            elif key['plane deformation'] == 'plane strain':
                
                D = (E/((1 - 2*v)*(1+v)))*np.array([[(1-v),     v,       0],
                                                    [    v, (1-v),       0],
                                                    [    0,     0, 0.5*(1-2*v)]])
            
            # Defining volume factor
            Vf = t          
                        
            # Initializing local stiffness matrix
            k = np.zeros([elementDofs, elementDofs])
            
            
            # Loop for numerical integration
            
            for i in range(gauss_points_number):
                    
                    # Compute integration points location
                    int_p_vect = int_p[:,i]
                
                    # Computing shape functions derivatives vector                     
                    (detJ, dN, N) = FEM_utilities.shape_functions(elementPoints, int_p_vect, 
                                                                  nodesInElement, mesh.d, elType)
    
                    B = np.array([[dN[0,0],       0, dN[0,1],       0, dN[0,2],       0], 
                                  [      0, dN[1,0],       0, dN[1,1],       0, dN[1,2]],
                                  [dN[1,0], dN[0,0], dN[1,1], dN[0,1], dN[1,2], dN[0,2]]])
                    
                    k += 0.5*detJ * (B.T @ D @ B) * Vf * w[i]
                     
            return k
        else:
            raise FemEngineError("stiffness matrix evaluation key '{}' is not defined!".format(evaluation))

#----------------------------------------------------------

def mass_matrix(e, mesh, MaterialSets):
    # WORK IN PROGRESS -> we'll need these for future eigenvalues retrievals for 
    # modal analysis
    
    """
    Return the local mass matrix.
    
    Parameters
    ----------
    e:           element number
    mesh:        object containing mesh info
    material:    object containing material info
    utilities:   object containing integration quadrature scheme info
                 and shape functions derivatives
    
    Returns
    ---------
    k : local stiffness matrix
    
    """
    # Loading element type and material set number
    elType    = mesh.elementType[e]
    elMatSet  = mesh.elementMaterialTag[e]
    
    key = MaterialSets[str(elMatSet)]   
    
    # Loading stiffness matrix evaluation infos
    if key['stiffness matrix']["evaluation"] == 'numerical integration':
        evaluation          = "numerical integration"
        domain              = key['stiffness matrix']["domain"]
        rule                = key['stiffness matrix']["rule"]
        gauss_points_number = key['stiffness matrix']["points"]
        
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
    if mesh.d == 1:
        strainComponents = 1
    elif mesh.d == 2:
        strainComponents = 3
    elif mesh.d == 3:
        strainComponents = 6
    else:
        raise FemEngineError("Error, spatial dimension different from 1, 2, 3: mesh.d = {}\n".format(mesh.d))
        
    D = np.zeros([strainComponents, strainComponents])
    
    if elType == 'diff_bar':

            # Loading 'stiffnes matrix' properties
            ro =  key["mass properties"]["density"] 

            # Loading points coordinates
            n1 = mesh.elements[e][0]
            n2 = mesh.elements[e][1]
                        
            # Computing element length                    
            L = mesh.points[n2,0] - mesh.points[n1,0]
            
            if L < 0:
                
                raise FemEngineError("Error: bar lenght is negative!")
                                
            else:
                
                return ((ro*L)/6)*np.array([(2, 1), (1, 2)]) 


#----------------------------------------------------------

def DofMap(e, mesh):
    """
    Return a vector containing the global dof numbers associated\n
    with the local dofs for the given element\n
    
    Parameters
    ----------
    e : element number\n
    mesh: instance of Mesh() class containing mesh infos
    
    Returns
    ---------
    dof : dof map for the given element
    """   

    dofsNode = mesh.dofsNode
    nodesInElement = len(mesh.elements[e])
       
    # dof vector generation via list comprehension
    dofs = [mesh.elements[e][i]*dofsNode + j for i in range(nodesInElement) for j in range(dofsNode)]
    return dofs


#----------------------------------------------------------

def assemble(K,k,dof):
    """
    Assemble in the global stiffness matrix "K" the local stiffness 
    matrix "k" associated with the given element associated to the dofs 
    map "dof"
    
    Parameters
    ----------
    K  : Global stiffness matrix \n
    k  : local stiffness matrix  \n
    dof: dof map for the given element  \n
        
    Returns
    -------
    K  : Assembled global stiffness matrix
    """
    # Computing the element number of degrees of freedom
    elementDofs = len(dof)
    
    if elementDofs <= 2:
        
        for i in range(elementDofs):
            ii = dof[i]
            for j in range(elementDofs):
                jj = dof[j]
                K[ii,jj] += k[i,j]
    else:
        K[np.ix_(dof,dof)] += k

    return K

# -----------------------------------------------------------------------------

def stress_bar(e, mesh, MaterialSets, U):
    """Compute the axial stress in a 1D bar in a 2D space"""
      
    #TODO: make this part of a general function (as stiffness_matrix)

    elMat = mesh.elementMaterialTag[e]
    key = MaterialSets[str(elMat)]    
    elType = mesh.elementType[e]   
    evaluation = key['stiffness matrix']["evaluation"]
    
    # extract system dofs associated with element
    dofs = DofMap(e,mesh) 

    if evaluation == 'closed form':    
          
        if elType == "spring":
            print("spring: you should not be here")
            sys.exit()
    
        elif elType == "bar":
            print("bar: you should not be here")
            sys.exit()
        
        elif elType == "bar2d":
            Young = key['elastic properties'  ]['Young modulus']
           
            n1 = mesh.elements[e][0]
            n2 = mesh.elements[e][1]
            x1 = mesh.points[n1,0]
            y1 = mesh.points[n1,1]
            x2 = mesh.points[n2,0]
            y2 = mesh.points[n2,1]
            L = np.sqrt( pow(x2-x1,2) + pow(y2-y1,2) )
            
            if L < 0.0:
                print("stiffness_matrix(): closed form: Oh dear, bar",e,"has a negative length!")
                sys.exit()

            # define sine and cosine of angle between local and global axes
            c = (x2-x1)/L
            s = (y2-y1)/L
            
            # compute elongation, strain, and stress            
            elongation = [-c, -s, c, s] @ U[dofs]
            axial_strain = elongation / L
            axial_stress = Young * axial_strain

            return axial_stress
        
        else:
            print("stress_bar(): closed form: I don't know elType =",elType)
            sys.exit()

    else:
        
        print("you should not be here")
        sys.exit()
            
# -----------------------------------------------------------------------------
"""
The followings are the expressions for the displacement and the stress field in 
the x-y plane related to the Kirsh plate under plane strain hypothesis.
"""
def u_x(x, y, R, E, v, sig):
    
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
    
    return U_x

def u_y(x, y, R, E, v, sig):
    
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
    
    return U_y

def s_x(x, y, R, sig): 
    if x != 0:       
        t = np.arctan2(y,x) #(y/x)
    else:
        t = np.pi/2
    r = np.sqrt(x**2 + y**2)
    
    S_x = sig*(1 - (R**2/r**2)*(1.5*np.cos(2*t) + np.cos(4*t)) + 
               1.5*(R**4/r**4)*np.cos(4*t))
    
    return S_x

def s_y(x, y, R, sig):
    
    if x != 0:       
        t = np.arctan2(y,x) #(y/x)
    else:
        t = np.pi/2
    r = np.sqrt(x**2 + y**2)
    
    S_y = sig*( -(R**2/r**2)*(0.5*np.cos(2*t) - np.cos(4*t)) - 
               1.5*(R**4/r**4)*np.cos(4*t))
    
    return S_y

def s_xy(x, y, R, sig):
    
    if x != 0:       
        t = np.arctan2(y,x) #(y/x)
    else:
        t = np.pi/2
    r = np.sqrt(x**2 + y**2)
    
    S_xy = sig*( -(R**2/r**2)*(0.5*np.sin(2*t) + np.sin(4*t)) + 
               1.5*(R**4/r**4)*np.sin(4*t))
    
    return S_xy
