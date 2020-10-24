# -*- coding: utf-8 -*-
"""
Created on Wed May 20 20:19:44 2020

@author: Umberto Bettinardi
         bettinardi96@gmail.com
"""
import numpy as np

#----------------------------------------------------------

class PostProcessorError(Exception):
    def __init__(self, text):
        self.text = text
        Exception.__init__(self, text)

#----------------------------------------------------------

def run_post(mesh, MaterialSets, U_vector):
    
    # Retrieving element number
    ElementNumber = len(mesh.elements) 
        
    # Retrieving strain size
    if mesh.d == 1:
        strainComponents = 1
    elif mesh.d == 2:
        strainComponents = 3
    elif mesh.d == 3:
        strainComponents = 6
    # else:
    #     raise FemE("Error, spatial dimension different from 1, 2, 3: mesh.d = {}\n".format(mesh.d))
        
    D = np.zeros([strainComponents, strainComponents])
    
    stress      = np.zeros([ElementNumber, strainComponents])
    von_mises   = np.zeros(len(stress))
    elem_area   = np.zeros(len(stress))
    
    for e in range(len(mesh.elements)):
        
        #dofs = int(mesh.dofsNode*len(mesh.elements[e]))
        
        # Loading element type and material set number
        elMatSet  = mesh.elementMaterialTag[e]
        key = MaterialSets[str(elMatSet)]
        
        # Loading thikness t, Young's modulus E and Poisson ration v.
        E = key['elastic properties'  ]['Young modulus']
        v = key['elastic properties'  ]['poisson ratio']
        
        
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
                                                [    0,     0, 0.5*(1-2*v)]])
            
        U = U_vector[mesh.elements[e]].reshape([6,1])
        
        sigma = D @ B @ U
        
        stress[e,:] = sigma.reshape(1,3)
       
        von_mises[e] = np.sqrt(sigma[0,0]**2 + sigma[1,0]**2 - sigma[0,0]*sigma[1,0] + 3*sigma[2,0]**2)
        
    return stress, von_mises

def test():
    print("test")