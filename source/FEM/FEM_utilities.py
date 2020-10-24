# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 17:43:32 2020

@author: Umberto Bettinardi
         bettinardi96@gmail.com
"""
import sys
import numpy as np


class QuadratureRuleError(Exception):
    def __init__(self, text):
        self.text = text
        Exception.__init__(self, text)
        
class InversionError(Exception):
    def __init__(self, text):
        self.text = text
        Exception.__init__(self, text)
        
class FemEngineError(Exception):
    def __init__(self, text):
        self.text = text
        Exception.__init__(self, text)

class Integration():
    
    quad_scheme = {
        
        'Gauss-Legendre': { 'line' : {
                                    1: np.array([
                                                [2.],
                                                [0.]]),
                                    2: np.array([
                                                [-0.5773502691896258,
                                                 0.5773502691896258],
                                                [1.,1.]]),
                                    3: np.array([
                                                [-0.7745966692414834, 
                                                 0., 
                                                 0.7745966692414834],
                                                [0.5555555555555556, 
                                                 0.8888888888888888, 
                                                 0.5555555555555556]]), 
                                    4: np.array([
                                                [-0.8611363115940526, 
                                                 -0.33998104358485626, 
                                                 0.33998104358485626, 
                                                 0.8611363115940526],
                                                [0.34785484513745385,
                                                 0.6521451548625462,
                                                 0.6521451548625462,
                                                 0.34785484513745385]])},
                    'rectangular' : {
                                    1: np.array([
                                                [0.],
                                                [2.]]),
                                    2: np.array([
                                                [-0.5773502691896258,
                                                 0.5773502691896258],
                                                [1.,1.]]),
                                    3: np.array([
                                                [-0.7745966692414834, 
                                                 0., 
                                                 0.7745966692414834],
                                                [0.5555555555555556, 
                                                 0.8888888888888888, 
                                                 0.5555555555555556]]), 
                                    4: np.array([
                                                [-0.8611363115940526, 
                                                 -0.33998104358485626, 
                                                 0.33998104358485626, 
                                                 0.8611363115940526],
                                                [0.34785484513745385,
                                                 0.6521451548625462,
                                                 0.6521451548625462,
                                                 0.34785484513745385]])},
                     'triangular' : {
                                    1: np.array([
                                                [1.],
                                                [0.3333333333333333],
                                                [0.3333333333333333]]),
                                    3: np.array([
                                                [0.3333333333333333,  
                                                 0.3333333333333333,  
                                                 0.3333333333333333],
                                                [0., 0.5, 0.5],
                                                [0.5, 0., 0.5]])}
                             }
        }
        
    
    
    def __init__(self):
        pass
    
    def quadrature_rule(self, rule, domain, n_points):
        
        if rule in self.quad_scheme.keys():
            
            if domain in self.quad_scheme[rule].keys():
                
                if n_points in self.quad_scheme[rule][domain].keys():
                    
                    if domain == "triangular":
                        weights = self.quad_scheme[rule][domain][n_points][0,:]
                        points  = self.quad_scheme[rule][domain][n_points][1:,:]
                    
                    else:                    
                        weights = self.quad_scheme[rule][domain][n_points][1,:]
                        points  = self.quad_scheme[rule][domain][n_points][0,:]
                    
                    return (weights, points)
                else:
                    
                    err_msg = "'{}' quadrature for a ".format(rule)
                    err_msg +="'{}' domain hasn't been defined for ".format(domain)
                    err_msg+= "{} points of integration!".format(n_points)
                    
                    raise QuadratureRuleError(err_msg)
                    
                    # print('Error: ', err_msg)
                    # sys.exit()
                    
            else:
                    
                    err_msg = "'{}' quadrature for a ".format(rule)
                    err_msg +="'{}' domain hasn't been defined!".format(domain)
                    
                    raise QuadratureRuleError(err_msg)
                    
                    # print('Error: ', err_msg)
                    # sys.exit()
        else:
                    
                    err_msg = "'{}' quadrature hasn't been defined!".format(rule)
                    
                    raise QuadratureRuleError(err_msg)
                    
                    # print('Error: ', err_msg)
                    # sys.exit()

def shape_functions(elementPoints, integrationPoint, NodesElement, d, ElType):
    """
    Return the shape function derivatives vector evalueted
    at the quadrature scheme integration point. \n
    The conversion from the [-1, 1] domain to the [x_i, x_j] domain
    is done internally in this function

    Parameters
    ----------
    elementPoints : Element coordinates xi, xj for a certain element \n
    integrationPoint : Shape function derivatives vector evaluation point in natural coordinates.\n
    NodesElement: number of nodes in the element.\n
    d: spatial dimension of the problem, e.g. 1 for 1D, 2 for 2D,...\n

    Returns
    -------
    dN: shape function vector

    """
        
    
    # Vector initialisation
    N  = np.zeros( [d, NodesElement] )
    dN = np.zeros( [d, NodesElement] )
    J  = np.zeros( [d,d] )
    
    
    
    if d == 1:
        
        xi  = integrationPoint
        
        if NodesElement == 2:
            # Shape functions evaluation
            N[0,0] = 0.5*(1-xi)
            N[0,1] = 0.5*(1+xi)
            
            # Shape functions derivatives evaluation
            dN[0,0] = -0.5
            dN[0,1] =  0.5
        else:
            # Shape functions evaluation
            N[0,0] = 0.5*xi*(xi - 1)
            N[0,1] = 1 - xi**2       
            N[0,2] = 0.5*xi*(xi + 1)
            
            # Shape functions derivatives evaluation
            dN[0,0] = xi - 0.5
            dN[0,1] = - 2*xi
            dN[0,2] = xi + 0.5
            
    elif d == 2:
        
        xi  = integrationPoint[0]
        eta = integrationPoint[1]
        
        if ElType == 'quad':
        
            N[0,0] = 0.25*(1-xi)*(1-eta)
            N[0,1] = 0.25*(1+xi)*(1-eta)
            N[0,2] = 0.25*(1+xi)*(1+eta)
            N[0,3] = 0.25*(1-xi)*(1+eta)
            N[1,0] = 0.25*(1-xi)*(1-eta)
            N[1,1] = 0.25*(1+xi)*(1-eta)
            N[1,2] = 0.25*(1+xi)*(1+eta)
            N[1,3] = 0.25*(1-xi)*(1+eta)
            
            dN[0,0] = -0.25*(1-eta)
            dN[0,1] = +0.25*(1-eta)
            dN[0,2] = +0.25*(1+eta)
            dN[0,3] = -0.25*(1+eta)
            dN[1,0] = -0.25*(1-xi)
            dN[1,1] = -0.25*(1+xi)
            dN[1,2] = +0.25*(1+xi)
            dN[1,3] = +0.25*(1-xi)
            
        elif ElType == 'triangle':
            
            N[0,0] = xi
            N[0,1] = eta
            N[0,2] = 1 - xi - eta
            N[1,:] = N[0,:]
            
            dN[0,0] =  1
            dN[0,1] =  0
            dN[0,2] = -1
            dN[1,0] =  0
            dN[1,1] =  1
            dN[1,2] = -1
        
        elif ElType == 'triangle6':
            N[0,0] = xi*(2*xi - 1)
            N[0,1] = eta*(2*eta - 1)
            N[0,2] = (1 - xi - eta)*(1 - 2*xi - 2*eta)
            N[0,3] = 4*xi*eta
            N[0,4] = 4*eta*(1 - xi - eta)
            N[0,5] = 4*xi*(1 - xi - eta)
            N[1,:] = N[0,:]
            
            dN[0,0] = 4*xi -1 
            dN[0,2] = 4*eta + 4*xi -3
            dN[0,3] = 4*eta
            dN[0,4] = -4*eta
            dN[0,5] = 4 -4*eta -8*xi
            dN[1,1] = 4*eta -1 
            dN[1,2] = 4*eta + 4*xi -3
            dN[1,3] = 4*xi
            dN[1,4] = 4 -4*xi -8*eta
            dN[1,5] = -4*xi
        else:    
            raise FemEngineError("Unknown element type!")
            
        
    else:    
        raise FemEngineError("Problem dimension mesh.d must be 1 or 2!")
    
    
    # Jacobian matrix evaluation
    J = np.dot(dN,elementPoints)

    if d <= 2:
        
        detJ = determinant(J)
        invJ = inverse(J)
    else:
        invJ = np.linalg.inv(J)
        detJ = np.linalg.det(J)
    
    if detJ < 0:
        raise FemEngineError("The Jacobian |J| < 0!")
    else:    
        
        # dN has to be expressed in terms of the global coordinate system (x,y,z)
        # B(x,y) = inv(J)*B(xi,eta)  
            
        dN = np.dot(invJ, dN)
        
        return (detJ, dN, N)
    

        
def determinant(matrix):
    """
    Compute the determinant of the given matrix.
    The matrix 

    Parameters
    ----------
    matrix 

    Returns
    -------
    det: Determinant

    """   
    det = 0
    n = np.shape(matrix)[1]
    m = np.shape(matrix)[1]
    
    if n == m:        
    
        if n == 1:
            return matrix[0,0]
        else:
            
            for i in range(n):
                minor = np.array( [ matrix[1:n, j] for j in range(m) if j != i] )
                det += (-1)**(i) * matrix[0,i] * determinant(minor)
            
            return det
            
    else:
        raise InversionError("The matrix isn't square!")
                

def inverse(matrix):
    """
    Compute the determinant of the given matrix.
    The matrix 

    Parameters
    ----------
    matrix 

    Returns
    -------
    det: Determinant

    """   
    det = 0
    n = np.shape(matrix)[1]
    m = np.shape(matrix)[1]
    
    if n == m:        
    
        if n == 1:
            return 1/matrix[0,0]
        elif n == 2:
            det = determinant(matrix)
            
            return (1/det) * np.array( [[ matrix[1,1], -matrix[0,1] ],
                                        [-matrix[1,0],  matrix[0,0] ]] )
            
            
        else:
            raise InversionError("Code for matrix bigger than 3x3 hasn't been developed yet!")
            
    else:
        raise InversionError("The matrix isn't square!")         
        
def rearrange(A, v1, v2):
    """
    Rearrange matrix A along rows following v1 indexes and\n
    along columns following v2

    Parameters
    ----------
    A : TYPE
        DESCRIPTION.
    v1 : TYPE
        DESCRIPTION.
    v2 : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    
    (m,n) = np.shape(A)
    
    for i in range(m):
        b = A[i, :]
        A[i,:] = b[v1]
    
    for j in range(n):
        b = A[:, j]
        A[:,j] = b[v2]
        
    return A

def element_table(mesh):
    """
    Returns the element table for the given mesh\n
    using an uniform MatSet and ElType

    Parameters
    ----------
    MatSet : Material set.\n
    ElType : Element type.\n
    mesh : mesh object.\n

    Returns
    -------
    Element table.

    """
    elem = np.zeros([mesh.Elements,mesh.NodesElement], dtype = int)     #= np.array([[5,8,0,1],
	        
    for i in range(mesh.Elements):
        for j in range(int(mesh.NodesElement)):
                  
            elem[i,j] = (mesh.NodesElement-1)*i+(j)
            
    return elem
    
    
def main():        
    print("------------------------------------------------------------------------")
    print("                              TESTING CODE")    
    print("------------------------------------------------------------------------\n")
    
    elem_coordinates = np.array([3, 6])
    NodesElement = 2
    d = 1
    eval_point = 4
    
    
    (detJ, dN, N) = shape_functions(elem_coordinates, eval_point, NodesElement, d)
    
    dN_exp = np.zeros([2,1])
    N_exp = np.zeros([2,1])
    
    print("shape_functions() results...\n")
    print("detJ:\n{}\n".format(detJ))
    print("dN:\n{}\n".format(dN))
    print("N:\n{}\n".format(N))
    
    L = elem_coordinates[-1] - elem_coordinates[0]
    
    x_val = eval_point*L*0.5 + 0.5*(elem_coordinates[-1] + elem_coordinates[0])
    
    N_exp[0] = 1 - x_val/L
    N_exp[-1] = x_val/L
    dN_exp[0] = -1/L
    dN_exp[-1] = 1/L
    detJ_exp = L/2
    
    print("Expected results...\n")
    print("detJ:\n{}\n".format(detJ_exp))
    print("dN:\n{}\n".format(dN_exp))
    print("N:\n{}\n".format(N_exp))
    
    A = np.array([[11,12,13,14,15,16],
     [21,22,23,24,25,26],
     [31,32,33,34,35,36],
     [41,42,43,44,45,46],
     [51,52,53,54,55,56],
     [61,62,63,64,65,66]])

    v1 = [5, 4, 3, 2, 1, 0]
    v2 = np.array(range(6))
    B = rearrange(A, v1, v2)
    print(B)

if __name__ == "__main__":
    main()


