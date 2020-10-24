# -*- coding: utf-8 -*-
"""
Created on Sat May 16 19:31:37 2020

@author: Umberto Bettinardi
         bettinardi96@gmail.com
"""


import matplotlib.pyplot as plt
import numpy as np

def plot_points(mesh,U,n, mag, show = False):
    U_print = U.copy().reshape(len(mesh.points),mesh.d)
    old_points = mesh.points.copy()
    new_points = old_points.copy()
    
    for i in range(len(new_points)):
        
        new_points[i,0] += U_print[i,0]*mag
        new_points[i,1] += U_print[i,1]*mag
    
    plt.figure(n)
    plt.title("Plate with ciruclar hole under\n"+r"uniform tension $\sigma$")
    plt.xlabel('x [cm]')
    plt.ylabel('y [cm]')
    plt.axis('equal')
    plt.grid(which = 'both')
    plt.plot(old_points[:,0],old_points[:,1],'+b', markersize = 0.6)
    plt.plot(new_points[:,0],new_points[:,1],'+r', markersize = 0.6)
    plt.savefig('piastra.png', dpi = 1200)
    plt.legend(['Undeformed', 'Deformed'])
    
    if show == True:
        plt.show()