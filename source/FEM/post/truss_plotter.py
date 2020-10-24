#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 29 23:01:17 2020

@author: angelo.simone@unipd.it
"""

# Shamelessy copied (and heavily modified) from https://github.com/nasseralkmim/sapy/blob/master/sapy/plotter.py
# code retrieved 29 April 2020

# TODO: too many repetions?

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

import FEM.FEM_engine as FEM_engine

# -----------------------------------------------------------------------------

def undeformed_mesh(mesh, U, MaterialSets, Procedures):
    points = mesh.points
    _plot_mesh(mesh,points,"undeformed")    

# -----------------------------------------------------------------------------

def deformed_mesh(mesh, U, MaterialSets, Procedures):
    magnification = float(Procedures["postprocessing"]["deformed mesh"]["magnification factor"])
    points = _add_displacement(mesh, U, magnification)
    _plot_mesh(mesh,points,"deformed")    

# -----------------------------------------------------------------------------

def axial_stress(mesh, U, MaterialSets, Procedures):
    points = mesh.points
    Stress = []
    for e in range(len(mesh.elements)):
        stress = FEM_engine.stress_bar(e,mesh,MaterialSets,U)
        Stress.append(stress)
    _plot_mesh(mesh,points,"axial stress",axial_stress = Stress)

# -----------------------------------------------------------------------------

def _window(name):
    return plt.figure(name)

# -----------------------------------------------------------------------------

def _add_displacement(mesh, U, magnification):
    points = np.copy(mesh.points)
    for n in range(len(mesh.points)):
        for d in range(mesh.d): # mesh.d = spatial dimension -> number of displacement dofs
            dof = n*mesh.dofsNode + d
            points[n, d] += U[dof]*magnification
    return points            

# -----------------------------------------------------------------------------

def _plot_mesh(mesh, points, title, **kwargs):
    header = _window(title)
    axes = header.add_subplot(111, aspect='equal')
    
    if title == "undeformed" or title == "axial stress":
        _draw(points, mesh.elements, axes, color='black')
    elif title == "deformed":
        _draw(mesh.points, mesh.elements, axes, color='black')
        _draw(points, mesh.elements, axes, color='red')
    
    _add_node_label(points, axes)    
    _add_element_label(points, mesh.elements, axes, **kwargs)
    header.tight_layout()

    # return axes

# -----------------------------------------------------------------------------

def _draw(points, elements, axes, color):
    axes.set_xlabel('x')
    axes.set_ylabel('y')

    # draw nodes
    for p in points:
        axes.scatter(p[0], p[1], c='k', alpha=0.3, marker='s')

    # draw edges
    for elementNodes in elements:
        x1 = points[elementNodes[0]][0]
        y1 = points[elementNodes[0]][1] 
        x2 = points[elementNodes[1]][0]
        y2 = points[elementNodes[1]][1]
        xs = [x1, x2]
        ys = [y1, y2]
        edge = Line2D(xs, ys, linewidth=1.0, color=color)
        axes.add_line(edge)

# -----------------------------------------------------------------------------

def _add_node_label(points, axes):
    for nodeNumber, nodeCoordinates in enumerate(points):
        axes.text(nodeCoordinates[0], nodeCoordinates[1], str(nodeNumber), color='b', size=10)

# -----------------------------------------------------------------------------

def _add_element_label(points, elements, axes, **kwargs):    
    for elementNumber, elementNodes in enumerate(elements):
        x1 = points[elementNodes[0]][0]
        y1 = points[elementNodes[0]][1] 
        x2 = points[elementNodes[1]][0]
        y2 = points[elementNodes[1]][1]
        t=0.4 # label position along normalized bar length
        x = (1-t)*x1 + t*x2
        y = (1-t)*y1 + t*y2

        if("axial_stress" in kwargs):
            Q=kwargs['axial_stress']
            if Q[elementNumber] < 0:
                axes.text(x, y, str(np.round_(Q[elementNumber], 1)), color='r', size=10)
            else:
                axes.text(x, y, str(np.round_(Q[elementNumber], 1)), color='g', size=10)
        else:
            axes.text(x, y, str(elementNumber), color='g', size=10)
       
# -----------------------------------------------------------------------------

# NOTE:

# - Function names starting with a signle underscore are supposed to be used 
#   only in this module, but you can still call them from other modules.
#   To quote PEP-8:
#   _single_leading_underscore: weak "internal use" indicator. 
#   E.g. from M import * does not import objects whose name starts 
#   with an underscore.

# - We make use of optional parameters by means of keyword arguments (kwargs)
# https://stackoverflow.com/questions/1769403/what-is-the-purpose-and-use-of-kwargs
# https://stackoverflow.com/questions/43279256/multiple-optional-arguments-python
# https://www.earthdatascience.org/courses/intro-to-earth-data-science/write-efficient-python-code/functions-modular-code/write-functions-with-multiple-and-optional-parameters-in-python/
# https://scipy-lectures.org/intro/language/functions.html (1.2.4.6. Variable number of parameters)
