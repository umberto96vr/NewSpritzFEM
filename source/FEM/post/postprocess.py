#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 14:44:21 2020

@author: angelo.simone@unipd.it
"""

import matplotlib.pyplot as plt
import sys

# specify where to look for modules
sys.path.append("../../source")

import FEM.post.truss_plotter as tp

def run(mesh, U, MaterialSets, Procedures):
    
    keys = Procedures["postprocessing"]

    # In-place replacement of white space in key with underscore 
    # for all the keys in parameters["postprocessing"]
    # This enables the use of keywords with spaces to be translated into function names: axial stress -> axial_stress
    # https://www.geeksforgeeks.org/python-remove-spaces-from-dictionary-keys/
    keys = {x.replace(' ', '_'): v
               for x, v in keys.items()} 

    for key in keys:
        exec ("tp."+key+"(mesh, U, MaterialSets, Procedures)") 

    plt.show() #TODO: this is apparently not need in spyder but necessary if code
               #      is executed from ipython ...but in that case there is an error.
