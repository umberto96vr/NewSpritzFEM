"""
Created on Fri Mar 20 14:45:31 2020

@author: Umberto Bettinardi
         bettinardi96@gmail.com
"""
import os
import sys

import meshio

class MeshEngineError(Exception):
    def __init__(self, text):
        self.text = text
        Exception.__init__(self, text)
        
#---------------------------------------------------------------------------------------------

class Mesh():
    """ Contains information about the discretization.\n
        Parameters: 
        ----------
        mesh.d
        \t This scalar defines the spatial dimension (d = 1, 2, 3). This attribute 
        \t should be defined by the user because GMSH works with a xyz coordinate 
        \t system (the spatial dimension cannot be retrieved from the coordinate array).
        mesh.dofsNode  
        \tThis scalar defined the number of degrees of freedom per 
        \tnode and is assumed to be the same for all nodes. This attribute 
        \tshould be defined by the user because in a coupled problem 
        \tthe number of dofs per node cannot be retrieved from the spatial dimension.
        mesh.elements 
        \tThe elements are defined in terms of their node numbers in the 
        \telement connectivity array. Each row of this array defines an element.
        \tElement nodes are ordered in a counter-clockwise fashion.
        mesh.elementMaterialTag
        \t[list]\n
        \tEach element is assigned to a material
        mesh.elementType
        \t[list of strings] \n
        \tThis list of strings defines the element type.\n
        \tAvailable elements are "bar", "bar2d", "diff_bar", "bar_pullout", "triangle", "quad".
        mesh.points
        \t[numpy array] \n
        \tCoordinates in xyz format of all mesh points (= nodes).

        """
        
    #------------------------------------------------------------------------------------------
    
    def __init__(self):
        
        self.elementMaterialTag = None 
        self.elementType        = None
        self.Nodes              = None
        self.dofsNode           = None
        self.elements           = None 
        self.points             = None
        self.d                  = None
      
       
    def __str__(self):
        
        str_info = "------------------------------------------\n"
        str_info += "             MESH INFO\n"
        str_info += "------------------------------------------\n"
        str_info += "Number of elements:          {}\n".format(self.Elements)
        str_info += "Number of nodes:             {}\n".format(self.Nodes)
        str_info += "Number of nodes per elment:  {}\n".format(self.NodesElement)
        str_info += "Number of dofs per node:     {}\n\n".format(self.dofsNode)
        str_info += "Connectivity table:\n\n"
        str_info += "| Material |Element type  |   Node i    |Node j\n"
        for i in range(self.Elements):
            
            str_info += "|   {}      |        {}     |     {}       |     {}\n".format(
                self.elements[i,0],self.elements[i,1],self.elements[i,2],self.elements[i,3])
        return str_info


def get_key(my_dict,val): 
    """ Function to return key for any value. """
    
    # This function returns the key if the first item in the array value 
    # of a dictionary is equal to val. If my_dict contains 
    # 'Fixed': array([667,   0]), get_key(my_dict,667) returns Fixed
    
    for key, value in my_dict.items(): 
         if val == value[0]: 
             return key 
    
    print("\n value",val,"doesn't exist as \'key\': array([value, 0]) in\n", my_dict)
    sys.exit()
  
# -----------------------------------------------------------------------------

def GMSH(mesh_file, rewrite = True):

    # TODO: make check on spatial dimension and element type
    # TODO: make check on GMSH file versions
    # TODO: generalize to hybrid meshes with quad and triangle
    #
    # Known limitations: 
    # This function has been tested for 3-node triangular elements, and 
    # it could work for four-node quadrilateral elements. It does not have
    # checks for one- or three-dimensional meshes; it does not work with hybrid meshes (e.g., triangular and quadrilateral elements)
    
    # remove existing msh file and doesn't complain if the file is absent 
    if rewrite == True:
        try:
            os.remove(mesh_file+".msh")
        except OSError:
            pass
        
        # run GMSH to generate the mesh
        os.system("gmsh -2 " + mesh_file + ".geo" + " -o " + mesh_file + ".msh")
        
        # create a mesh object
        mesh = meshio.read(mesh_file+".msh")
        
    else:
        # create a mesh object
        mesh = meshio.read(mesh_file+".msh")
    
    # check if the mesh object contains attributes needed by pyFEM
    # - pyFEM_MeshAttributes is a list of all the mesh attributes needed by pyFEM
    # - we are going to reuse the attribute points and add the other attribute from pyFEM_MeshAttributes
    pyFEM_MeshAttributes = ["d", "dofsNode", "elements", "elementMaterialTag", "elementType", "points"]

    for attribute in pyFEM_MeshAttributes:
        if attribute in dir(mesh):
            if attribute == "points":
                pass
            else:
                print("Error: meshio already contains the attribute",attribute)
                print("       ...do something!")
                sys.exit()

    # add the missing attributes from pyFEM_MeshAttributes
    
    # Note: it is assumed that the mesh is two-dimensional and that the
    # domain is dicretized with triangular elements and that there are
    # two degrees of freedom per node (i.e., this is a plain equilibrium problem)
    
    mesh.d = 2 
    mesh.dofsNode = 2 
    mesh.elements = []
    mesh.elementMaterialTag = []
    mesh.elementType = []
    meshing = False
    
    quad = False
    try:
        dummy = mesh.cell_data_dict['gmsh:physical']['quad']
        quad = True
    except KeyError:
        pass#raise MeshEngineError("No quadrilateral elements in mesh")
    
    triangle = False
    try:
        dummy = mesh.cell_data_dict['gmsh:physical']['triangle']
        triangle = True
    except KeyError:
        pass#raise MeshEngineError("No triangular elements in mesh")
        
    triangle6 = False
    try:
        dummy = mesh.cell_data_dict['gmsh:physical']['triangle6']
        triangle6 = True
    except KeyError:
        pass#raise MeshEngineError("No 6-node triangular elements in mesh")

    if quad:
        meshing = True
        quads = len(mesh.cell_data_dict["gmsh:physical"]["quad"])
        for q in range(quads):
            mesh.elements.append(mesh.cells_dict["quad"][q])
            materialTag=mesh.cell_data_dict["gmsh:physical"]["quad"][q]           
            key = get_key(mesh.field_data, materialTag)
            mesh.elementMaterialTag.append(key)            
            mesh.elementType.append("quad")
        
    if triangle: 
        meshing = True
        triangles = len(mesh.cell_data_dict["gmsh:physical"]["triangle"])
        for t in range(triangles):
            mesh.elements.append(mesh.cells_dict["triangle"][t])
            materialTag=mesh.cell_data_dict["gmsh:physical"]["triangle"][t]          
            key = get_key(mesh.field_data, materialTag)
            mesh.elementMaterialTag.append(key)            
            mesh.elementType.append("triangle")
    
    if triangle6: 
        meshing = True
        triangles6 = len(mesh.cell_data_dict["gmsh:physical"]["triangle6"])
        for t in range(triangles6):
            mesh.elements.append(mesh.cells_dict["triangle6"][t])
            materialTag=mesh.cell_data_dict["gmsh:physical"]["triangle6"][t]          
            key = get_key(mesh.field_data, materialTag)
            mesh.elementMaterialTag.append(key)            
            mesh.elementType.append("triangle6")
    
    if not meshing:
        raise MeshEngineError("something went wrong: could not extract mesh data")
        sys.exit()
        
    # TODO: ...check that all the necessary attributes have been defined in a correct manner

    return mesh