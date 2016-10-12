#/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 11 10:13:00 2015

@author: agross
"""

import warnings
import numpy as np
import vtk
import vtk.util.numpy_support

class vtkToNumpy:
    """
    reads vtk file to numpy format for plotting
    """
    def __init__(self, filename=None):
        """
        initialize new instance with file name (optimal)
        """
        self.filename=filename
        if isinstance(filename,str):
            self.update()    
        else:
            self.clear()

    def filename(self, newname=False):
        """
        set file name(optinal) and return file name
        """
        if isinstance(newname, str):
            self.clear()
            self.filename = newname
        return self.filename

    def clear(self):
        """
        clear all data
        """
        self.reader = None
        self.field_names = []
        self.triangles_np = np.array([])
        self.nodes_np = np.array([])
        self.fields_np = np.asarray([[]])
        self.clipped_nodes = np.array([])
        self.clipped_triangles = np.array([])
        self.clipped_fields = np.asarray([[]])

    def update(self):
        """
        re-read data file and convert
        """
        self.clear()
        self.reader=vtk.vtkPolyDataReader()
        self.reader.SetFileName(self.filename)
        self.reader.Update()
        output=self.reader.GetOutput()
        
        # Get the coordinates of nodes in the triangulated mesh
        nodes_vtk_array= output.GetPoints().GetData()
        
        # Get the scalar field to be plotted
        self.field_names = [output.GetPointData().GetArrayName(i) for i in range(output.GetPointData().GetNumberOfArrays())]
        field_vtk_arrays = [output.GetPointData().GetArray(n) for n in self.field_names]
        
        # Get cells (= triangles, other types currently not supported)
        polys = output.GetPolys()
        
        #Convert vtkPolys to [n1,n2,n3] - numpy format and create node to tri map
        idList = vtk.vtkIdList()
        self.triangles_np=[]
        #triToNode=[None]*polys.GetNumberOfCells()
        for i in range(polys.GetNumberOfCells()):
            polys.GetNextCell(idList)
            if int(idList.GetNumberOfIds()) is not 3:
                warnings.warn("Element is not a triangle, ignoring.", SyntaxWarning, stacklevel=2)
            else:
                self.triangles_np.append([idList.GetId(i) for i in range(3)])
        self.triangles_np = np.asarray(self.triangles_np)
        
        #Get the coordinates of the nodes and their scalar values
        self.nodes_np = vtk.util.numpy_support.vtk_to_numpy(nodes_vtk_array)
        self.fields_np = [vtk.util.numpy_support.vtk_to_numpy(field_vtk_arrays[i]) for i in range(len(field_vtk_arrays))]

    def boxclip(self, xmin, xmax):
        """
        clip fields geometrically in box
        """
        if len(xmin) is not 3 or len(xmax) is not 3:
            raise ValueError("boxclip bound must be of size 3")
        x_min, y_min, z_min = map(float,xmin)
        x_max, y_max, z_max = map(float,xmax)
        mask=np.zeros(len(self.nodes_np))
        for idx, m in enumerate(mask):
            if ( self.nodes_np[idx,0]>x_min and self.nodes_np[idx,0]<x_max
            and self.nodes_np[idx,1]>y_min and self.nodes_np[idx,1]<y_max
            and self.nodes_np[idx,2]>z_min and self.nodes_np[idx,2]<z_max ):
                mask[idx] = 1.0
        
        self.clip(mask)
        return mask

    def sphereclip(self, centre, radius):
        """
        clip fields geometrically in box
        """
        if len(centre) is not 3:
            raise ValueError("sphereclip centre must be of size 3")
        x, y, z = map(float,centre)
        rsquare = np.square(float(radius))
        mask=np.zeros(len(self.nodes_np))
        for idx, m in enumerate(mask):
            distsquare = ( (self.nodes_np[idx,0]-x)**2.0 +
                           (self.nodes_np[idx,1]-y)**2.0 +
                           (self.nodes_np[idx,2]-z)**2.0 )
            if not distsquare-rsquare > 0.0:
                mask[idx] = 1.0
        self.clip(mask)
        return mask

    def clip(self, mask):
        """
        clip fields geometrically with clipmask
        all triangles with at least one node with a mask value larger then 0.5
        will be preseved
        """
        if len(mask)-len(self.nodes_np) != 0:
            raise ValueError("clip map mask ("+str(len(mask))+") must be of same size as nodes ("+str(len(self.nodes_np))+")")

        # generate new indices for nodes
        trianglemask=[0]*len(self.triangles_np)
        extended_mask=np.zeros(len(mask))
        mask=map(float,mask)
        for i, m in enumerate(trianglemask):
            for j in self.triangles_np[i]:
                if mask[j]>0.5:
                    trianglemask[i]=1.0
            if trianglemask[i] > 0.5:
                for j in self.triangles_np[i]:
                    extended_mask[j]=1.0
        
        # map old node numbers to new and create new points
        mapOldToNew = -1.0*np.ones(len(mask))
        self.clipped_nodes = []
        num=0
        self.clipped_fields = [[] for _ in range(len(self.fields_np))]
    
        for i, m in enumerate(mapOldToNew):
            if extended_mask[i] > 0.5:
                mapOldToNew[i]=num
                num+=1
                self.clipped_nodes.append(list(self.nodes_np[i]))
                for f in range(len(self.fields_np)):
                    self.clipped_fields[f].append(self.fields_np[f][i])
        
        self.clipped_nodes=np.array(self.clipped_nodes)
        self.clipped_triangles = []
        for i, m in enumerate(trianglemask):
            if m > 0.5:
                oldTri=self.triangles_np[i]
                tri = [mapOldToNew[j] for j in oldTri]
                if -1.0 in tri:
                    ValueError("trying to map vertices of tri "+str(i)+" which is not in mapped")
                self.clipped_triangles.append(tri)
                
        self.clipped_triangles = np.asarray(self.clipped_triangles)
        self.clipped_fields= [np.asarray(self.clipped_fields[i]) for i in range(len(self.clipped_fields))]