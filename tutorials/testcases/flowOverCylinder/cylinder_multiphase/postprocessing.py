#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 10 11:23:16 2015

@author: agross
"""

import numpy as np

"""
import matplotlib as mpl
mpl.use('pgf')

from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

def figsize(scale):
    fig_width_pt = 469.755                          # Get this from LaTeX using \the\textwidth
    inches_per_pt = 1.0/72.27                       # Convert pt to inch
    golden_mean = (np.sqrt(5.0)-1.0)/2.0            # Aesthetic ratio (you could change this)
    fig_width = fig_width_pt*inches_per_pt*scale    # width in inches
    fig_height = fig_width*golden_mean              # height in inches
    fig_size = [fig_width,fig_height]
    return fig_size

pgf_with_latex = {                      # setup matplotlib to use latex for output
    "pgf.texsystem": "pdflatex",        # change this if using xetex or lautex
    "text.usetex": True,                # use LaTeX to write all text
    "font.family": "serif",
    "font.serif": [],                   # blank entries should cause plots to inherit fonts from the document
    "font.sans-serif": [], 
    "font.monospace": [], 
    "axes.labelsize": 10,               # LaTeX default is 10pt font.
    "text.fontsize": 10, 
    "legend.fontsize": 8,               # Make the legend/label fonts a little smaller
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "figure.figsize": figsize(0.9),     # default fig size of 0.9 textwidth
    "pgf.preamble": [
        r"\usepackage[utf8x]{inputenc}",    # use utf8 fonts becasue your computer can handle it :)
        r"\usepackage[T1]{fontenc}",        # plots will be generated using this preamble
        r"\usepackage{pgfplots}",           # force pgf plot load
        r"\usepgfplotslibrary{external}",   # load externalization fpr pgf plots
        r"\tikzexternalize"                 # load externalization fpr pgf plots
        ]
    }
mpl.rcParams.update(pgf_with_latex)


# I make my own newfig and savefig functions
def newfig(width):
    plt.close('all')
    plt.clf()
    fig = plt.figure(figsize=figsize(width))
    ax = fig.add_subplot(111)
    return fig, ax

def savefig(filename):
    print("saving "+str(filename))
    #plt.savefig('{}.pgf'.format(filename))
    plt.savefig('{}.pdf'.format(filename.replace(".","-")))
"""

import matplotlib.pyplot as plt
from vtkToNumpy import vtkToNumpy

timesteps = ["33.324"]

U=2.63
rho=1000

fs=0.6

norm=np.pi/(0.5*rho*U*U)

for timestep in timesteps:
    data = []
    f = open("postProcessing/sets/"+timestep+"/cylinderPoints_pMean_p_rgh_alpha.waterMean_patchMagSf.xy")
    for l in f:
        data.append(map(float,l.split()))
    data = np.array(data)
    #data = np.genfromtxt("postProcessing/sets/"+timestep+"/cylinderPoints_pMean_alpha.waterMean.xy")
    sorting = np.argsort(data[:,2])
    coords = data[:,0:3][sorting]
    p = data[:,3][sorting]
    p_rgh = data[:,4][sorting]
    alpha = data[:,5][sorting]
    magSf = data[:,6][sorting]
    #data = np.genfromtxt("postProcessing/sets/"+timestep+"/cylinderPoints_wallShearStress_patchSf.xy")
    data=[]
    f = open("postProcessing/sets/"+timestep+"/cylinderPoints_wallShearStress_patchSf.xy")
    for l in f:
        data.append(map(float,l.split()))
    data = np.array(data)
    sorting = np.argsort(data[:,2])
    wallShearStress = data[:,3:6][sorting]
    Sf = data[:,6:9][sorting]
    old_z = -999
    Ff = []
    Fp = []
    A = []
    z = []
    for i in range(len(p)):
        delta=coords[i,2]-old_z
        deltaFf = -wallShearStress[i,:]*alpha[i]*rho*magSf[i]
        deltaFp = Sf[i,:]*p_rgh[i]
        if np.abs(delta) > 1e-5:
            old_z=coords[i,2]
            z.append(old_z)
            A.append(magSf[i])
            Ff.append(deltaFf)
            Fp.append(deltaFp)
        else:
            Ff[-1] += deltaFf
            Fp[-1] += deltaFp
            A[-1] += magSf[i]
    Ff=np.array(Ff)
    Fp=np.array(Fp)
    A=np.array(A)
    cf = np.array(Ff*norm/np.array([A,A,A]).T)
    cp = np.array(Fp*norm/np.array([A,A,A]).T)
    cd = cf+cp

    plt.plot(cf[:,0],z,label="Friction")
    plt.plot(cp[:,0],z,label="Pressure")
    plt.plot(cd[:,0],z,label="Drag")
    plt.show()

    vtkName = "postProcessing/surfaces/"+timestep+"/freeSurfaceMean.vtk"
    data = vtkToNumpy(vtkName)

    data.boxclip([-4, 0,-999],[4,4,999])

    x,y,z = data.clipped_nodes[:,0],data.clipped_nodes[:,1],data.clipped_nodes[:,2]

    clipped_tris = data.clipped_triangles

    delta_z = 0.03
    n_zmin = int(np.floor(np.min(z)/delta_z))
    n_zmax = int(np.ceil(np.max(z)/delta_z))
    
    colorLevels = np.linspace(n_zmin*delta_z,n_zmax*delta_z,n_zmax-n_zmin+1)

#    ticks = np.linspace(minF,maxF,nTicks,endpoint=False)
    
#    ax.tricontour(x,y,trianglesn,clipped_field_np,colorLevels,linewidths=0.5,colors="k")
#    cntr = ax.tricontourf(x,y,trianglesn,clipped_field_np,colorLevels)

    #plt.tripcolor(x,y,clipped_tris,z)
    CS = plt.tricontour(x,y,clipped_tris,z,colorLevels,colors='k')
    plt.tricontourf(x,y,clipped_tris,z,colorLevels,cm=plt.cm.jet)
    plt.clabel(CS, colorLevels[n_zmin%4+1::4],  # label every second level
           inline=1,
           fmt='%1.2f', colors='k')
    plt.xlim((-4,4))
    plt.ylim((0,4))
    plt.axes().set_aspect('equal')
    plt.show()
    #fig, ax  = newfig(fs)
    #ax.set_ylabel("$y$")
    #ax.set_xlabel("$z$")
    #cmap=mpl.cm.jet
    #im=ax.tripcolor(cy,cz,clipped_tris,clipped_p,cmap=cmap)
    #fig.colorbar(im)
    #ax.plot((centre_p[1]),(centre_p[2]),marker,color='black')
    #ax.plot((centre_U[1]),(centre_U[2]),marker,color='white')
    #ax.plot((centre_LO[1]),(centre_LO[2]),marker,color='#808080')
    #ax.set_aspect('equal', 'datalim')
    #savefig(result_prefix+section+"_pressure")
    
#plt.plot(magSf)
#plt.show()

