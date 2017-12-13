# -*- coding: utf-8 -*-
import numpy as np                                                                        
from enthought.mayavi.mlab import *
from fileread import *
from pmeshdeform import *
from scipy.io import *
from scipy.interpolate import Rbf
from scipy.spatial import KDTree
from scipy.stats.mstats import find_repeats
from ptrimeshalt import *
from rotatescaletrans import *
import time
LIM = 8
NC,Tri,LM1 = readmsh('BSkull.msh')
NC2,Tri2,LM2 = readmsh('BSkull2.msh')
NC = np.ma.load('BSkullNCS_landm.txt')
Tri2 = np.c_[Tri2[:,0],Tri2[:,2],Tri2[:,1]]
NC2 = np.c_[-NC[:,0],NC[:,1:3]]
NCE = np.ma.load('ElasNodes_SymmetricOriginalSettings.txt')
CONV = np.ma.load('ElasConv_SymmetricOriginalSettings.txt')
NCE1 = np.ma.load('ElasNodes_SymmetricAllBeta0.txt')
CONV1 = np.ma.load('ElasConv_SymmetricAllBeta0.txt')
NCE2 = np.ma.load('ElasNodes_Symmetric2stage1Beta100kmax5.txt')
CONV2 = np.ma.load('ElasConv_Symmetric2stage1Beta100kmax5.txt')
NCE3 = np.ma.load('ElasNodes_Symmetric4stageBeta100.txt')
CONV3 = np.ma.load('ElasConv_Symmetric4stageBeta100.txt')
NCE4 = np.ma.load('TempElasNodes_Iter11_TimeWed Jul 14 12:57:41 2010.txt')
#NCE2 = np.ma.load('ElasNodes1Jul.txt')
#NCE3 = np.ma.load('ElasNodes12Jul.txt')


NeigB = np.ma.load('BSkullNeig2layers.txt')
NeigT = NeigB


figure(fgcolor=(0,0,0),bgcolor=(1,1,1))
#base = triangular_mesh(NC[:,0],NC[:,2],NC[:,1],Tri,color=(1,0,0),opacity=0)
targ = triangular_mesh(NC2[:,0],NC2[:,2],NC2[:,1],Tri2,color=(0,0,0),opacity=0)
#elas = triangular_mesh(NCE[:,0],NCE[:,2],NCE[:,1],Tri,color=(1,0,0),opacity=0)
#elas1 = triangular_mesh(NCE1[:,0],NCE1[:,2],NCE1[:,1],Tri,color=(1,0,0),opacity=0)
#elas2 = triangular_mesh(NCE2[:,0],NCE2[:,2],NCE2[:,1],Tri,color=(1,0,0),opacity=0)
#elas3 = triangular_mesh(NCE3[:,0],NCE3[:,2],NCE3[:,1],Tri,color=(1,0,0),opacity=0)
elas4 = triangular_mesh(NCE4[:,0],NCE4[:,2],NCE4[:,1],Tri,color=(1,0,0),opacity=0)
#cutbase = pipeline.scalar_cut_plane(base,plane_orientation='y_axes',color=(1,0,0))
cuttarg = pipeline.scalar_cut_plane(targ,plane_orientation='y_axes',color=(0,0,0),name='Target')
#cutelas = pipeline.scalar_cut_plane(elas,plane_orientation='y_axes',color=(1,0,0),name='OriginalSettings')
#cutelas1 = pipeline.scalar_cut_plane(elas1,plane_orientation='y_axes',color=(1,0,0),name='OriginalCourse')
#cutelas2 = pipeline.scalar_cut_plane(elas2,plane_orientation='y_axes',color=(1,0,0),name='Stage1Beta100kmax5')
#cutelas3 = pipeline.scalar_cut_plane(elas3,plane_orientation='y_axes',color=(1,0,0),name='Stage2Beta100kmax100')
cutelas4 = pipeline.scalar_cut_plane(elas4,plane_orientation='y_axes',color=(1,0,0),name='Temp@iter10Stage2Beta150kmax100')
#triangular_mesh(NCE2[:,0],NCE2[:,2],NCE2[:,1],Tri,color=(0,1,0),opacity=0.2)
#triangular_mesh(NCE3[:,0],NCE3[:,2],NCE3[:,1],Tri,color=(0,1,1),opacity=0.2)
#points3d(NC[LM1,0],NC[LM1,2],NC[LM1,1],color=(0,0,0),scale_factor=5)
#points3d(NC2[LM2,0],NC2[LM2,2],NC2[LM2,1],color=(0,1,0),scale_factor=5)