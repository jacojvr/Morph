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
NC,Tri= readmsh('BSkull.msh')
NC2,Tri2 = readmsh('BSkull2.msh')
Tri2 = np.c_[Tri2[:,0],Tri2[:,2],Tri2[:,1]]
NC = np.ma.load('BSkullNCS_landm.txt')
NC2 = np.c_[-NC[:,0],NC[:,1:3]]

NCB1 = np.ma.load('ElasNodes_SymmetricOriginalSettings.txt')
NCB2 = np.ma.load('ElasNodes_Symmetric4stageBeta100.txt')
#NeigB = np.ma.load('BSkullNeig2layers.txt')
#NeigT = NeigB

##keep1 = NodeselectFP(NC,Tri,NeigB,np.array([200,1000]))
#keepT = np.ma.load('BaseKeepN_5_200_1000.txt')
##keep2 = NodeselectFP(NC,Tri,NeigB,np.array([50,100,1000]))
#keepB = np.ma.load('BaseKeepN_5_50_100_1000.txt')

#remove1=find_repeats(keepB)[0]
#remove2=find_repeats(keepT)[0]
#for i in remove1:
  #rows=np.where(keepB==i)[0]
  #keepB=np.r_[keepB[0:rows[0,],],keepB[rows[0,]+1:keepB.size],]
#for i in remove2:
  #rows=np.where(keepT==i)[0]
  #keepT=np.r_[keepT[0:rows[0,],],keepT[rows[0,]+1:keepT.size],]
  
##NC = NC*10
##NC2 = NC2*10

##NeigB = np.ma.load('1PrognNodeEigenval.txt')
##NeigT = NeigB

####lmkeep=np.where((NC2[LM2[:,0],0]<1)|(NC2[LM2[:,0],2]>-50))[0]
####LMB=LMB[lmkeep,]
####LM2=LM2[lmkeep,]
####node4SFP = np.where((NC2[:,0]<0)|(NC2[:,2]>-50))[0]
####nodes = np.where((NC2[keep2,0]<0)|(NC2[keep2,0]>-50))[0]
####keep2=keep2[nodes,]
####k=0
#SFPTT=np.array([range(0,Tri2.size/3)])#np.zeros((Tri2.size/3,))
####for i in range(0,tri2.size/3):
  ####if find_repeats(np.r_[tri2[i,].T,node4SFP])[0].size==3:
    ####SFPTT[k,]=i
    ####k=k+1
####SFPTT=SFPTT[range(0,k),]
####SFPTT=np.array(range(0,tri2.size/3))#SFPTT,int)
    
####NC = np.ma.load('ElasNodes29Jun.txt')
####NCE2 = np.ma.load('ElasticNodes14Jun14h30.txt')
##NCE=NCE*10
####NC2=NCE2*10
#beta = 150
#ELASNodes,CONV = elasticsurf(NC,Tri,LM1,NC2,Tri2,LM2,keepB,NeigB,keepT,NeigT,SFPTT,100,beta,np.array([range(1,8)]))#+range(11,16)]))
#np.ma.dump(ELASNodes,'ElasNodes_Symmetric4stageBeta'+str(beta)+'.txt')
#np.ma.dump(CONV,'ElasConv_Symmetric4stageBeta'+str(beta)+'.txt')

figure(fgcolor=(0,0,0),bgcolor=(1,1,1))
basE = triangular_mesh(NC[:,0],NC[:,2],NC[:,1],Tri,color=(1,1,0),representation='wireframe',opacity=0)
targeT = triangular_mesh(NC2[:,0],NC2[:,2],NC2[:,1],Tri2,color=(1,0,0),representation='wireframe',opacity=0)
plt0 = triangular_mesh(NCB1[:,0],NCB1[:,2],NCB1[:,1],Tri,color=(0,0,1),representation='wireframe',opacity=0)
plt1 = triangular_mesh(NCB2[:,0],NCB2[:,2],NCB2[:,1],Tri,color=(0,0,1),representation='wireframe',opacity=0)

cutpT = pipeline.scalar_cut_plane(targeT,color=(1,0,0),line_width=5,plane_orientation='y_axes',name='Target')
cutpT.implicit_plane.origin = (0,20,-40)
cutpT.implicit_plane.widget.enabled = False
cutpB = pipeline.scalar_cut_plane(basE,color=(0,0,0),line_width=3,plane_orientation='y_axes',name='Base')
cutpB.implicit_plane.origin = (0,20,-40)
cutpB.implicit_plane.widget.enabled = False
cutp0 = pipeline.scalar_cut_plane(plt0,color=(0,0,0),line_width=3,plane_orientation='y_axes',name='iter_00')
cutp0.implicit_plane.origin = (0,20,-40)#(20,0,0)
cutp0.implicit_plane.widget.enabled = False
cutp1 = pipeline.scalar_cut_plane(plt1,color=(0,0,0),line_width=3,plane_orientation='y_axes',name='iter_01')
cutp1.implicit_plane.origin = (0,20,-40)#(20,0,0)
cutp1.implicit_plane.widget.enabled = False

##triangular_mesh(NCE2[:,0],NCE2[:,2],NCE2[:,1],Tri,color=(0,0,1),representation='wireframe',opacity=0.2)
#points3d(NCE[LM1,0],NCE[LM1,2],NCE[LM1,1],color=(0,0,0),scale_factor=5)
#points3d(NC2[LM2,0],NC2[LM2,2],NC2[LM2,1],color=(0,1,0),scale_factor=5)
