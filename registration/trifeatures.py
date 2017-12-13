# -*- coding: utf-8 -*-

import numpy as np
import scipy as sp
from scipy.stats.mstats import find_repeats
import pprocess
from ptrimeshalt import *

LIM=8

def crestlines(NC,Tri):
  N1 = NC.shape[0]
  NPP1 = N1/LIM
  T1 = Tri.shape[0]
  TPP1 = T1/LIM
  # Determine Triangle unit normal
  TBC = np.c_[np.sum(np.c_[NC[Tri[:,0],0],NC[Tri[:,1],0],NC[Tri[:,2],0]],1)/3,
    np.sum(np.c_[NC[Tri[:,0],1],NC[Tri[:,1],1],NC[Tri[:,2],1]],1)/3,np.sum(np.c_[NC[Tri[:,0],2],NC[Tri[:,1],2],NC[Tri[:,2],2]],1)/3]
  TNORM = np.cross(NC[Tri[:,1],:]-NC[Tri[:,0],:],NC[Tri[:,2],:]-NC[Tri[:,0],:])
  TNORM = (TNORM.T/(np.ones((3,1))*np.sqrt(np.sum(np.array([TNORM*TNORM]),2)))).T
  VNORM = vrtxnormal(NC,Tri,TBC,TNORM)
  
  # Determine the principal curvature directions at each node
  PDIRECT=np.zeros(NC.shape)
  for ni in range(0,N1):
    PDIRECT[i,]=principalDirection(ni,NC,Tri,VNORM)
  #results = pprocess.Map(limit=LIM)
  #calc = results.manage(pprocess.MakeParallel(principalDirection))
  #for j in range(0,LIM):                                                                                                                                        
    #calc(np.array(range(0,NPP1))+j*NPP1,NC,Tri,VNORM)
  #for j in range(0,LIM):
    #PDIRECT[np.array(range(0,NPP1))+j*NPP1,:] = results[j]
  #PDIRECT[range(LIM*NPP1,N1),:]=principalDirection(range(LIM*NPP1,N1),NC,Tri,VNORM)


def principalDirection(ni,NC,Tri,VNORM):
  # Determine principal curvature direction for node "ni" on unstructured surface mesh.
  neighb1 = Nneighbours(ni,NC,Tri,1)[1]
  PlaneN = np.sum(VNORM[neighb1,:],0)
  PlaneN = PlaneN/(np.ones((3,))*np.sqrt(np.sum(PlaneN*PlaneN)))
  DistNeigh = np.ones((neighb1.size,1))*NC[ni,:]-NC[neighb1,:]
  VertDist = (np.ones((neighb1.size,1))*PlaneN)*DistNeigh
  NCPlaneNeighb = NC[neighb3,:]+VertDist*PlaneN
  NC = NCPlaneNeighb
  DistNeigh = np.ones((neighb1.size,1))*NC[ni,:]-NC[neighb1,:]
  # From Meyer et al., determine the mean curvature normal operator (Laplace-Beltrami) K(x_i)
  neigT = np.r_[neighb1,neighb1[0]]
  # determine alpha and beta angles
  v1=(DistNeigh.T/(np.ones((3,1))*np.sqrt(np.sum(np.array([DistNeigh*DistNeigh]),2)))).T
  v2=NC[neigT[1:neigT.size],]-NC[neigT[0:neigT.size-1],]
  v2=(v2.T/(np.ones((3,1))*np.sqrt(np.sum(np.array([v2*v2]),2)))).T
  Aij = np.arccos(np.sum(v1*v2,1))
  Bij = np.arccos(np.sum(v1*-np.r_[v2[neighb1.size-1,].reshape((1,3)),v2[range(0,neighb1.size-1)]],1))
  Amix2 = 5
  AijBij = np.power(np.tan(np.r_[Aij[neighb1.size-1,],Aij[range(0,neighb1.size-1)]]),-1)+np.power(np.tan(np.r_[Bij[range(1,neighb1.size),],Bij[0,]]),-1)
  Kxi = np.sum(AijBij*DistNeigh,0)/Amix2
  return Aij,
    