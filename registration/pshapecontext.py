# -*- coding: utf-8 -*-

import numpy as np
import scipy as sp
from scipy.spatial import KDTree
from scipy.interpolate import Rbf	#Radial basis function module
from scipy.stats.mstats import find_repeats
from scipy.spatial.distance import cdist
from ptrimeshalt import *
import pprocess
import time
LIM = 8

def FeatPoints(Kmax,Kmin,NC,Tri,radius,alpha=0.5,beta=0.5):
  # Get feature points from computed maximum and minnimum feature curvature
  print 'Get Triangle and Vertex normals'
  VNORM = vertexnormal(NC,Tri)
  print 'Use Shape Index for Feature Point extraction'
  SI = 0.5 - np.arctan((Kmax+Kmin)/(Kmax-Kmin))/np.pi
  KDTnc = KDTree(NC,5)
  N1 = NC.shape[0]
  NPP1 = N1/LIM
  FeatPoints = np.array([np.nan]*NC.shape[0])
  results = pprocess.Map(limit=LIM)
  calc = results.manage(pprocess.MakeParallel(IsFeat))
  for j in range(0,LIM):
      calc(np.array(range(0,NPP1))+j*NPP1,NC,VNORM,KDTnc,radius,SI,alpha,beta,FeatPoints)
  for j in range(0,LIM):
      FeatPoints[np.array(range(0,NPP1))+j*NPP1] = results[j]
  FeatPoints[np.array(range(LIM*NPP1,N1))] = IsFeat(np.array(range(LIM*NPP1,N1)),NC,VNORM,KDTnc,radius,SI,alpha,beta,FeatPoints)
  FeatPoints = FeatPoints[FeatPoints>=0]
  return np.array(FeatPoints,int)
  

def IsFeat(nodes,NC,VNORM,KDTnc,radius,SI,alpha,beta,FeatPoints):
  for i in nodes:
    neighb = np.array(KDTnc.query_ball_point(NC[i,],radius),int)
    #if SI[i]<0.5:
      #neighb = neighb[SI[neighb,0]<0.5]
    #else:
      #neighb = neighb[SI[neighb,0]>0.5]
    neighb = neighb[np.dot(VNORM[neighb,],VNORM[i,].reshape((3,)))>0]
    mu = np.sum(SI[neighb])/len(neighb)
    cond1 = (SI[i]==np.max(SI[neighb]))&(SI[i]>=(1+alpha)*mu)
    cond2 = (SI[i]==np.min(SI[neighb]))&(SI[i]<=(1-beta)*mu)
    if cond1|cond2:
      FeatPoints[i] = i
  return FeatPoints[nodes]

def ShapeHist(NC,Tri,FeatPoints,radius,thetaB=12,phiB=12,rhoB=6):
  # Construct Shape Context histogram for Feature points
  N1 = FeatPoints.size
  NPP1 = N1/LIM
  print 'Get Triangle and Vertex normals'
  VNORM = vertexnormal(NC,Tri)
  KDTnc = KDTree(NC,5)
  print 'Set up Polar Histogram for given points'
  PolarHist = np.zeros((N1,12,12,6))
  results = pprocess.Map(limit=LIM)
  calc = results.manage(pprocess.MakeParallel(HistP))
  for j in range(0,LIM):
      calc(np.array(range(0,NPP1))+j*NPP1,FeatPoints,NC,VNORM,KDTnc,radius,wd,thetaB,phiB,rhoB)
  for j in range(0,LIM):
      PolarHist[np.array(range(0,NPP1))+j*NPP1,] = results[j]
  PolarHist[np.array(range(LIM*NPP1,N1)),] = HistP(np.array(range(LIM*NPP1,N1)),FeatPoints,NC,VNORM,KDTnc,radius,wd,thetaB,phiB,rhoB)
  return PolarHist
  
def HistP(nodes,FeatPoints,NC,VNORM,KDTnc,radius,dd):
  PHW = np.array([]*6)
  for inc in nodes:
    i = FeatPoints[inc]
    wd = VNORM[i].reshape((1,3))
    ud = np.cross(wd,np.cross(dd,wd))
    vd = np.cross(wd,ud)
    neighb =  np.array(KDTnc.query_ball_point(NC[i,],radius),int)
    NNC = NC[neighb[neighb<>i,],]-NC[i,]
    NNC_u,NNC_v,NNC_w = np.dot(NNC,ud.reshape((3,))),np.dot(NNC,vd.reshape((3,))),np.dot(NNC,wd.reshape((3,)))
    rho = cdist(NNC,np.array([[0,0,0]]),'euclidean')[:,0]
    theta,phi = np.arctan(NNC_v/NNC_u),np.arccos(NNC_w/rho)
    Logrho = np.log(rho)
    PolarHist[inc,] = sp.histogramdd(np.c_[theta,phi,Logrho],bins = (thetaB,phiB,rhoB))[0]*PHW
  return PolarHist[nodes,]
  
  
  