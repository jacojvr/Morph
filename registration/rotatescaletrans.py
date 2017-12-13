# -*- coding: utf-8 -*-
import numpy as np
from fileread import *
import scipy as sp
from scipy.optimize import *

NCB,LMB,LMBI = [],[],[]

def SymmNC(NC,LM1,LM2):
  global NCB,LMB,LMBI
  NCB = NC - np.ones((NC.size/3,1))*np.sum(NC[LM1[:,0],],0)/LM1.size
  NCBI = np.c_[-NC[:,0],NC[:,1:3]]
  LMB=LM1
  LMBI=LM2
  xf = fmin_powell(costSymmetry,np.array([0,0,1,0,0,0,0]))
  directionvec = xf[0:3]
  theta = xf[3,]
  NCS = np.c_[NCB[:,0]+xf[4,],NCB[:,1]+xf[5,],NCB[:,2]+xf[6,]]
  x,y,z=directionvec[0:3,]/np.sqrt(np.sum(np.power(directionvec,2)))
  omcost = 1-np.cos(theta)
  sint = np.sin(theta)
  RMat = np.array([[1+omcost*(x*x-1), -z*sint+omcost*x*y, y*sint+omcost*x*z],
  [z*sint+omcost*x*y, 1+omcost*(y*y-1), -x*sint+omcost*y*z],
  [-y*sint+omcost*x*z, x*sint+omcost*y*z, 1+omcost*(z*z-1)]])
  NCS = np.dot(NCS,RMat)
  return NCS
  
def costSymmetry(XO):
  directionvec = XO[0:3]
  theta = XO[3,]
  global NCB,LMB,LMBI
  NCS = np.c_[NCB[:,0]+XO[4,],NCB[:,1]+XO[5,],NCB[:,2]+XO[6,]]
  x,y,z=directionvec[0:3,]/np.sqrt(np.sum(np.power(directionvec,2)))
  omcost = 1-np.cos(theta)
  sint = np.sin(theta)
  RMat = np.array([[1+omcost*(x*x-1), -z*sint+omcost*x*y, y*sint+omcost*x*z],
  [z*sint+omcost*x*y, 1+omcost*(y*y-1), -x*sint+omcost*y*z],
  [-y*sint+omcost*x*z, x*sint+omcost*y*z, 1+omcost*(z*z-1)]])
  NCS = np.dot(NCS,RMat)
  NCSI = np.c_[-NCS[:,0],NCS[:,1:3]]
  Points1 = NCS[LMB,]
  Points2 = NCSI[LMBI,]
  return np.sum(np.sqrt(sum(np.power(Points1-Points2,2),1)))
    