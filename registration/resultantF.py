# -*- coding: utf-8 -*-
import numpy as np                                                                                   
from enthought.mayavi.mlab import *
from fileread import *
from pmeshdeform import *
import ptetmeshdeform as ptet
from scipy.io import *
from scipy.interpolate import Rbf
from scipy.spatial import KDTree
from scipy.stats.mstats import find_repeats
from ptrimeshalt import *
from ptrifeatures import *
from plinedeform import *
import pshapecontext as sh
import pmeshquality as qu
import writeFEBio as wf
import time


alpha = 0 # alpha*Progn+(1-alpha)*Orthogn
fnameNodes = 'SkullPO_NC000_10'
fileName = 'SkullPO_000_10'
doWhat = 'M' #or 'I' ?
Tet = np.ma.load('SkullPO_TetUse000_10')
NC = np.ma.load(fnameNodes)
ndl = np.ma.load('SkullAverOuterBC')




PDIR = np.array([[-0.09,0.22,-0.97],[-0.36,-0.12,-0.91],[0.46,0.39,-0.78]])
ODIR = np.array([[-0.06,0.18,-0.97],[-0.15,-0.18,-0.96],[0.65,0.12,-0.75]])
MusclDir = alpha*PDIR + (1-alpha)*ODIR
MusclDir = MusclDir*(np.ones((3,1))*np.sum(MusclDir*MusclDir,1)).T

EQnc = np.zeros((25,3))
nrOfNds = np.zeros((25,))
for i in range(25):
  nodeU = ndl[ndl[:,1]==i,0]
  EQnc[i,] = np.sum(NC[nodeU,],0)/nodeU.size 
  nrOfNds[i] = nodeU.size 

NC = NC-EQnc[0,:]
EQnc = EQnc-EQnc[0,:]

MandibleNC = np.array([np.sum(EQnc[[5,6,7],],0)/3,np.sum(EQnc[[15,16,17],],0)/3])

TACT = np.array([0.86, 0.7, 0.56, 0.52, 0.44, 0.47, 0.49])
TempNweight = np.r_[TACT,TACT]
TempNweight = TempNweight*nrOfNds[[8,9,10,11,12,13,14,18,19,20,21,22,23,24]]
TempNweight[0:7] = TempNweight[0:7]/np.sum(TempNweight[0:7])
TempNweight[7:14] = TempNweight[7:14]/np.sum(TempNweight[7:14])

TempForceL = (-EQnc[[8,9,10,11,12,13,14],]+MandibleNC[0,])*(np.ones((3,1))*TempNweight[0:7]).T
TempForceR = (-EQnc[[18,19,20,21,22,23,24],]+MandibleNC[1,])*(np.ones((3,1))*TempNweight[7:14]).T

TempForceL = 410*TempForceL/np.sqrt(np.sum(np.sum(TempForceL,0)*np.sum(TempForceL,0)))
TempForceR = 170*TempForceR/np.sqrt(np.sum(np.sum(TempForceR,0)*np.sum(TempForceR,0)))

LefT = np.array([230,100,250])
RighT = np.array([108,36,50])

MusclFL = np.r_[(np.ones((3,1))*LefT).T*MusclDir,TempForceL]
MusclFR = np.r_[(np.ones((3,1))*RighT).T*np.c_[-MusclDir[:,0],MusclDir[:,1:]],TempForceR]
MD = EQnc[5:,]
Force = np.r_[MusclFL,MusclFR]

if doWhat == 'M':
  RD = EQnc[[1,3,4],]
else:
  RD = EQnc[[2,3,4],]
  
MATR = np.zeros((6,6))
MATR[0,0],MATR[1,1],MATR[1,2],MATR[2,3],MATR[2,4],MATR[2,5] = 1,1,1,1,1,1
MATR[3,1],MATR[3,2],MATR[3,3],MATR[3,4],MATR[3,5]=-RD[1,2],-RD[2,2],RD[1,1],RD[2,1],RD[0,1]
MATR[4,0],MATR[4,3],MATR[4,4],MATR[4,5]=(RD[1,2]+RD[2,2])/2,-RD[1,0],-RD[2,0],-RD[0,0]
MATR[5,0],MATR[5,1],MATR[5,2]=-(RD[1,1]+RD[2,1])/2,RD[1,0],RD[2,0]
MatINV = np.linalg.pinv(MATR)

RESUL = np.zeros((6,1))
RESUL[0],RESUL[1],RESUL[2]=-np.sum(Force[:,0]),-np.sum(Force[:,1]),-np.sum(Force[:,2])
RESUL[3] = np.sum(Force[:,1]*MD[:,2])-np.sum(Force[:,2]*MD[:,1])
RESUL[4] = np.sum(Force[:,2]*MD[:,0])-np.sum(Force[:,0]*MD[:,2])
RESUL[5] = np.sum(Force[:,0]*MD[:,1])-np.sum(Force[:,1]*MD[:,0])

Freact = np.linalg.solve(MATR,RESUL)
FAE = np.zeros((3,3))
FAE[0,2] = Freact[5]
FAE[1,],FAE[2,] = np.array([Freact[0,0]/2, Freact[1,0], Freact[3,0] ]),np.array([Freact[0,0]/2, Freact[2,0], Freact[4,0] ])

Force = np.r_[FAE,Force]
nodeConst = ndl[ndl[:,1]==0,0]
nrOfLCs = 67
nodeLists = [[]]*67
if doWhat == 'M':
  nodeLists[0] = ndl[ndl[:,1]==1,0]
else:
  nodeLists[0] = ndl[ndl[:,1]==2,0]
for i in range(22):
  nodeLists[3*(i+1)-2],nodeLists[3*(i+1)-1],nodeLists[3*(i+1)] = ndl[ndl[:,1]==i+3,0],ndl[ndl[:,1]==i+3,0],ndl[ndl[:,1]==i+3,0]

Fdir = ['z']+['x','y','z']*22
scaleValue = np.zeros((67,))
if doWhat == 'M':
  scaleValue[0] = 1./nrOfNds[1]
else:
  scaleValue[0] = 1./nrOfNds[2]
for i in range(22):
  scaleValue[3*(i+1)-2],scaleValue[3*(i+1)-1],scaleValue[3*(i+1)] = 1./nrOfNds[i+3],1./nrOfNds[i+3],1./nrOfNds[i+3]
ForceVal = np.zeros((67,))
ForceVal[0] = Force[0,2]
for i in range(22):
  ForceVal[3*(i+1)-2],ForceVal[3*(i+1)-1],ForceVal[3*(i+1)] = Force[i+1,0],Force[i+1,1],Force[i+1,2]
LoadC = [[]]*67
for i in range(67):
  LoadC[i] = [Fdir[i], i+1, scaleValue[i], ForceVal[i]]
  
wf.NodeForce(NC,Tet,nodeConst,nodeLists,LoadC,fileName)