# -*- coding: utf-8 -*-

import numpy as np
import scipy as sp
import scipy.optimize as opt
from scipy.spatial import KDTree
from scipy.interpolate import Rbf	#Radial basis function module
from scipy.stats.mstats import find_repeats
from scipy.spatial.distance import cdist
from ptrimeshalt import *
import pysparse as ps
import pprocess
import time

LIM = 8

def OptINT(inner,outer,NC,Tet,neighbListTet,iterations,boundPenalty=10):
  N1 = NC.shape[0]
  NPP1 = N1/LIM
  NCcur = np.r_[NC]
  T1 = Tet.shape[0]
  TPP1 = T1/LIM
  print 'Submesh Optimization on ',inner.size,' of ',NC.shape[0],' nodes'
  for i in range(iterations):
    print '	Iteration ',i+1
    NCprev = np.r_[NCcur]
    EQ,delt,Sn2,Sig = elemQual_mu(np.array(range(Tet.shape[0])),NCprev,Tet)
    np.ma.dump(EQ,'ElQual'+str(i+1))
    #EQ = np.ma.load('ElQual1')
    print '			Max ',np.max(EQ),'; Min ',np.min(EQ),'; Average ',np.average(EQ)
    print '			number of degenerate elements : ',np.where(Sig<0)[0].size
    print '		build optimization matrix G'
    global Gmat
    Gmat = ps.spmatrix.ll_mat(NC.shape[0]+outer.size,NC.shape[0])
    Gmat.put(1,np.array(range(NC.shape[0])),np.array(range(NC.shape[0])))
    rows = np.array(range(outer.size))+NC.shape[0]
    Gmat.put(1,rows,outer)
    print '			number of initial non-zeros: ',Gmat.nnz
    EQinv = 1/EQ
    GMparts = [[]]*LIM
    results = pprocess.Map(limit=LIM)
    calc = results.manage(pprocess.MakeParallel(US_Gmat))
    for j in range(LIM):
      calc(np.array(range(0,NPP1))+j*NPP1,NCprev,Tet,neighbListTet,EQinv,Sig,j)
    for j in range(LIM):
      vals,rows,cols = results[j]
    #for j in range(LIM):
      Gmat.put(vals,rows,cols)
    if N1>LIM*NPP1:
      vals,rows,cols = US_Gmat(np.array(range(LIM*NPP1,N1)),NCprev,Tet,neighbListTet,EQinv,Sig,1)
      Gmat.put(vals,rows,cols)
    print '			final number of non-zeros: ',Gmat.nnz
    Gmat.export_mtx('Gmat'+str(i+1))
    print '		do optimization using conjugate gradient'
    GTG = ps.spmatrix.dot(Gmat,Gmat)
    gVec = ps.spmatrix.ll_mat(NC.shape[0]+outer.size,3)
    for j in outer:
      rows = np.array([j,j,j])+NC.shape[0]
      cols = np.array([0,1,2])
      vals = NC[j,]
      gVec.put(vals,rows,cols)
    GTg = ps.spmatrix.dot(Gmat,gVec)
    GTgNUMPY = np.zeros(NC.shape)
    GTgNUMPY[GTg.keys()] = GTg.values()
    mu = GTg.norm('inf')
    mu=mu*mu
    if mu<1:
      mu=1
    GTG.update_add_at(mu*np.ones((NC.shape[0],)),np.array(range(NC.shape[0])),np.array(range(NC.shape[0])))
    RHS = GTgNUMPY+mu*NCprev
    X,Y,Z = np.zeros((NC.shape[0],)),np.zeros((NC.shape[0],)),np.zeros((NC.shape[0],))
    print 'X'
    infx,itex,resx = ps.itsolvers.cgs(GTG,RHS[:,0],X,0.00000000001,100)
    print 'Y'
    infy,itey,resy = ps.itsolvers.cgs(GTG,RHS[:,1],Y,0.00000000001,100)
    print 'Z'
    infz,itez,resz = ps.itsolvers.cgs(GTG,RHS[:,2],Z,0.00000000001,100)
    NCcur = np.c_[X,Y,Z]
    np.ma.dump(NCcur,'NCopt_temp'+str(i+1))
  return NCcur

def OptEQ(NC,Tet,move,ITER=10,degen=0.15):
  for i in range(ITER):
    print 'ITERATION ',i+1
    NCprev = np.r_[NC]
    Jac = DetTJ(np.array(range(Tet.shape[0])),NC,Tet)
    np.ma.dump(Jac,'Jac'+str(i+1))
    TETS = np.where(Jac<degen)[0]
    print ' Number of degenerate triangles : ',TETS.size
    EQ,deltPrescr,Sn2,Sig = elemQual_mu(np.array(range(Tet.shape[0])),NC,Tet)
    np.ma.dump(EQ,'ElQual'+str(i+1))
    print ' Element quality : '
    print '		Max ',np.max(EQ),'; Min ',np.min(EQ),'; Average ',np.average(EQ)
    if TETS.size>0:
      #nodes = Tet[TETS,:].reshape((4*TETS.size,))
      #move = np.array(find_repeats(nodes)[0],int)
      #move = np.array(range(NC.shape[0]))
      #for i in range(NC.shape[0]):
	#if np.where(nodes==i)[0].size==0:
	  #move = move[move<>i]
      print ' Number of nodes to move : ',move.size
      NCopt = opt.fmin_cg(EQweight,np.r_[NC[move,]],args=(NCprev,move,Tet,EQ,deltPrescr))[0]
      NC[move,] = NCopt.reshape((move.size,3))
      np.ma.dump(NC,'NCopt_'+str(i+1))
    else:
      return NC
  return NC

def EQweight(NCinner,NC,inner,Tet,EQprev,deltPrescr):
  EQ=np.r_[EQprev]
  NCcur = np.r_[NC]
  NCDIFF = NCinner.reshape((inner.size,3))-NC[inner,]
  NCcur[inner,] = NCinner.reshape((inner.size,3))
  nodes = inner[np.where(np.sum(NCDIFF*NCDIFF,1)>0)[0]]
  tets=np.array([])
  if nodes.size>0:
    for i in nodes:
      tets = np.r_[tets,np.where(Tet==i)[0]]
    #tets = np.array(tets,int)
    tets = np.array(find_repeats(np.r_[tets,tets])[0],int)
  if tets.size>0:
    EQ[tets] = elemQual_mu(tets,NCcur,Tet,0,deltPrescr)[0]
  COST = np.sum(1/EQ)
  return COST 
    
def tetVec(nodes,inner,NC,neighbListTet,QList):
  VecPart = np.zeros((nodes.size,3))
  for ii in range(nodes.size):
    inc = nodes[ii]
    i = inner[inc]
    VecPart[ii,] = NC[i,] - np.sum(NC[neighbListTet[i],]*np.c_[QList[inc],QList[inc],QList[inc]],0)
  return VecPart
    
def DetTJ(elem,NC,Tet):
  Jacobian = np.zeros((elem.size,))
  for i in range(elem.size):
    el = elem[i]
    Vect = NC[Tet[el,1:],]-NC[Tet[el,0],]
    Jacobian[i] = np.dot(np.cross(Vect[0,],Vect[1,]),Vect[2,])/8
  return Jacobian
  
def elemQual_mu(elem,NC,Tet,disp=1,deltPresc = 0):
  if disp==1:
    print 'Determine element Quality [mu]'
  gamma = 0.0000001
  delta = 0
  T1 = elem.shape[0]
  TPP1 = T1/LIM
  Sn2 = np.zeros((elem.shape[0],))
  Sig = np.zeros((elem.shape[0],))
  if elem.size>LIM:
    results = pprocess.Map(limit=LIM)
    calc = results.manage(pprocess.MakeParallel(Sn2Sig))
    for j in range(0,LIM):
	calc(elem[np.array(range(0,TPP1))]+j*TPP1,NC,Tet)
    for j in range(0,LIM):
	Sn2[j*TPP1:(1+j)*TPP1],Sig[j*TPP1:(1+j)*TPP1] = results[j]
  if np.array(range(LIM*TPP1,T1)).size>0:
    Sn2[LIM*TPP1:T1],Sig[LIM*TPP1:T1]=Sn2Sig(elem[np.array(range(LIM*TPP1,T1))],NC,Tet)
  if np.min(Sig)<gamma:
    if disp==1:
      print '	minnimum Sig = ',np.min(Sig)
    delta = np.sqrt(gamma*(gamma-np.min(Sig)))
  if delta<deltPresc:
    delta=deltPresc
  if disp==1:
    print '	delta = ',delta
  h_sig = (Sig + np.sqrt(Sig*Sig + 4*delta*delta))/2
  return 3*np.power(h_sig,2./3)/Sn2,delta,Sn2,Sig 
    
def Sn2Sig(elem,NC,Tet):
  Sn2 = np.zeros((elem.size,))
  sig = np.zeros((elem.size,))
  W = np.linalg.inv(np.matrix([[1,0.5,0.5],[0,np.sqrt(3)/2,np.sqrt(3)/6],[0,0,np.sqrt(2./3)]]))
  for i in range(elem.size):
    el = elem[i]
    Sel = (NC[Tet[el,1:],]-NC[Tet[el,0],]).T*W
    Sn2[i] = np.trace(Sel.T*Sel)
    sig[i] = np.linalg.det(Sel)
  return Sn2,sig
  
def US_Gmat(nSubM,NC,Tet,neighbListTet,elemQualInv,Sig,showPercent=0):
  Gmatprt = ps.spmatrix.ll_mat(NC.shape[0],NC.shape[0])
  strt = 0
  W = np.linalg.inv(np.matrix([[1,0.5,0.5],[0,np.sqrt(3)/2,np.sqrt(3)/6],[0,0,np.sqrt(2./3)]]))
  for inc in range(nSubM.shape[0]):
  ###for i in nSubM:
    if (showPercent == 0):
      if (inc*100/nSubM.shape[0]>strt):
	strt = inc*100/nSubM.shape[0]
	if np.mod(strt,5)==0:
	  print strt,' %'
    i = nSubM[inc]
    cols = neighbListTet[i]
    rows = np.array(i*np.ones(cols.shape),int)
    tets = np.where(Tet==i)[0]
    TetSub_i = Tet[tets,]
    TriSub_i = TetSub_i[TetSub_i<>i]
    TriSub_i = TriSub_i.reshape((TriSub_i.size/3,3))
    Weights = elemQualInv[tets]#np.zeros((TriSub_i.shape[0],))
    Row_i = np.zeros((cols.shape))#NC.shape[0],))
    Centr = NC[i,]
    neighbNC = NC[cols,]
    for j in range(TriSub_i.shape[0]):
      tri = TriSub_i[j,]
      tri_NC = NC[tri,]
      barryC = np.sum(tri_NC,0)/3
      normv = np.cross(tri_NC[1,]-tri_NC[0,],tri_NC[2,]-tri_NC[0,])
      if ((np.dot(normv,Centr-barryC)<0)&(Sig[tets[j]]>0))|((np.dot(normv,Centr-barryC)>0)&(Sig[tets[j]]<0)):
	tri = TriSub_i[j,[0,2,1]]
	normv = -normv
	tri_NC = NC[tri,]
      NCnewj = np.sum(tri_NC,0)/3+opt.fmin(optPos,np.array([0.5]),args=(tri_NC,normv,W),disp=0)[0]*normv
      Distj = neighbNC-NCnewj
      Distj = np.sum(Distj*Distj,1)
      Row_i = Row_i+(Weights[j]*Weights[j])/(Distj*Distj)
    #Row_i = Row_i[cols]
    values = -Row_i/np.sum(Row_i)
    Gmatprt.put(values,rows,cols)
  return Gmatprt.values(),Gmatprt.keys()[0],Gmatprt.keys()[1]
      
def optPos(alph,tri_NC,norm,W):
  pt = np.sum(tri_NC,0)/3 + alph*norm
  Sel = (tri_NC[[2,1,0],]-pt).T*W
  Sn2 = np.trace(Sel.T*Sel)
  sig = np.linalg.det(Sel)
  return Sn2/(3*np.power(sig,2./3))