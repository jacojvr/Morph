# -*- coding: utf-8 -*-
'''
THIS PROGRAM/ALGORITHM IS NOT YET COMPLETE/IN ANY WAY FUNCTIONAL
'''
import numpy as np
import scipy as sp
from scipy.spatial import KDTree
from scipy.interpolate import Rbf	#Radial basis function module
import time

def elastictri(NCB,ConnectB,LMB,NCT,ConnectT,LMT,Sampling,k_max):
  # Elastic surface registration where triangles are moved to a target surface:
  #inputs:
  # NCB,NCT: nodal coordinates of base and target surfaces
  # ConnectB;ConnectT: Base&target connectivity
  # LMB,LMT: Base and target landmarks that have to correspond
  # Sampling: sampling rate (number of nodes to use...equally spaced at the moment)
  # 	if sampling is chosen as "0", all nodes are taken into account
  # k_max: maximum number of iterations...
  ######## THE OTHER SETTINGS ARE SET INTERNAL TO CODE #########
  t_start = time.clock()
  ConnectB = np.array(ConnectB-1,int)
  ConnectT = np.array(ConnectT-1,int)
  # Settings for elastic surface registration:
  m=10 # nearest neighbour parameter
  alph=0.5 # normilization factor
  gamm=2 # smoothing parameter1
  sigm0=10 # smoothing parameter2
  f=1.0715 # smoothing parameter3
  #k_max=2 # max iterations
  Tol=0.0001 # stopping criteria

  # determine N1,N2,T1 and T2:
  N1=np.array(np.shape(NCB))[0,]
  N2=np.array(np.shape(NCT))[0,]
  T1=np.array(np.shape(ConnectB))[0,]
  T2=np.array(np.shape(ConnectT))[0,]

  #Compute target-mesh triangle centroids:
  print "determining centroids and normal directions of target surface triangles "
  S_2_centr = np.c_[np.sum(np.c_[NCT[ConnectT[:,0],0],NCT[ConnectT[:,1],0],NCT[ConnectT[:,2],0]],1)/3,
    np.sum(np.c_[NCT[ConnectT[:,0],1],NCT[ConnectT[:,1],1],NCT[ConnectT[:,2],1]],1)/3,np.sum(np.c_[NCT[ConnectT[:,0],2],NCT[ConnectT[:,1],2],NCT[ConnectT[:,2],2]],1)/3]
  S_2_normvec = np.cross(NCT[ConnectT[:,1],:]-NCT[ConnectT[:,0],:],NCT[ConnectT[:,2],:]-NCT[ConnectT[:,0],:])
  S_2_normvec = (S_2_normvec.T/(np.ones((3,1))*np.sqrt(np.sum(np.array([S_2_normvec*S_2_normvec]),2)))).T
  print "determining kd-trees of target surface centroids and nodal coordinates"
  KDT_TC = KDTree(S_2_centr,m)
  KDT_TN = KDTree(NCT,m)

  print "initialize iterations, deforming mesh surface and displacement matrices"
  k = 1
  W_km1 = NCB
  if Sampling ==1:
    SamplingB=T1
    SamplingT=T2
  else:
    SamplingB=Sampling
    SamplingT=Sampling

  LMB = np.array(np.ceil(np.linspace(1,T1,SamplingB))-1,int)
  LMT = np.array(np.ceil(np.linspace(1,T2,SamplingT))-1,int)
  D1 = np.zeros((SamplingB,3))
  D2 = np.zeros((SamplingT,3))
  DS1 = np.zeros((N1,3))
  DS2 = np.zeros((N1,3))

  # MAIN MESH DEFORMATION ALGORITHM:
  while k<=k_max:
    print
    print "MESH DEFORMATION ITERATION",k
    print "	determining centroids and normal directions of deforming mesh"
    W_km1_centr = np.c_[np.sum(np.c_[W_km1[ConnectB[:,0],0],W_km1[ConnectB[:,1],0],W_km1[ConnectB[:,2],0]],1)/3,
      np.sum(np.c_[W_km1[ConnectB[:,0],1],W_km1[ConnectB[:,1],1],W_km1[ConnectB[:,2],1]],1)/3,np.sum(np.c_[W_km1[ConnectB[:,0],2],W_km1[ConnectB[:,1],2],W_km1[ConnectB[:,2],2]],1)/3]
    W_km1_normvec = np.cross(W_km1[ConnectB[:,1],:]-W_km1[ConnectB[:,0],:],W_km1[ConnectB[:,2],:]-W_km1[ConnectB[:,0],:])
    W_km1_normvec = (S_2_normvec.T/(np.ones((3,1))*np.sqrt(np.sum(np.array([S_2_normvec*S_2_normvec]),2)))).T
    print "	determining kd-tree of current deforming surface centroids and nodal coordinates"
    KDT_KC = KDTree(W_km1_centr,m)
    KDT_KN = KDTree(W_km1,m)

    print "	setting up D1(i,d)"
    tic = time.clock()
    for i in range(0,SamplingB):
      nn=LMB[i]
      D1[i,:] = Ndisp(nn,W_km1,NCT,ConnectT,S_2_centr,KDT_TC,KDT_TN)
    print time.clock()-tic
    
    print "	setting up D2(j,c)"
    tic = time.clock()
    for j in range(0,SamplingT): 
      nn=LMT[j]
      D2[j,:] = Ndisp(nn,NCT,W_km1,ConnectB,W_km1_centr,KDT_KC,KDT_KN)
    
    # Compute displacement update for each node using suggested Gaussian radial basis function:
    print time.clock()-tic
    print "	determining smoothed displacement field"
    tic=time.clock()
    

    # Deform mesh using scipy.interpolate.Rbf(*args):
    #nKD=np.r_[W_km1[LMB,],NCT[LMT,]+D2/1.2]
    #Disp=np.r_[D1,-D2]/10
    #rbfx = Rbf(nKD[:,0],nKD[:,1],nKD[:,2],Disp[:,0])
    #rbfy = Rbf(nKD[:,0],nKD[:,1],nKD[:,2],Disp[:,1]) 
    #rbfz = Rbf(nKD[:,0],nKD[:,1],nKD[:,2],Disp[:,2]) 
    #DSx = rbfx(W_km1[:,0],W_km1[:,1],W_km1[:,2])
    #DSy = rbfy(W_km1[:,0],W_km1[:,1],W_km1[:,2])
    #DSz = rbfz(W_km1[:,0],W_km1[:,1],W_km1[:,2])
    #DS = np.c_[DSx,DSy,DSz]
    

    # Deform mesh using known displacements and thin plate spline RBF:
    #DS1 = RBFmorph(W_km1,W_km1[LMB[:,],:],D1*0.2)
    #DS2 = RBFmorph(W_km1,(NCT[LMT[:,],:]+D2*0.9),-D2*0.2)
    #DS = DS1+DS2

    sigma_k2 = np.power(np.power(f,-k)*sigm0,2)
    for nodes in range(0,N1):
      G1node = np.array([W_km1[nodes,:]]).T*np.ones((1,SamplingB))-W_km1[LMB,:].T
      G1node = np.exp(-np.sum(G1node*G1node,0)/sigma_k2)
      G2node = np.array([W_km1[nodes,:]]).T*np.ones((1,SamplingT))-NCT[LMT,:].T-D2.T
      G2node = np.exp(-np.sum(G2node*G2node,0)/sigma_k2)
      DS[nodes,:] = (np.sum(D1.T*G1node,1)/np.sum(G1node,0)-np.sum(D2.T*G2node,1)/np.sum(G2node,0))/gamm
    
    # Deform mesh:
    print time.clock()-tic
    print "	deforming mesh (update of W_{k-1})"
    W_km1 = W_km1+DS

    # determine wheter convergence is acheived
    TotalMorph = np.sum(np.sqrt(np.sum(DS*DS,1)),0)
    print "	total displacement for current deformation iteration:"
    print TotalMorph
    if TotalMorph<Tol:
      print
      print "SOLUTION CONVERGED to within specified tolerance of:",Tol
      k = k_max+5
    else:
      print "problem not yet converged at iteration",k
      k = k+1;
      
    if k==k_max+1:
      print
      print "SOLUTION TERMINATED: maximum iterations,(",k_max,") reached"
    print
    print "TOTAL TIME FOR ELASTIC SURFACE REGISTRATION : ",time.clock()-t_start,"seconds"
    return W_km1
