# -*- coding: utf-8 -*-

import numpy as np
import scipy as sp
from scipy.spatial import KDTree
from scipy.interpolate import Rbf	#Radial basis function module
from scipy.stats.mstats import find_repeats
from trimeshalt import *
import time

def elasticsurf(NCB,ConnectB,LandmB,NCT,ConnectT,LandmT,FDNB,SFPNT,SFPTT,k_max):
  # Elastic surface registration:
  #inputs:
  # NCB,NCT: nodal coordinates of base and target surfaces
  # ConnectB;ConnectT: Base&target connectivity
  # LandmB,LandmT landmarks that have to have a 1 to 1 correspondence (input 0 no landmarks are present)
  # FDN: Feature dependant nodes on Base-mesh (indices in NCB)
  # SFPNT&SFPTT: Selective Feature preserving nodes and triangles (indices in NCT and ConnectT) on target mesh.
  # k_max: maximum number of iterations
  #	1) The base mesh is first morphed into the target using 1 to 1 correspondence with target mesh Landmarks.
  #	2) Elastic surface registration is done using FDNB&SFPNT to target and base mesh respectively
  #		if FDNB is displaced to a triangle other than SFPTT, displacement of concern is discarded in order
  #		to retain only selected features.
  ######## ADDITIONAL SETTINGS REQUIRED ARE SET INTERNAL TO CODE #########
  print
  print "SELECTIVE MESH MORPHING ALGORITHM USING ELASTIC SURFACE REGISTRATION"
  print "	-G.J.J.v.Rensburg - 22/04/2010-"
  t_start = time.clock()
  ConnectB = np.array(ConnectB,int)
  ConnectT = np.array(ConnectT,int)
  LandmB = np.array(LandmB[:,0],int)		# do -1 later to be consistent with python indexing, first need to do other "temporary landmarks"& check that they dont fall on actual landmark positions!
  LandmT = np.array(LandmT[:,0],int)

  # Settings for elastic surface registration:
  m=10 # nearest neighbour parameter
  alph=0.5 # normilization factor
  gamm=2 # smoothing parameter1
  sigm0=10 # smoothing parameter2
  f=1.0715 # smoothing parameter3
  Tol=0.0001 # stopping criteria

  # determine N1,N2,T1 and T2:
  N1=np.array(np.shape(NCB))[0,]
  N2=np.array(np.shape(NCT))[0,]
  T1=np.array(np.shape(ConnectB))[0,]
  T2=np.array(np.shape(ConnectT))[0,]
  NL=np.array(np.shape(LandmB))[0,]

  ################################     INITIALIZE & NODES OF CONCERN:    #################################
  ########################################################################################################
  print
  print "INITIALIZE SURFACE DEFORMATION"
  k = 1
  CONV = np.zeros((k_max,1))
  print " 	enquire nodes where required displacement is checked"
  # remove Landmarks from FDNB and SFPNT:
  for i in range(0,NL):
    if find_repeats(np.r_[FDNB,LandmB[i,]])[0].size>0:
      r=np.where(FDNB==LandmB[i,])[0]
      FDNB = np.r_[FDNB[0:r,],FDNB[r+1:FDNB.size,]]
    if find_repeats(np.r_[SFPNT,LandmT[i,]])[0].size>0:
      r=np.where(SFPNT==LandmT[i,])[0]
      SFPNT = np.r_[SFPNT[0:r,],SFPNT[r+1:SFPNT.size,]]
  SamplingB=FDNB.size
  SamplingT=SFPNT.size
  LMB = np.r_[FDNB,LandmB]	# Last NL entries are reserved for Landmarks that HAVE TO FIT points on the target mesh
  LMT = np.r_[SFPNT,LandmT]
  print "	Compute known displacement for Base_Landmark to Target_Landmark"
  knownC = NCB[LandmB,]
  knownD = NCT[LandmT,]-knownC
  print	
  print "COARSE SURFACE REGISTRATION"
  print "	using landmark displacements to deform using RBF"
  W_km1 = RBFmorph(NCB,NCB[LandmB,],NCT[LandmT,]-NCB[LandmB])
  
  ################################    MAIN MESH DEFORMATION ALGORITHM:   #################################
  ########################################################################################################


  print
  print "ELASTIC SURFACE REGISTRATION"  
  #Compute target-mesh triangle centroids:
  print "determining centroids of target surface triangles"
  S_2_centr = np.c_[np.sum(np.c_[NCT[ConnectT[:,0],0],NCT[ConnectT[:,1],0],NCT[ConnectT[:,2],0]],1)/3,
    np.sum(np.c_[NCT[ConnectT[:,0],1],NCT[ConnectT[:,1],1],NCT[ConnectT[:,2],1]],1)/3,np.sum(np.c_[NCT[ConnectT[:,0],2],NCT[ConnectT[:,1],2],NCT[ConnectT[:,2],2]],1)/3]
  print "determining kd-trees of target surface centroids and nodal coordinates"
  KDT_TC = KDTree(S_2_centr,m)
  KDT_TN = KDTree(NCT,m)
  while k<=k_max:
    D1 = np.zeros((SamplingB,3))
    D2 = np.zeros((SamplingT,3))
    DS = np.zeros((N1,3))
    print
    print "MESH DEFORMATION ITERATION",k
    print "determining known displacement of landmarks"
    knownD = NCT[LandmT,]-W_km1[LandmB,]
    print "	determining centroids of deforming mesh"
    W_km1_centr = np.c_[np.sum(np.c_[W_km1[ConnectB[:,0],0],W_km1[ConnectB[:,1],0],W_km1[ConnectB[:,2],0]],1)/3,
      np.sum(np.c_[W_km1[ConnectB[:,0],1],W_km1[ConnectB[:,1],1],W_km1[ConnectB[:,2],1]],1)/3,np.sum(np.c_[W_km1[ConnectB[:,0],2],W_km1[ConnectB[:,1],2],W_km1[ConnectB[:,2],2]],1)/3]
    print "	determining kd-tree of current deforming surface centroids and nodal coordinates"
    KDT_KC = KDTree(W_km1_centr,m)
    KDT_KN = KDTree(W_km1,m)

    print "	setting up D1(i,d)"
    tic = time.clock()
    for i in range(0,SamplingB):
      nn=LMB[i]
      # query kd-tree for closest triangle to node:
      ncl = KDT_TC.query(W_km1[nn,:])[1]
      # check if 
      if np.where(SFPTT==ncl)[0].size>0:
	## determine target triangle normal vector:
	#tnorm = np.cross(NCT[ConnectT[ncl,1],:]-NCT[ConnectT[ncl,0],:],NCT[ConnectT[ncl,2],:]-NCT[ConnectT[ncl,0],:])
	#tnorm = tnorm/np.sqrt(np.sum(tnorm*tnorm))
	## determine current vertex normal
	#vnorm=vertexnormal(W_km1,ConnectB[np.where(ConnectB==nn)[0],])
	#if np.dot(vnorm[nn,],tnorm)>0:	#check for correlation between closest triangle directional normal&base-mesh curvature
	  ## move to triangle / closest node
	D1[i,:] = Ndisp(W_km1[nn,:],ncl,NCT,ConnectT,S_2_centr,KDT_TN)
    D1=np.r_[D1,knownD]
    print "			",time.clock()-tic," seconds"
    
    print "	setting up D2(j,c)"
    tic = time.clock()
    for j in range(0,SamplingT):
      nn=LMT[j]
      ncl = KDT_KC.query(NCT[nn,:])[1]
      #tnorm = np.cross(W_km1[ConnectB[ncl,1],:]-W_km1[ConnectB[ncl,0],:],W_km1[ConnectB[ncl,2],:]-W_km1[ConnectB[ncl,0],:])
      #tnorm = tnorm/np.sqrt(np.sum(tnorm*tnorm))
      #vnorm=vertexnormal(NCT,ConnectT[np.where(ConnectT==nn)[0],])
      #if np.dot(vnorm[nn,],tnorm)>0:
      D2[j,:] = Ndisp(NCT[nn,:],ncl,W_km1,ConnectB,W_km1_centr,KDT_KN)
      #else:
	#D2[j,:] = np.array([0,0,0])
    D2 = np.r_[D2,-knownD]
    print "			",time.clock()-tic," seconds"
   
    # Compute displacement update for each node using suggested Gaussian radial basis function:
    print "	determining smoothed displacement field"
    tic=time.clock()
    # Deform mesh using Gaussian smoothing as suggested in paper by R.Bryan et al.
    sigma_k2 = np.power(np.power(f,-k)*sigm0,2)
    for nodes in range(0,N1):
      G1node = np.array([W_km1[nodes,:]]).T*np.ones((1,SamplingB+NL))-W_km1[LMB,:].T
      G1node = np.exp(-np.sum(G1node*G1node,0)/sigma_k2)
      G2node = np.array([W_km1[nodes,:]]).T*np.ones((1,SamplingT+NL))-NCT[LMT,:].T-D2.T
      G2node = np.exp(-np.sum(G2node*G2node,0)/sigma_k2)
      DS[nodes,:] = (np.sum(D1.T*G1node,1)/np.sum(G1node,0)-np.sum(D2.T*G2node,1)/np.sum(G2node,0))/gamm
      
    
    print "			",time.clock()-tic," seconds"
    # determine wheter convergence is acheived
    TotalMorph = np.sum(np.sqrt(np.sum(DS*DS,1)),0)
    print "	total displacement for current deformation iteration:"
    print TotalMorph
    if (TotalMorph<Tol):
      print
      print "CONVERGED SOLUTION OBTAINED"
      CONV[k-1,0] = TotalMorph
      k = k_max*10+1
      W_km1 = W_km1+DS
    elif (k<10)|(TotalMorph<=CONV[0,]):
      print "problem not yet converged at iteration",k
      CONV[k-1,0] = TotalMorph
      k = k+1;
      # Deform mesh:
      print "	deforming mesh (update of W_{k-1})"
      W_km1 = W_km1+DS
    else:
      print "PROBLEM DIVERGING"
      k=k_max*10-1

    if (k>1)&(np.mod(k-1,10)==0):
      print
      print "Do 3 iterations of Laplacian Smoothing to improve element quality"
      tic = time.clock()
      W_km1 = LaplacianSmooth(W_km1,ConnectB,3)
      print "		Smoothing done in ",time.clock()-tic," seconds"
      print "COARSE SURFACE REGISTRATION"
      print "	using landmark displacements to deform using RBF"
      W_km1 = RBFmorph(W_km1,W_km1[LandmB,],NCT[LandmT,]-W_km1[LandmB,])
    print
      
  if k==k_max+1:
    print
    print "SOLUTION TERMINATED: maximum iterations,(",k_max,") reached"
  print
  print "TOTAL TIME FOR ELASTIC SURFACE REGISTRATION : ",time.clock()-t_start,"seconds"
  return W_km1,CONV


def LaplacianSmooth(NCL,Tri,it):
  print "LAPLACIAN SMOOTHING:"
  for i in range(0,it):
    print "	iteration ",i+1
    for j in range(0,NCL.size/3):
      adj=Nneighbours(j,NCL,Tri,1)
      NCL[j,]=np.sum(NCL[adj,],0)/(adj.size)
  return NCL

def RBFmorph(NCB,knownC,knownD):
  # Determine displacement of Base nodes NCB by displacing nodes at coordinates 
  # knownC by known landmark displacement
  ####rbfx = Rbf(knownC[:,0],knownC[:,1],knownC[:,2],knownD[:,0])
  ####rbfy = Rbf(knownC[:,0],knownC[:,1],knownC[:,2],knownD[:,1])
  ####rbfz = Rbf(knownC[:,0],knownC[:,1],knownC[:,2],knownD[:,2])
  ####DSx = rbfx(NCB[:,0],NCB[:,1],NCB[:,2])
  ####DSy = rbfy(NCB[:,0],NCB[:,1],NCB[:,2])
  ####DSz = rbfz(NCB[:,0],NCB[:,1],NCB[:,2])
  ####NCI = NCB + np.c_[DSx,DSy,DSz]
  nrB = np.array([knownC.shape])[0,0]
  shapeN = np.array([NCB.shape])
  nrN = shapeN[0,0]
  MATR = np.zeros((nrB,nrB))
  knownD = np.r_[knownD,np.zeros((4,3))]
  Pb = np.c_[np.ones((nrB,1)),knownC]

  for i in range(0,nrB):
    # Thin plate spline:
    TEMP = np.array([np.sqrt(np.sum(np.power((np.ones((nrB,1))*knownC[i,:]-knownC),2),1))])
    TEMP = TEMP*TEMP*np.log10(TEMP)
    MATR[i,:] = TEMP.reshape((1,nrB))
    MATR[i,i] = 0

  MATR = np.c_[np.r_[MATR,Pb.T],np.r_[Pb,np.zeros((4,4))]]
  Coeff = np.array(np.linalg.solve(MATR,knownD))
  ResDisp = np.zeros((nrN,3))

  for i in range(0,nrN):
    TEMP = np.sqrt(np.sum(np.power((np.ones((nrB,1))*NCB[i,:]-knownC),2),1))
    TEMPlog10 = np.log10(TEMP)
    TEMPlog10 = sp.select([TEMPlog10<-1000,TEMPlog10>=0],[-1000,TEMPlog10])
    TEMP = TEMP*TEMP*TEMPlog10
    ResDisp[i,:] = np.sum(Coeff*np.r_[TEMP.reshape((nrB,1))*np.ones((1,3)),np.array([[1,1,1]]),np.transpose(np.ones((3,1))*NCB[i,:])],0)
  
  NCI = NCB+ResDisp
  return NCI


def Ndisp(nc,ncl,NCT,triT,Tcentr,KDT_TN):
  # Determine nodal displacement for node on Base to Target closest triangle
  # inputs:
  # nc: nodal coordinates on base mesh	ncl:target triangle	NCB&NCT: nodal coordinates of base and target mesh
  # triT: target connectivity 	Tcentr: centroids of target triangles
  # KDT_TN: kd-tree representation of target nodes
  # determine intersection point on closest triangle:
  #		perpendicular nodal vector from node to plane in which triangle lies:
  v0 = NCT[triT[ncl,1],:]-NCT[triT[ncl,0],:]
  v1 = NCT[triT[ncl,2],:]-NCT[triT[ncl,0],:]
  nvec = np.cross(v0,v1)
  nvec = nvec/np.sqrt(np.sum(nvec*nvec))
  # 	Intersection point:
  G = nc+np.dot(nvec,Tcentr[ncl,:]-nc)*nvec
  #		Check for G inside triangle ncl (done using Barycentric technique)
  v2 = G-NCT[triT[ncl,0],:]
  if barycent(v0,v1,v2):
    Disp = G-nc
  else:
    ndist,ncl = KDT_TN.query(nc)
    Disp = NCT[ncl,:]-nc
  return Disp



def barycent(v0,v1,v2):	
  # check wheter point is inside a triangle using vectors 
  #(first 2 constructed from triangle points and third from point to check and "base node" on triangle
  dot00 = np.dot(v0,v0)
  dot01 = np.dot(v0,v1)
  dot02 = np.dot(v0,v2)
  dot11 = np.dot(v1,v1)
  dot12 = np.dot(v1,v2)
  # Barycentric coordinates:
  invDenom = 1/(dot00*dot11 - dot01*dot01)
  u = (dot11*dot02-dot01*dot12)/invDenom
  v = (dot00*dot12-dot01*dot02)/invDenom
  return (u>0)&(v>0)&(u+v<1)
