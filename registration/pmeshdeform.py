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

LIM=24

def elasticsurf(NCB,ConnectB,LandmB,LandmB_NC,AllowableBI,NCT,ConnectT,AllowableT,UseN_B,UseN_T,k_max,USENORMALS,gamm=2,sigm0=10,f=1.0715):
  # Elastic surface registration:
  #inputs:
  # NCB,NCT: nodal coordinates of base and target surfaces
  # ConnectB;ConnectT: Base&target connectivity
  # LandmB,LandmB_NC landmarks that have to have a 1 to 1 correspondence (input 0 no landmarks are present)
  # UseN_B & AllowableB: Feature dependant nodes on Base-mesh (indices in NCB) and allowable triangles to match.
  # UseN_T & AllowableT: Selective Feature preserving nodes and triangles (indices in NCT and ConnectT) on target mesh.
  # k_max: maximum number of iterations
  ######## ADDITIONAL SETTINGS REQUIRED ARE SET INTERNAL TO CODE #########
  print
  print "SELECTIVE MESH MORPHING ALGORITHM USING ELASTIC SURFACE REGISTRATION"
  print "	-G.J.J.v.Rensburg - 22/04/2010-"
  t_start = time.clock()
  ConnectB = np.array(ConnectB,int)
  ConnectT = np.array(ConnectT,int)
  #LandmB = np.array(LandmB[:,0],int)		# do -1 later to be consistent with python indexing, first need to do other "temporary landmarks"& check that they dont fall on actual landmark positions!

  # Settings for elastic surface registration:
  m=20 # nearest neighbour parameter
  alph=0.5 # normilization factor
  #gamm=2 # smoothing parameter1
  #sigm0=10 # smoothing parameter2
  #f=1.0715 # smoothing parameter3
  Tol=0.0001 # stopping criteria

  # determine N1,N2,T1 and T2:
  N1=NCB.shape[0]
  N2=NCT.shape[0]
  T1=ConnectB.shape[0]
  T2=ConnectT.shape[0]
  NL=LandmB.shape[0]
  # For parallel programming devide Nr of computations by number of parallel processes (LIM)
  NPP1 = N1/LIM
  NPP2 = N2/LIM

  ################################     INITIALIZE & NODES OF CONCERN:    #################################
  ########################################################################################################
  print
  print
  print "Set up 1-ring neighbor list for all points on the generic mesh"
  #neighbList = [[0]]*N1
  #results = pprocess.Map(limit=LIM)
  #calc = results.manage(pprocess.MakeParallel(Get1neigh))
  #for j in range(0,LIM):
      #calc(np.array(range(0,NPP1))+j*NPP1,NCB,ConnectB)
  #for j in range(0,LIM):
      #neighbList[j*NPP1:(1+j)*NPP1] = results[j]
  #neighbList[LIM*NPP1:N1]=Get1neigh(np.array(range(LIM*NPP1,N1)),NCB,ConnectB)
  #np.ma.dump(neighbList,'SkullSurf_neighbList')
  neighbList = np.ma.load('SkullSurf_neighbList')
  
  print
  print "INITIALIZE SURFACE DEFORMATION"
  CONV = []
  print " 	enquire nodes where required displacement is checked"
  ###remove Landmarks from FDNB and SFPNT:
  #for i in range(0,NL):
    #if find_repeats(np.r_[UseN_B,LandmB[i,]])[0].size>0:
      #r=np.where(UseN_B==LandmB[i,])[0]
      #UseN_B = np.r_[UseN_B[0:r,],UseN_B[r+1:UseN_B.size,]]
  SamplingB=UseN_B.size
  SamplingT=UseN_T.size
  ## Full list of nodes used in Surface registration:
  LMB = np.r_[UseN_B]#,LandmB]	# Last NL entries are reserved for Landmarks that HAVE TO FIT points on the target mesh
  LMT = np.r_[UseN_T]


  # For parallel programming devide Nr of computations by number of parallel processes (LIM)
  SBPP = SamplingB/LIM
  STPP = SamplingT/LIM
  FMorph = 0

  
  print	
  print "COARSE SURFACE REGISTRATION"
  #print "	Compute known displacement for Base_Landmarks "
  #knownC = NCB[LandmB,]
  #knownD = LandmB_NC-knownC
  ####print "	using landmark displacements to deform using RBF"
  ####W_km1 = RBFmorph(NCB,knownC,knownD)
  ####tic = time.clock()
  ####W_km1 = MeshSmooth(W_km1,neighbList,10)
  ####print "		Smoothing done in ",time.clock()-tic," seconds"
  ####np.ma.dump(W_km1,'TempElasNodes_Iter'+str(k-1)+'_Time'+time.ctime())
  #print 'Smooth Gaussian Weight deformation to align Landmarks to target positions'
  #k=0
  #Err = 2
  #W_km1 = np.r_[NCB]
  #while (k<100)|(Err>Tol):
    #k=k+1
    #print 'Iteration : ',k
    #DS = np.zeros((N1,3))
    #knownC = W_km1[LandmB,]
    #knownD = LandmB_NC-knownC
    #knownD[np.isnan(knownD)]=0
    ## Deform mesh using Gaussian smoothing as suggested in paper by R.Bryan et al.
    #sigma_k2 = np.power(np.power(f,-k)*20,2)
    #results = pprocess.Map(limit=LIM)
    #calc = results.manage(pprocess.MakeParallel(GaussianSmooth))
    #for j in range(0,LIM):
      #calc(np.array(range(0,NPP1))+j*NPP1,W_km1,knownC,knownD,sigma_k2,gamm)
    #for j in range(0,LIM):
      #DS[np.array(range(0,NPP1))+j*NPP1,:] = results[j]
    #DS[range(LIM*NPP1,N1),:]=GaussianSmooth(np.array(range(LIM*NPP1,N1)),W_km1,knownC,knownD,sigma_k2,gamm)
    #DS[np.isnan(DS)]=0
    #W_km1 = W_km1+DS
    #Err = np.sum(np.sqrt(np.sum(DS*DS,1)),0)/N1
  #W_km1 = MeshSmooth(W_km1,neighbList,10)
  #np.ma.dump(W_km1,'TempElasNodes_Iter0_TimeWedMar14_2011_20')
  ###np.ma.dump(W_km1,'TempElasNodes_Iter0_Time'+time.ctime())
  W_km1 = NCB

  ################################    MAIN MESH DEFORMATION ALGORITHM:   #################################
  ########################################################################################################
  k=1
  print
  print "ELASTIC SURFACE REGISTRATION"
  print "determine vertex normals of target surface"
  #Compute target-mesh triangle centroids:
  print "determining centroids of target surface triangles"
  S_2_centr = np.c_[np.sum(np.c_[NCT[ConnectT[:,0],0],NCT[ConnectT[:,1],0],NCT[ConnectT[:,2],0]],1)/3,
    np.sum(np.c_[NCT[ConnectT[:,0],1],NCT[ConnectT[:,1],1],NCT[ConnectT[:,2],1]],1)/3,np.sum(np.c_[NCT[ConnectT[:,0],2],NCT[ConnectT[:,1],2],NCT[ConnectT[:,2],2]],1)/3]
  print "determine triangle and vertex normals of target surface"
  TNORMT = np.cross(NCT[ConnectT[:,1],:]-NCT[ConnectT[:,0],:],NCT[ConnectT[:,2],:]-NCT[ConnectT[:,0],:])
  TNORMT = (TNORMT.T/(np.ones((3,1))*np.sqrt(np.sum(np.array([TNORMT*TNORMT]),2)))).T
  VNORMT = vrtxnormal(NCT,ConnectT,S_2_centr,TNORMT)
 
  print "determining kd-trees of target surface centroids and nodal coordinates"
  KDT_TC = KDTree(S_2_centr,m)
  KDT_TN = KDTree(NCT,m)
  
  print 'initialize absolute Gaussian weight for final displacement to preserve element quality'
  GW = np.ones((SamplingB+SamplingT,1))
  
  while k<=k_max:
    D1 = np.zeros((SamplingB,3))
    D2 = np.zeros((SamplingT,3))
    DS = np.zeros((N1,3))
    AllowableB = np.r_[AllowableBI]
    print
    print "MESH DEFORMATION ITERATION",k
    print "	determining known displacement of landmarks"
    if NL>0:
      knownD = LandmB_NC-W_km1[LandmB,]
    print "	determining centroids of deforming mesh"
    W_km1_centr = np.c_[np.sum(np.c_[W_km1[ConnectB[:,0],0],W_km1[ConnectB[:,1],0],W_km1[ConnectB[:,2],0]],1)/3,
      np.sum(np.c_[W_km1[ConnectB[:,0],1],W_km1[ConnectB[:,1],1],W_km1[ConnectB[:,2],1]],1)/3,np.sum(np.c_[W_km1[ConnectB[:,0],2],W_km1[ConnectB[:,1],2],W_km1[ConnectB[:,2],2]],1)/3]
    print "	determine triangle and vertex normals of deforming surface"
    TNORMB = np.cross(W_km1[ConnectB[:,1],:]-W_km1[ConnectB[:,0],:],W_km1[ConnectB[:,2],:]-W_km1[ConnectB[:,0],:])
    TNORMB = (TNORMB.T/(np.ones((3,1))*np.sqrt(np.sum(np.array([TNORMB*TNORMB]),2)))).T
    VNORMB = vrtxnormal(W_km1,ConnectB,W_km1_centr,TNORMB)

    print "	determining kd-tree of current deforming surface centroids and nodal coordinates"
    KDT_KC = KDTree(W_km1_centr,m)
    KDT_KN = KDTree(W_km1,m)
    #if find_repeats(np.r_[USENORMALS,k])[0].size>0:
    #print " ###	Use triangle and vertex normals in setting up point correspondence"
    print "		setting up D1(i,d)"
    tic = time.clock()
    results = pprocess.Map(limit=LIM)
    calc = results.manage(pprocess.MakeParallel(DsetupNorm))
    for j in range(0,LIM):
      calc(np.array(range(0,SBPP))+j*SBPP,W_km1,VNORMB,NCT,TNORMT,VNORMT,ConnectT,S_2_centr,AllowableT,LMB,D1)
    for j in range(0,LIM):
      D1[np.array(range(0,SBPP))+j*SBPP,:] = results[j]
    D1[range(LIM*SBPP,SamplingB),:]=DsetupNorm(range(LIM*SBPP,SamplingB),W_km1,VNORMB,NCT,TNORMT,VNORMT,ConnectT,S_2_centr,AllowableT,LMB,D1)
    #D1=np.r_[D1,knownD]
    print "			",time.clock()-tic," seconds"
    print "		update allowable triangles on generic mesh:"
    remP = D1[:,0]+D1[:,1]+D1[:,2]==0
    removeP = LMB[remP]
    print "			unregistered points on generic mesh: ",removeP.size
    print "			number of original generic triangles allowed: ",AllowableB.shape[0]
    for rp in removeP:
      rowsNo = np.where(AllowableB==rp)[0]
      rowsNo.sort
      for rr in rowsNo[::-1]:
	AllowableB = AllowableB[np.where(range(AllowableB.shape[0])<>rr)[0],]
    print "			number of generic triangles allowed for current iteration: ",AllowableB.shape[0]
    if find_repeats(np.r_[USENORMALS,k])[0].size>0:
      print " ###	Use triangle and vertex normals in setting up point correspondence"
      print "		setting up D2(j,c)"
      tic = time.clock()
      results = pprocess.Map(limit=LIM)
      calc = results.manage(pprocess.MakeParallel(DsetupNorm))
      for j in range(0,LIM):
	calc(np.array(range(0,STPP))+j*STPP,NCT,VNORMT,W_km1,TNORMB,VNORMB,ConnectB,W_km1_centr,AllowableB,LMT,D2)
      for j in range(0,LIM):
	D2[np.array(range(0,STPP))+j*STPP,:] = results[j]
      D2[range(LIM*STPP,SamplingT),:]=DsetupNorm(range(LIM*STPP,SamplingT),NCT,VNORMT,W_km1,TNORMB,VNORMB,ConnectB,W_km1_centr,AllowableB,LMT,D2)
      print "			",time.clock()-tic," seconds"
    else:
      print "	Simple closest point search iteration "
      #print "		setting up D1(i,d)"
      #tic = time.clock()
      #results = pprocess.Map(limit=LIM)
      #calc = results.manage(pprocess.MakeParallel(Dsetup))
      #for j in range(0,LIM):
	#calc(np.array(range(0,SBPP))+j*SBPP,W_km1,NCT,ConnectT,S_2_centr,AllowableT,LMB,D1,KDT_TC,KDT_TN)
      #for j in range(0,LIM):
	#D1[np.array(range(0,SBPP))+j*SBPP,:] = results[j]
      #D1[range(LIM*SBPP,SamplingB),:]=Dsetup(range(LIM*SBPP,SamplingB),W_km1,NCT,ConnectT,S_2_centr,AllowableT,LMB,D1,KDT_TC,KDT_TN)
      ##D1=np.r_[D1,knownD]
      #print "			",time.clock()-tic," seconds"
      #remP = D1[:,0]+D1[:,1]+D1[:,2]==0
      #removeP = LMB[remP]
      #print "			unregistered points on generic mesh: ",removeP.size
      #print "			number of original generic triangles allowed: ",AllowableB.shape[0]
      #for rp in removeP:
	#rowsNo = np.where(AllowableB==rp)[0]
	#rowsNo.sort
	#for rr in rowsNo[::-1]:
	  #AllowableB = AllowableB[np.where(range(AllowableB.shape[0])<>rr)[0],]
      #print "			number of generic triangles allowed for current iteration: ",AllowableB.shape[0]
      print "		setting up D2(j,c)"
      tic = time.clock()
      results = pprocess.Map(limit=LIM)
      calc = results.manage(pprocess.MakeParallel(Dsetup))
      for j in range(0,LIM):
	calc(np.array(range(0,STPP))+j*STPP,NCT,W_km1,ConnectB,W_km1_centr,AllowableB,LMT,D2,KDT_KC,KDT_KN)
      for j in range(0,LIM):
	D2[np.array(range(0,STPP))+j*STPP,:] = results[j]
      D2[range(LIM*STPP,SamplingT),:]=Dsetup(range(LIM*STPP,SamplingT),NCT,W_km1,ConnectB,W_km1_centr,AllowableB,LMT,D2,KDT_KC,KDT_KN)
      print "			",time.clock()-tic," seconds"
   
    # Compute displacement update for each node using suggested Gaussian radial basis function:
    print "	determining smoothed displacement field"
    
    tic=time.clock()
    NCp = np.r_[W_km1[LMB,:],NCT[LMT,:]+D2]
    DD = np.r_[D1,-D2]
    # Mask Nan and Inf values if any:
    DD[np.isnan(DD)]=0
    DD[np.isinf(DD)]=0
    #keepP = DD[:,0]+DD[:,1]+DD[:,2]<>0
    #print keepP
    #NCp,DD = NCp[keepP,:],DD[keepP,:]
    #KDTp = KDTree(NCp,5)
    # Deform mesh using Gaussian smoothing as suggested in paper by R.Bryan et al.
    sigma_k2 = np.power(np.power(f,-k)*sigm0,2)
    results = pprocess.Map(limit=LIM)
    calc = results.manage(pprocess.MakeParallel(GaussianSmooth))
    for j in range(0,LIM):
      calc(np.array(range(0,NPP1))+j*NPP1,W_km1,NCp,DD,sigma_k2,gamm)
    for j in range(0,LIM):
      DS[np.array(range(0,NPP1))+j*NPP1,:] = results[j]
    DS[range(LIM*NPP1,N1),:]=GaussianSmooth(np.array(range(LIM*NPP1,N1)),W_km1,NCp,DD,sigma_k2,gamm)
    print "			",time.clock()-tic," seconds"

    
    # Mask Nan and Inf if any:
    DS[np.isnan(DS)]=0
    DS[np.isinf(DS)]=0
    
    #print 'Check if current iteration reduces element quality to below allowable and stiffen mesh accordingly'
    print
    print
    print 'Convergence History'
    print CONV
    print
    print
    
    # Determine Jacobian of all elements and if unsattisfied apply stiffening (Decrease GW <1) untill this doesn't happen
    
    # determine wheter convergence is acheived
    #TotalMorph = np.sum(np.sqrt(np.sum(DS*DS,1)),0)/NCB.shape[0]    
    TotalMorph = np.sum(np.sqrt(np.sum(DD*DD,1)))/(DD.size/3)
    CONV = CONV + [TotalMorph]
    FMorph = (k==1)*TotalMorph+FMorph
    print "	average nodal displacement for current deformation iteration:"
    print TotalMorph
    if (TotalMorph<Tol):
      print
      print "CONVERGED SOLUTION OBTAINED"
      #CONV = CONV + [TotalMorph]
      k = k_max*10+1
      W_km1 = W_km1+DS
    elif (k<10)|(TotalMorph<10*FMorph):
      print "problem not yet converged at iteration",k
      #CONV = CONV + [TotalMorph]
      k = k+1;
      # Deform mesh:
      print "	deforming mesh (update of W_{k-1})"
      W_km1 = W_km1+DS
      #np.ma.dump(W_km1,'Femur2NC_'+str(k))
    else:
      print "PROBLEM DIVERGING"
      k=k_max*10-1
    
    #np.ma.dump(W_km1,'TempElasNodes_Iter'+str(k-1)+'_Time'+time.ctime())
    if (k>2)&(np.mod(k-1,5)==0):
      print
      #np.ma.dump(W_km1,'TempElasNodes_Iter'+str(k-1)+'_Time'+time.ctime())
      #W_km1 = RBFmorph(W_km1,W_km1[LandmB,],LandmB_NC-W_km1[LandmB,])
      tic = time.clock()
      W_km1 = MeshSmooth(W_km1,neighbList,10)
      np.ma.dump(W_km1,'SkullUnique2_gamm'+str(gamm)+'_sigN'+str(sigm0)+'_iter'+str(k-1))
      print "		Smoothing done in ",time.clock()-tic," seconds"
    #print "COARSE SURFACE REGISTRATION"
    #print "	using landmark displacements to deform using RBF"
    #W_km1 = RBFmorph(W_km1,W_km1[LandmB,],LandmB_NC-W_km1[LandmB,])
  print
      
  if k==k_max+1:
    print
    print "SOLUTION TERMINATED: maximum iterations,(",k_max,") reached"
  print
  print "TOTAL TIME FOR ELASTIC SURFACE REGISTRATION : ",time.clock()-t_start,"seconds"
  CONV = np.array(CONV)
  return W_km1,CONV

def DsetupNorm(places,NCB,VNORMB,NCT,TNORMT,VNORMT,ConnectT,CentrTT,AllowableT,LMB,D1): 
  for i in places:
    nn=LMB[i]
    keep=np.where((np.dot(TNORMT,VNORMB[nn,].reshape((3,)))>0.2))[0]
    if keep.size>0:
      distA = cdist(CentrTT[keep,],NCB[nn,].reshape((1,3)),'euclidean')
      ncl = keep[np.where(distA==np.min(distA))[0]][0]
      if np.where(AllowableT==ncl)[0].size>0:
	D1[i,:] = NdispNORM(NCB[nn,],VNORMB[nn,].reshape((3,)),ncl,NCT,ConnectT,CentrTT,VNORMT)
  return D1[places,:]
  
def Dsetup(places,NCB,NCT,ConnectT,CentrTT,AllowableT,LMB,D1,KDT_TC,KDT_TN):
  for i in places:
    nn=LMB[i]
    # query kd-tree for closest triangle to node:
    ncl = KDT_TC.query(NCB[nn,])[1]
    # check if 
    if np.where(AllowableT==ncl)[0].size>0:
      D1[i,:] = Ndisp(NCB[nn,],ncl,NCT,ConnectT,CentrTT,KDT_TN)
  return D1[places,:]

#def MeshSmooth(nodes,NCL,Tri):
  #for i in nodes:
    #adj=Nneighbours(i,NCL,Tri,1)[1]
    #NCL[i,:] = np.sum(NCL[adj,],0)/(adj.size)
  #return NCL[nodes,:]
def MeshSmooth(NC,neighbList,Iter):
  print "Two-stage Taubin Smoothing"
  N1 = NC.shape[0]
  NPP1 = N1/LIM
  for i in range(Iter):
    print "	Iteration ",i+1
    print "		Umbrella-Operator Step"
    Umb = -np.zeros((NC.shape))
    results = pprocess.Map(limit=LIM)
    calc = results.manage(pprocess.MakeParallel(MeshUmbrella))
    for j in range(0,LIM):
	calc(np.array(range(0,NPP1))+j*NPP1,NC,neighbList,Umb)
    for j in range(0,LIM):
	Umb[np.array(range(0,NPP1))+j*NPP1,:] = results[j]
    Umb[range(LIM*NPP1,N1),:]=MeshUmbrella(range(LIM*NPP1,N1),NC,neighbList,Umb)
    print "		Smoothing Step"
    results = pprocess.Map(limit=LIM)
    calc = results.manage(pprocess.MakeParallel(MeshTaubin))
    for j in range(0,LIM):
	calc(np.array(range(0,NPP1))+j*NPP1,NC,neighbList,Umb)
    for j in range(0,LIM):
	NC[np.array(range(0,NPP1))+j*NPP1,:] = results[j]
    NC[range(LIM*NPP1,N1),:]=MeshTaubin(range(LIM*NPP1,N1),NC,neighbList,Umb)
  return NC
    
    
def Get1neigh(nodes,NC,Tri):
  Nz = nodes.size
  neighbSeg = [[0]]*Nz
  for i in range(Nz):
    neighbSeg[i]=Nneighbours(nodes[i],NC,Tri,1)[1]
  return neighbSeg
  
def MeshUmbrella(nodes,NC,neighbList,Umb):
  for i in nodes:
    adj=neighbList[i]
    Umb[i,:] = np.sum(NC[adj,],0)/(adj.size) - NC[i,]
  return Umb[nodes,:]

def MeshTaubin(nodes,NC,neighbList,Umb):
  for i in nodes:
    adj=neighbList[i]
    NC[i,:] = NC[i,:] + Umb[i,:]*0.235 -0.265*(np.sum(Umb[adj,],0))/(adj.size)
  return NC[nodes,:]

#def GaussianSmooth(nodes,W_km1,KDTp,DD,sigma_k2,gamm):
  #DSs=np.zeros((nodes.size,3))
  #for j in range(nodes.size):
    #i = nodes[j]
    #G1node,neighb = KDTp.query(W_km1[i,:],200)
    #neighb=np.array(neighb,int)
    #G1node = (np.ones((3,1))*np.exp(-G1node.T/sigma_k2)).T
    #DSs[j,:] = np.sum(DD[neighb,:]*G1node,0)/(np.sum(G1node[:,0].T)*gamm)
  #return DSs#[nodes,:]


def GaussianSmooth(nodes,W_km1,NCp,DD,sigma_k2,gamm):
  DSs=np.zeros((nodes.size,3))
  for j in range(nodes.size):
    i = nodes[j]
    G1node = cdist(NCp,W_km1[i,:].reshape((1,3)),'euclidean')
    G1node = (np.ones((3,1))*np.exp(-G1node.T/sigma_k2)).T
    DSs[j,:] = np.sum(DD*G1node,0)/(np.sum(G1node[:,0].T)*gamm)
  return DSs#[nodes,:]


def RBFmorph(NCB,knownC,knownD):
  # Determine displacement of Base nodes NCB by displacing nodes at coordinates 
  # knownC by known landmark displacement
  #rbfx = Rbf(knownC[:,0],knownC[:,1],knownC[:,2],knownD[:,0],function='gaussian')
  #rbfy = Rbf(knownC[:,0],knownC[:,1],knownC[:,2],knownD[:,1],function='gaussian')
  #rbfz = Rbf(knownC[:,0],knownC[:,1],knownC[:,2],knownD[:,2],function='gaussian')
  #DSx = rbfx(NCB[:,0],NCB[:,1],NCB[:,2])
  #DSy = rbfy(NCB[:,0],NCB[:,1],NCB[:,2])
  #DSz = rbfz(NCB[:,0],NCB[:,1],NCB[:,2])
  #print np.c_[DSx,DSy,DSz]
  #print NCB
  #NCI = np.c_[DSx,DSy,DSz]+NCB
  nrB = np.array([knownC.shape])[0,0]
  shapeN = np.array([NCB.shape])
  nrN = shapeN[0,0]
  MATR = np.zeros((nrB,nrB))
  knownD = np.r_[knownD,np.zeros((4,3))]
  Pb = np.c_[np.ones((nrB,1)),knownC]
  for i in range(0,nrB):
    TEMP = np.array([np.sqrt(np.sum(np.power((np.ones((nrB,1))*knownC[i,:]-knownC),2),1))])
    ##### Thin plate spline:
    TEMP = TEMP*TEMP*np.log10(TEMP)
    ##### Gauss weight:
    #TEMP = np.exp(-TEMP*TEMP)
    ##### CPC2
    #TEMP = np.power((1-TEMP),4)*(4*TEMP+1)
    MATR[i,:] = TEMP.reshape((1,nrB))
    MATR[i,i] = 0
  MATR = np.c_[np.r_[MATR,Pb.T],np.r_[Pb,np.zeros((4,4))]]
  #MATR = np.matrix(MATR)
  #MATR = np.linalg.pinv(MATR)
  #knownD = np.matrix(knownD)
  Coeff = np.array(np.linalg.solve(MATR,knownD))
  ResDisp = np.zeros((nrN,3))
  for i in range(0,nrN):
    TEMP = np.sqrt(np.sum(np.power((np.ones((nrB,1))*NCB[i,:]-knownC),2),1))
    ##### Thin Plate Spline:
    TEMPlog10 = np.log10(TEMP)
    TEMPlog10 = sp.select([TEMPlog10<-1000,TEMPlog10>=0],[-1000,TEMPlog10])
    TEMP = TEMP*TEMP*TEMPlog10
    ##### Gauss weight:
    #TEMP = np.exp(-TEMP*TEMP)
    ##### CPC2
    #TEMP = np.power((1-TEMP),4)*(4*TEMP+1)
    
    ResDisp[i,:] = np.sum(Coeff*np.r_[TEMP.reshape((nrB,1))*np.ones((1,3)),np.array([[1,1,1]]),np.transpose(np.ones((3,1))*NCB[i,:])],0)
  NCI = NCB+ResDisp
  return NCI

def NdispNORM(nc,nnNorm,ncl,NCT,triT,Tcentr,VNORMT):
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
  v2 = G-NCT[triT[ncl,0],]
  Disp = np.array([0,0,0])
  if barycent(v0,v1,v2):
    Disp = G-nc
  else:	
    #keep=np.where(np.dot(VNORMT,nnNorm)>0.5)[0]
    ### Get closest node in the targeted triangle
    keep = triT[ncl,]
    if keep.size>0:
      distA = cdist(NCT[keep,],nc.reshape((1,3)),'euclidean')
      ncl = keep[np.where(distA==np.min(distA))[0]][0]
      Disp = NCT[ncl,]-nc
  return Disp
  
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
  v2 = G-NCT[triT[ncl,0],]
  if barycent(v0,v1,v2):
    Disp = G-nc
  else:
    ncl = KDT_TN.query(nc)[1]
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
