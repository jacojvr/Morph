# -*- coding: utf-8 -*-

# Iterative closest point to perform rigid body transformations using line registration
import numpy as np
import scipy as sp
#import mdp
from scipy.optimize import *
from scipy.stats.mstats import find_repeats
from scipy.spatial import KDTree
from scipy.spatial.distance import cdist
import pprocess
from ptrimeshalt import *

LIM=24
#NC_B,LM_B,NC_T,LM_T = [],[],[],[]
#[RsegmB,VsegmB,RsegmT,VsegmT] = [np.array([[],[],[]]).reshape((0,3)),np.array([[],[],[]]).reshape((0,3)),
  #np.array([[],[],[]]).reshape((0,3)),np.array([[],[],[]]).reshape((0,3))]
#[RnodesB,VnodesB,RnodesT,VnodesT] = [np.array([[],[]]).reshape((0,2)),np.array([[],[]]).reshape((0,2)),
  #np.array([[],[]]).reshape((0,2)),np.array([[],[]]).reshape((0,2))]
def LineReg(NCB,TriB,RlinesB,VlinesB,NCT,TriT,RlinesT,VlinesT,percentReg,DistMax):
  N1 = NCB.shape[0]
  NPP1 = N1/LIM
  [RsegmB,VsegmB,RsegmT,VsegmT] = [np.array([[],[],[]]).reshape((0,3)),np.array([[],[],[]]).reshape((0,3)),
    np.array([[],[],[]]).reshape((0,3)),np.array([[],[],[]]).reshape((0,3))]
  [RnodesB,VnodesB,RnodesT,VnodesT] = [np.array([[],[]]).reshape((0,2)),np.array([[],[]]).reshape((0,2)),
    np.array([[],[]]).reshape((0,2)),np.array([[],[]]).reshape((0,2))]
  # Transform Lines into linesegments with i'th segment allocated as [NC_point1,NC_point2,Line_nr]
  #	and list of nodes allocated as [Nd_nr, Line_nr]
  for i in range(1,RlinesB[0]+1):
    Lsize=RlinesB[i].size
    RnodesB = np.array(np.r_[RnodesB,np.c_[RlinesB[i],np.ones((Lsize,1))*i]],int)
    RsegmB = np.array(np.r_[RsegmB,np.c_[RlinesB[i][0:Lsize-1],RlinesB[i][1:Lsize],np.ones((Lsize-1,1))*i]],int)
  for i in range(1,VlinesB[0]+1):
    Lsize=VlinesB[i].size
    VnodesB = np.array(np.r_[VnodesB,np.c_[VlinesB[i],np.ones((Lsize,1))*i]],int)
    VsegmB = np.array(np.r_[VsegmB,np.c_[VlinesB[i][0:Lsize-1],VlinesB[i][1:Lsize],np.ones((Lsize-1,1))*i]],int)
  for i in range(1,RlinesT[0]+1):
    Lsize=RlinesT[i].size
    RnodesT = np.array(np.r_[RnodesT,np.c_[RlinesT[i],np.ones((Lsize,1))*i]],int)
    RsegmT = np.array(np.r_[RsegmT,np.c_[RlinesT[i][0:Lsize-1],RlinesT[i][1:Lsize],np.ones((Lsize-1,1))*i]],int)
  for i in range(1,VlinesT[0]+1):
    Lsize=VlinesT[i].size
    VnodesT = np.array(np.r_[VnodesT,np.c_[VlinesT[i],np.ones((Lsize,1))*i]],int)
    VsegmT = np.array(np.r_[VsegmT,np.c_[VlinesT[i][0:Lsize-1],VlinesT[i][1:Lsize],np.ones((Lsize-1,1))*i]],int)
  # find average nodal coordinate of base linesegments
  RsegTNC,VsegTNC = (NCT[RsegmT[:,0],]+NCT[RsegmT[:,1],])/2,(NCT[VsegmT[:,0],]+NCT[VsegmT[:,1],])/2
  # set up k-d tree of nodes and segments in Base model lines
  print 'TARGET: Determine triangle centroids, triangle normals and weighted vertex normals'
  TBCT = np.c_[np.sum(np.c_[NCT[TriT[:,0],0],NCT[TriT[:,1],0],NCT[TriT[:,2],0]],1)/3,
    np.sum(np.c_[NCT[TriT[:,0],1],NCT[TriT[:,1],1],NCT[TriT[:,2],1]],1)/3,np.sum(np.c_[NCT[TriT[:,0],2],NCT[TriT[:,1],2],NCT[TriT[:,2],2]],1)/3]
  TNORMT = np.cross(NCT[TriT[:,1],:]-NCT[TriT[:,0],:],NCT[TriT[:,2],:]-NCT[TriT[:,0],:])
  TNORMT = (TNORMT.T/(np.ones((3,1))*np.sqrt(np.sum(np.array([TNORMT*TNORMT]),2)))).T
  VNORMT = vrtxnormal(NCT,TriT,TBCT,TNORMT)
  print 'TARGET: Determine segment normal directions'
  SegNormRT,SegNormVT = (VNORMT[RsegmT[:,0],]+VNORMT[RsegmT[:,1],]),(VNORMT[VsegmT[:,0],]+VNORMT[VsegmT[:,1],])
  SegNormRT = (SegNormRT.T/(np.ones((3,1))*np.sqrt(np.sum(np.array([SegNormRT*SegNormRT]),2)))).T
  SegNormVT = (SegNormVT.T/(np.ones((3,1))*np.sqrt(np.sum(np.array([SegNormVT*SegNormVT]),2)))).T
  
  
  # find average nodal coordinate of base linesegments
  RsegBNC,VsegBNC = (NCB[RsegmB[:,0],]+NCB[RsegmB[:,1],])/2,(NCB[VsegmB[:,0],]+NCB[VsegmB[:,1],])/2
  # set up k-d tree of nodes and segments in Base model lines
  print 'BASE: Determine triangle centroids, triangle normals and weighted vertex normals'
  TBCB = np.c_[np.sum(np.c_[NCB[TriB[:,0],0],NCB[TriB[:,1],0],NCB[TriB[:,2],0]],1)/3,
    np.sum(np.c_[NCB[TriB[:,0],1],NCB[TriB[:,1],1],NCB[TriB[:,2],1]],1)/3,np.sum(np.c_[NCB[TriB[:,0],2],NCB[TriB[:,1],2],NCB[TriB[:,2],2]],1)/3]
  TNORMB = np.cross(NCB[TriB[:,1],:]-NCB[TriB[:,0],:],NCB[TriB[:,2],:]-NCB[TriB[:,0],:])
  TNORMB = (TNORMB.T/(np.ones((3,1))*np.sqrt(np.sum(np.array([TNORMB*TNORMB]),2)))).T
  VNORMB = vrtxnormal(NCB,TriB,TBCB,TNORMB)
  print 'BASE: Determine segment normal directions'
  SegNormRB,SegNormVB = (VNORMB[RsegmB[:,0],]+VNORMB[RsegmB[:,1],]),(VNORMB[VsegmB[:,0],]+VNORMB[VsegmB[:,1],])
  SegNormRB = (SegNormRB.T/(np.ones((3,1))*np.sqrt(np.sum(np.array([SegNormRB*SegNormRB]),2)))).T
  SegNormVB = (SegNormVB.T/(np.ones((3,1))*np.sqrt(np.sum(np.array([SegNormVB*SegNormVB]),2)))).T
  
  deform,k = 2,0
  while (deform>0.0001)&(k<100):
    Supp = 50./(1+k/5)
    if Supp<10:
      Supp=10
    k=k+1
    print '		ITERATION ',k
    ## find average nodal coordinate of base linesegments
    #RsegBNC,VsegBNC = (NCB[RsegmB[:,0],]+NCB[RsegmB[:,1],])/2,(NCB[VsegmB[:,0],]+NCB[VsegmB[:,1],])/2
    ## set up k-d tree of nodes and segments in Base model lines
    #print 'BASE: Determine triangle centroids, triangle normals and weighted vertex normals'
    #TBCB = np.c_[np.sum(np.c_[NCB[TriB[:,0],0],NCB[TriB[:,1],0],NCB[TriB[:,2],0]],1)/3,
      #np.sum(np.c_[NCB[TriB[:,0],1],NCB[TriB[:,1],1],NCB[TriB[:,2],1]],1)/3,np.sum(np.c_[NCB[TriB[:,0],2],NCB[TriB[:,1],2],NCB[TriB[:,2],2]],1)/3]
    #TNORMB = np.cross(NCB[TriB[:,1],:]-NCB[TriB[:,0],:],NCB[TriB[:,2],:]-NCB[TriB[:,0],:])
    #TNORMB = (TNORMB.T/(np.ones((3,1))*np.sqrt(np.sum(np.array([TNORMB*TNORMB]),2)))).T
    #VNORMB = vrtxnormal(NCB,TriB,TBCB,TNORMB)
    #print 'BASE: Determine segment normal directions'
    #SegNormRB,SegNormVB = (VNORMB[RsegmB[:,0],]+VNORMB[RsegmB[:,1],]),(VNORMB[VsegmB[:,0],]+VNORMB[VsegmB[:,1],])
    #SegNormRB = (SegNormRB.T/(np.ones((3,1))*np.sqrt(np.sum(np.array([SegNormRB*SegNormRB]),2)))).T
    #SegNormVB = (SegNormVB.T/(np.ones((3,1))*np.sqrt(np.sum(np.array([SegNormVB*SegNormVB]),2)))).T
    
    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~REGISTRATION~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
    print 'Register feature lines'
    R12R = FeatReg(RnodesB,NCB,VNORMB,RsegmT,RsegTNC,SegNormRT,RnodesT,NCT,DistMax) #find registered mesh 1 to 2 Ridges
    R12V = FeatReg(VnodesB,NCB,VNORMB,VsegmT,VsegTNC,SegNormVT,VnodesT,NCT,DistMax) #find registered mesh 1 to 2 Valleys
    R21R = FeatReg(RnodesT,NCT,VNORMT,RsegmB,RsegBNC,SegNormRB,RnodesB,NCB,DistMax) #find registered mesh 2 to 1 Ridges
    R21V = FeatReg(VnodesT,NCT,VNORMT,VsegmB,VsegBNC,SegNormVB,VnodesB,NCB,DistMax) #find registered mesh 2 to 1 Valleys
    print 'Build topological map to disregard unmatched features'
    R12R,R21R=BuildMap(NCB,NCT,RnodesB,RnodesT,RsegmB,RsegmT,R12R,R21R,percentReg)
    R12V,R21V=BuildMap(NCB,NCT,VnodesB,VnodesT,VsegmB,VsegmT,R12V,R21V,percentReg)
    print 'Determine Base deformation required'
    B12R,T12R = NCB[np.array(R12R[:,0],int),],R12R[:,1:4]
    B12V,T12V = NCB[np.array(R12V[:,0],int),],R12V[:,1:4]
    B21R,T21R = R21R[:,1:4],NCT[np.array(R21R[:,0],int),]
    B21V,T21V = R21V[:,1:4],NCT[np.array(R21V[:,0],int),]
    Base_i = np.r_[B12R,B12V,B21R,B21V]
    Target_i = np.r_[T12R,T12V,T21R,T21V]
    DISP = Target_i - Base_i 
    DISP_Full = np.zeros(NCB.shape)
    print 'Determine smooth defromation'
    results = pprocess.Map(limit=LIM)
    calc = results.manage(pprocess.MakeParallel(GaussMove))
    for j in range(0,LIM):
      calc(np.array(range(0,NPP1))+j*NPP1,NCB,Base_i,DISP,DISP_Full,Supp)
    for j in range(0,LIM):
      DISP_Full[np.array(range(0,NPP1))+j*NPP1,:] = results[j]
    DISP_Full[range(LIM*NPP1,N1),:]=GaussMove(range(LIM*NPP1,N1),NCB,Base_i,DISP,DISP_Full,Supp)
    NCB = NCB+DISP_Full
    deform = np.sum(np.sqrt(np.sum(DISP_Full*DISP_Full,1)))/N1
    print '	TOTAL CURRENT DEFORMATION: ',deform
  # return which NODES within surface features are registered:
  keepB=np.r_[R12R[:,0],R12V[:,0]]
  keepT=np.r_[R21R[:,0],R21V[:,0]]
  keepB,keepT=np.array(keepB,int),np.array(keepT,int)
  return NCB,keepB,keepT
  
def FeatReg(nodesB,NCB,normB,segmINFO,segNCT,segNormT,nodesT,NCT,DistMax):
  registeredP = np.zeros((nodesB.shape[0],6))
  for i in range(0,nodesB.shape[0]):
    # find closest line segments
    keep=np.where(np.dot(segNormT,normB[nodesB[i,0],].reshape((3,)))>0.5)[0]
    if keep.size>0:
      distA = cdist(segNCT[keep,],NCB[nodesB[i,0],].reshape((1,3)),'euclidean')
      segm = keep[np.where(distA==np.min(distA))[0]]
      a,b = NCT[segmINFO[segm,0],]-NCT[segmINFO[segm,1],],NCT[segmINFO[segm,1],]-NCB[nodesB[i,0],]
      alpha = -np.sum(a*b)/np.sum(a*a)
      NCR = alpha*NCT[segmINFO[segm,0],]+(1-alpha)*NCT[segmINFO[segm,1],]
      dist = np.sqrt(np.sum(alpha*alpha*a*a + 2*alpha*a*b + b*b))
      registeredP[i,]=np.r_[nodesB[i,0],NCR.reshape((3,)),dist,segm]
  registeredP = registeredP[np.where(registeredP[:,5]<>0)[0]]
  registeredP = registeredP[np.where(registeredP[:,4]<DistMax)[0]]
  #print '	initial:  ',registeredP.size/6
  #for i in nodesT[:,0]:	#remove doule registeres final coordinates by keeping only closest register
    #pos = np.where((registeredP[:,1]==NCT[i,0])&(registeredP[:,2]==NCT[i,1])&(registeredP[:,3]==NCT[i,2]))[0]
    #if pos.size>1:
      #keep = np.where(((registeredP[:,1]<>NCT[i,0])|(registeredP[:,2]<>NCT[i,1])|(registeredP[:,3]<>NCT[i,2]))|(registeredP[:,4]==np.min(registeredP[pos,4])))[0]
      #registeredP = registeredP[keep,]
  print '	Nr of registered points:  ',registeredP.size/6
  return registeredP
  
def lineRST(NCB,RlinesB,VlinesB,NCT,RlinesT,VlinesT,UseFeat=0,UseScale=0,Use1Scale=1):
# Use lines of curvature on two models to determine rigid body transformation for best fit
# Takes as input the nodal coordinates of the two surface meshes as well as the ridge and valley lines of the two.
# Target mesh is then rotated, scaled and translated to best fit the corresponding lines of curvature on the Base mesh
  #global RsegmB,VsegmB,RsegmT,VsegmT
  #global RnodesB,VnodesB,RnodesT,VnodesT
  if UseFeat == 1:
    [RsegmB,VsegmB,RsegmT,VsegmT] = [np.array([[],[],[]]).reshape((0,3)),np.array([[],[],[]]).reshape((0,3)),
      np.array([[],[],[]]).reshape((0,3)),np.array([[],[],[]]).reshape((0,3))]
    [RnodesB,VnodesB,RnodesT,VnodesT] = [np.array([[],[]]).reshape((0,2)),np.array([[],[]]).reshape((0,2)),
      np.array([[],[]]).reshape((0,2)),np.array([[],[]]).reshape((0,2))]
    # Transform Lines into linesegments with i'th segment allocated as [NC_point1,NC_point2,Line_nr]
    #	and list of nodes allocated as [Nd_nr, Line_nr]
    for i in range(1,RlinesB[0]+1):
      Lsize=RlinesB[i].size
      RnodesB = np.array(np.r_[RnodesB,np.c_[RlinesB[i],np.ones((Lsize,1))*i]],int)
      RsegmB = np.array(np.r_[RsegmB,np.c_[RlinesB[i][0:Lsize-1],RlinesB[i][1:Lsize],np.ones((Lsize-1,1))*i]],int)
    for i in range(1,VlinesB[0]+1):
      Lsize=VlinesB[i].size
      VnodesB = np.array(np.r_[VnodesB,np.c_[VlinesB[i],np.ones((Lsize,1))*i]],int)
      VsegmB = np.array(np.r_[VsegmB,np.c_[VlinesB[i][0:Lsize-1],VlinesB[i][1:Lsize],np.ones((Lsize-1,1))*i]],int)
    for i in range(1,RlinesT[0]+1):
      Lsize=RlinesT[i].size
      RnodesT = np.array(np.r_[RnodesT,np.c_[RlinesT[i],np.ones((Lsize,1))*i]],int)
      RsegmT = np.array(np.r_[RsegmT,np.c_[RlinesT[i][0:Lsize-1],RlinesT[i][1:Lsize],np.ones((Lsize-1,1))*i]],int)
    for i in range(1,VlinesT[0]+1):
      Lsize=VlinesT[i].size
      VnodesT = np.array(np.r_[VnodesT,np.c_[VlinesT[i],np.ones((Lsize,1))*i]],int)
      VsegmT = np.array(np.r_[VsegmT,np.c_[VlinesT[i][0:Lsize-1],VlinesT[i][1:Lsize],np.ones((Lsize-1,1))*i]],int)
    # find average nodal coordinate of base linesegments
    RsegBNC,VsegBNC = (NCB[RsegmB[:,0],]+NCB[RsegmB[:,1],])/2,(NCB[VsegmB[:,0],]+NCB[VsegmB[:,1],])/2
    # set up k-d tree of nodes and segments in Base model lines
    kdt_RBS,kdt_VBS = KDTree(RsegBNC,5),KDTree(VsegBNC,5)
  else:
    N1,N2 = NCB.shape[0],NCT.shape[0]
    NPP1 = N1/LIM
    NPP2 = N2/LIM
    kdt_Base = KDTree(NCB,20)
  diff=2
  k=0
  Conv = np.zeros((100,))
  while (diff>0.000001)&(k<100):
    k=k+1
    print '		ITERATION ',k
    if UseFeat == 1:
      RsegTNC,VsegTNC = (NCT[RsegmT[:,0],]+NCT[RsegmT[:,1],])/2,(NCT[VsegmT[:,0],]+NCT[VsegmT[:,1],])/2
      kdt_RTS,kdt_VTS = KDTree(RsegTNC,5),KDTree(VsegTNC,5)
      R12R = LineICP(RnodesB,NCB,RnodesT,NCT,kdt_RTS,RsegmT) #find registered mesh 1 to 2 Ridges
      R12V = LineICP(VnodesB,NCB,VnodesT,NCT,kdt_VTS,VsegmT) #find registered mesh 1 to 2 Valleys
      R21R = LineICP(RnodesT,NCT,RnodesB,NCB,kdt_RBS,RsegmB) #find registered mesh 2 to 1 Ridges
      R21V = LineICP(VnodesT,NCT,VnodesB,NCB,kdt_VBS,VsegmB) #find registered mesh 2 to 1 Valleys
      
      # Determine translation required for best fit. Known that minnimum is at this position [D. Du et al]
      B12R,T12R = NCB[np.array(R12R[:,0],int),],R12R[:,1:4]
      B12V,T12V = NCB[np.array(R12V[:,0],int),],R12V[:,1:4]
      B21R,T21R = R21R[:,1:4],NCT[np.array(R21R[:,0],int),]
      B21V,T21V = R21V[:,1:4],NCT[np.array(R21V[:,0],int),]
      Base_i = np.r_[B12R,B12V,B21R,B21V]
      Target_i = np.r_[T12R,T12V,T21R,T21V]
    else:
      kdt_Targ = KDTree(NCT,20) 
      print 'k-d trees and closest point search'
      B2T,T2B = np.zeros((NCB.shape[0],1)),np.zeros((NCT.shape[0],1))
      print '	base'
      results = pprocess.Map(limit=LIM)
      calc = results.manage(pprocess.MakeParallel(FullICP))
      for j in range(0,LIM):
	calc(np.array(range(0,NPP1))+j*NPP1,NCB,kdt_Targ,B2T)
      for j in range(0,LIM):
	B2T[np.array(range(0,NPP1))+j*NPP1,:] = results[j]
      B2T[range(LIM*NPP1,N1),:]=FullICP(range(LIM*NPP1,N1),NCB,kdt_Targ,B2T)
      #for i in range(0,NCB.shape[0]):
	#B2T = B2T + [kdt_Targ.query(NCB[i,])[1]]
      print '	target'
      results = pprocess.Map(limit=LIM)
      calc = results.manage(pprocess.MakeParallel(FullICP))
      for j in range(0,LIM):
	calc(np.array(range(0,NPP2))+j*NPP2,NCT,kdt_Base,T2B)
      for j in range(0,LIM):
	T2B[np.array(range(0,NPP2))+j*NPP2,:] = results[j]
      T2B[range(LIM*NPP2,N2),:]=FullICP(range(LIM*NPP2,N2),NCT,kdt_Base,T2B)
      #for i in range(0,NCT.shape[0]):
	#T2B = T2B + [kdt_Base.query(NCT[i,])[1]]
      B2T,T2B = np.array(B2T,int).reshape((NCB.shape[0],)),np.array(T2B,int).reshape((NCT.shape[0],))
      Base_i = np.r_[NCB,NCB[T2B,]]
      Target_i = np.r_[NCT[B2T,],NCT]
    print 'translate'
    Distance = Base_i - Target_i
    Translate = np.sum(Distance,0)/Distance.shape[0]
    # Apply translation to Target nodal coordinates and Base to Target registered pairs
    NCTT = np.c_[NCT[:,0]+Translate[0],NCT[:,1]+Translate[1],NCT[:,2]+Translate[2]]
    Target_i = np.c_[Target_i[:,0]+Translate[0],Target_i[:,1]+Translate[1],Target_i[:,2]+Translate[2]]
    print '	translation = ',Translate
    # Determine Current Quadratic approximation:
    E1,E2,E3 = np.matrix([[0,-1,0],[1,0,0],[0,0,0]]),np.matrix([[0,0,0],[0,0,-1],[0,1,0]]),np.matrix([[0,0,1],[0,0,0],[-1,0,0]])
    D1,D2,D3 = np.matrix([[1,0,0],[0,0,0],[0,0,0]]),np.matrix([[0,0,0],[0,1,0],[0,0,0]]),np.matrix([[0,0,0],[0,0,0],[0,0,1]])
    IMat = np.matrix(np.eye(3))         
    #Ukm1,Skm1,Rkm1 = np.linalg.svd(AO)
    Ukm1,Skm1,Rkm1 = IMat,IMat,IMat#*1.3
    #k2,diff2=0,2
    #while (diff2>0.1)&(k2<100):
      #k2=k2+1
      #Bi1,Bi2,Bi3 = np.array(Ukm1*E1*Skm1*Rkm1*Target_i.T).T,np.array(Ukm1*E2*Skm1*Rkm1*Target_i.T).T,np.array(Ukm1*E3*Skm1*Rkm1*Target_i.T).T
      #Bi4,Bi5,Bi6 = np.array(Ukm1*Skm1*D1*Rkm1*Target_i.T).T,np.array(Ukm1*Skm1*D2*Rkm1*Target_i.T).T,np.array(Ukm1*Skm1*D3*Rkm1*Target_i.T).T
      #Bi7,Bi8,Bi9 = np.array(Ukm1*Skm1*Rkm1*E1*Target_i.T).T,np.array(Ukm1*Skm1*Rkm1*E2*Target_i.T).T,np.array(Ukm1*Skm1*Rkm1*E3*Target_i.T).T
      
      #Distance = np.array(Ukm1*Skm1*Rkm1*Target_i.T-Base_i.T).T
      ## Set up Hessian matrix:
      #Hes = np.matrix([[np.sum(Bi1*Bi1),np.sum(Bi1*Bi2),np.sum(Bi1*Bi3),np.sum(Bi1*Bi4),np.sum(Bi1*Bi5),np.sum(Bi1*Bi6),np.sum(Bi1*Bi7),np.sum(Bi1*Bi8),np.sum(Bi1*Bi9)],
	#[np.sum(Bi2*Bi1),np.sum(Bi2*Bi2),np.sum(Bi2*Bi3),np.sum(Bi2*Bi4),np.sum(Bi2*Bi5),np.sum(Bi2*Bi6),np.sum(Bi2*Bi7),np.sum(Bi2*Bi8),np.sum(Bi2*Bi9)],
	#[np.sum(Bi3*Bi1),np.sum(Bi3*Bi2),np.sum(Bi3*Bi3),np.sum(Bi3*Bi4),np.sum(Bi3*Bi5),np.sum(Bi3*Bi6),np.sum(Bi3*Bi7),np.sum(Bi3*Bi8),np.sum(Bi3*Bi9)],
	#[np.sum(Bi4*Bi1),np.sum(Bi4*Bi2),np.sum(Bi4*Bi3),np.sum(Bi4*Bi4),np.sum(Bi4*Bi5),np.sum(Bi4*Bi6),np.sum(Bi4*Bi7),np.sum(Bi4*Bi8),np.sum(Bi4*Bi9)],
	#[np.sum(Bi5*Bi1),np.sum(Bi5*Bi2),np.sum(Bi5*Bi3),np.sum(Bi5*Bi4),np.sum(Bi5*Bi5),np.sum(Bi5*Bi6),np.sum(Bi5*Bi7),np.sum(Bi5*Bi8),np.sum(Bi5*Bi9)],
	#[np.sum(Bi6*Bi1),np.sum(Bi6*Bi2),np.sum(Bi6*Bi3),np.sum(Bi6*Bi4),np.sum(Bi6*Bi5),np.sum(Bi6*Bi6),np.sum(Bi6*Bi7),np.sum(Bi6*Bi8),np.sum(Bi6*Bi9)],
	#[np.sum(Bi7*Bi1),np.sum(Bi7*Bi2),np.sum(Bi7*Bi3),np.sum(Bi7*Bi4),np.sum(Bi7*Bi5),np.sum(Bi7*Bi6),np.sum(Bi7*Bi7),np.sum(Bi7*Bi8),np.sum(Bi7*Bi9)],
	#[np.sum(Bi8*Bi1),np.sum(Bi8*Bi2),np.sum(Bi8*Bi3),np.sum(Bi8*Bi4),np.sum(Bi8*Bi5),np.sum(Bi8*Bi6),np.sum(Bi8*Bi7),np.sum(Bi8*Bi8),np.sum(Bi8*Bi9)],
	#[np.sum(Bi9*Bi1),np.sum(Bi9*Bi2),np.sum(Bi9*Bi3),np.sum(Bi9*Bi4),np.sum(Bi9*Bi5),np.sum(Bi9*Bi6),np.sum(Bi9*Bi7),np.sum(Bi9*Bi8),np.sum(Bi9*Bi9)]])
	
      ## set up fj
      #Fj = np.matrix([[np.sum(Bi1*Distance),np.sum(Bi2*Distance),np.sum(Bi3*Distance),np.sum(Bi4*Distance),np.sum(Bi5*Distance),
	#np.sum(Bi6*Distance),np.sum(Bi7*Distance),np.sum(Bi8*Distance),np.sum(Bi9*Distance)]]).T
    print 'find rotation and scale'
    Xf = fmin_powell(costF,np.array([0,0,0,1,1,1,0,0,0]),args=(Base_i,Target_i,Ukm1,Skm1,Rkm1,UseScale,Use1Scale))
    Xf = np.array(Xf).reshape((9,))
    print 	'u1..3,s1..3,r1..3 = ',Xf
    FvP = np.array(Ukm1*Skm1*Rkm1*Target_i.T-Base_i.T).T
    FvP=np.sum(FvP*FvP)
    Conv[k]=FvP
    Ukm1 = Ukm1+Ukm1*np.matrix(E1*Xf[0]+E2*Xf[1]+E3*Xf[2])
    if UseScale==1:
      Skm1 = Skm1+Skm1*np.matrix(D1*Xf[3]+D2*Xf[4]+D3*Xf[5])
    if Use1Scale==1:
      Skm1 = Skm1+Skm1*np.matrix(D1*Xf[3]+D2*Xf[3]+D3*Xf[3])
    print
    print Skm1
    print
    Rkm1 = Rkm1+Rkm1*np.matrix(E1*Xf[6]+E2*Xf[7]+E3*Xf[8])
    NCTT = np.array(Ukm1*Skm1*Rkm1*NCTT.T).T
    diff = np.sum((NCT-NCTT)*(NCT-NCTT))/NCT.shape[0]
    print '	Average difference between current and previous nodal coordinates:  ',diff
    np.ma.dump(NCTT,'Femur1NC_'+str(k))
    NCT = NCTT
  return NCT,Conv

def FullICP(nodes1,NC1,kdt_2,CPlist):
  for i in nodes1:
    CPlist[i,]=kdt_2.query(NC1[i,])[1]
  return CPlist[nodes1,]
    
    
def LineICP(Nodes1,NC1,Nodes2,NC2,kdt_LS2,LS2_info):
  registeredP = np.zeros((Nodes1.shape[0],6))
  for i in range(0,Nodes1.shape[0]):
    # find closest line segments
    segm = kdt_LS2.query(NC1[Nodes1[i,0],])[1]
    a,b = NC2[LS2_info[segm,0],]-NC2[LS2_info[segm,1],],NC2[LS2_info[segm,1],]-NC1[Nodes1[i,0],]
    alpha = -np.sum(a*b)/np.sum(a*a)
    NCR = alpha*NC2[LS2_info[segm,0],]+(1-alpha)*NC2[LS2_info[segm,1],]
    dist = np.sqrt(np.sum(alpha*alpha*a*a + 2*alpha*a*b + b*b))
    registeredP[i,]=np.r_[Nodes1[i,0],NCR,dist,segm]
  for i in Nodes2[:,0]:	#remove doule registeres final coordinates by keeping only closest register
    pos = np.where((registeredP[:,1]==NC2[i,0])&(registeredP[:,2]==NC2[i,1])&(registeredP[:,3]==NC2[i,2]))[0]
    if pos.size>1:
      keep = np.where(((registeredP[:,1]<>NC2[i,0])|(registeredP[:,2]<>NC2[i,1])|(registeredP[:,3]<>NC2[i,2]))|(registeredP[:,4]==np.min(registeredP[pos,4])))[0]
      registeredP = registeredP[keep,]
  return registeredP	# registeredP in form [node,coordinates after registration,distance] (nx5 matrix)
  
#def LineICPThresh(Reg12,Reg21,Segm1,Segm2,NC1,NC2):

def costF(Xf,Base_i,Target_i,Ukm1,Skm1,Rkm1,UseScale,Use1Scale):
  E1,E2,E3 = np.matrix([[0,-1,0],[1,0,0],[0,0,0]]),np.matrix([[0,0,0],[0,0,-1],[0,1,0]]),np.matrix([[0,0,1],[0,0,0],[-1,0,0]])
  D1,D2,D3 = np.matrix([[1,0,0],[0,0,0],[0,0,0]]),np.matrix([[0,0,0],[0,1,0],[0,0,0]]),np.matrix([[0,0,0],[0,0,0],[0,0,1]])
  Ukm1 = Ukm1+Ukm1*np.matrix(E1*Xf[0]+E2*Xf[1]+E3*Xf[2])
  if UseScale==1:
    Skm1 = Skm1+Skm1*np.matrix(D1*Xf[3]+D2*Xf[4]+D3*Xf[5])
  if Use1Scale==1:
    Skm1 = Skm1+Skm1*np.matrix(D1*Xf[3]+D2*Xf[3]+D3*Xf[3])
  Rkm1 = Rkm1+Rkm1*np.matrix(E1*Xf[6]+E2*Xf[7]+E3*Xf[8])
  Distance = np.array(Ukm1*Skm1*Rkm1*Target_i.T-Base_i.T).T
  #C = np.matrix(XO).T
  #Fc = C.T*Hes*C + 2*Fj.T*C + np.sum(np.sum(Distance*Distance,1))
  Fc = np.sum(Distance*Distance)
  #print Hes
  #print Fc
  #print C.T,
  #print
  return Fc#,Fj#np.array(Hes)
  
  
def BuildMap(NC1,NC2,nodes1,nodes2,segm1,segm2,Reg1,Reg2,percentReg):
  # takes in Nodal coordinates, nodes [nd nr, line nr]; segments [nd1, nd2, line nr] and registered results [nd registered, deformed NC, distance, segment on opposite features]
  # Build map where "allowable register" is if the registered line also registers back
  B2T,T2B = np.zeros((Reg1.shape[0],2)),np.zeros((Reg2.shape[0],2))
  #print Reg1.shape,Reg2.shape
  for i in range(0,Reg1.shape[0]):
    B2T[i,] = np.array([nodes1[np.where(nodes1[:,0]==Reg1[i,0])[0][0],1],segm2[Reg1[i,5],2]])
  for i in range(0,Reg2.shape[0]):
    T2B[i,] = np.array([nodes2[np.where(nodes2[:,0]==Reg2[i,0])[0][0],1],segm1[Reg2[i,5],2]])
  
  B2Ttop,T2Btop = np.zeros((Reg1.shape[0],))*False,np.zeros((Reg2.shape[0],))*False
  
  lnBnr,lnTnr = np.max(nodes1[:,1]),np.max(nodes2[:,1])
  for i in range(0,lnBnr):
    check = np.where(B2T[:,0]==i+1)[0]
    if check.size > 0:
      for j in check:
	OppB = np.where(B2T[check,1]==B2T[j,1])[0].size
	checkOpp = np.where(T2B[:,0]==B2T[j,1])[0]
	if checkOpp.size>0:
	  topCorr = np.where(T2B[checkOpp,1]==B2T[j,0])[0].size
	  if topCorr>0:
	    pT2Bj,pB2Tj = (np.array(topCorr,float)/checkOpp.size)*100>percentReg,(np.array(OppB,float)/check.size)*100>percentReg
	    B2Ttop[j]=pT2Bj|pB2Tj
  for i in range(0,lnTnr):
    check = np.where(T2B[:,0]==i+1)[0]
    if check.size > 0:
      for j in check:
	OppB = np.where(T2B[check,1]==T2B[j,1])[0].size
	checkOpp = np.where(B2T[:,0]==T2B[j,1])[0]
	if checkOpp.size>0:
	  topCorr = np.where(B2T[checkOpp,1]==T2B[j,0])[0].size
	  if topCorr>0:
	    pB2Tj,pT2Bj = (np.array(topCorr,float)/checkOpp.size)*100>percentReg,(np.array(OppB,float)/check.size)*100>percentReg
	    T2Btop[j]=pT2Bj|pB2Tj
	    
  Reg12 = Reg1[np.where(B2Ttop)[0],]
  Reg21 = Reg2[np.where(T2Btop)[0],]
  print '	Nr of registered points:  ',np.where(B2Ttop)[0].size,' & ',np.where(T2Btop)[0].size
  return Reg12,Reg21
  
def GaussMove(nodes,NCB,NC_BF,DISP_BFEAT,DISP_BFull,Supp):
  for i in nodes:
      G1node = cdist(NC_BF,NCB[i,:].reshape((1,3)),'euclidean')
      G1node = np.ones((3,1))*np.exp(-G1node.T/Supp)
      DISP_BFull[i,:] = np.sum(DISP_BFEAT.T*G1node,1)/np.sum(G1node[0,:],0)
  return DISP_BFull[nodes,:]
  
def AllowableNaT(NC,Tconnect,FeatAreaNodes,AllowNodes,Rlines,Vlines,keepLNnodes,doTri=1):
  AllowTri = []
  RemoveNodes = np.array([])
  LineNodes = np.array([])
  print "Get line nodes and set up kd-tree"
  for i in range(1,Rlines[0]+1):
    LineNodes = np.array(np.r_[LineNodes,Rlines[i]],int)
  for i in range(1,Vlines[0]+1):
    LineNodes = np.array(np.r_[LineNodes,Vlines[i]],int)
  kdtLN = KDTree(NC[LineNodes,],5)
  print "Find allowable and remove additional nodes"
  for i in FeatAreaNodes:
    cn = LineNodes[np.array(kdtLN.query(NC[i,])[1],int)]
    if np.where(keepLNnodes==cn)[0].size>0:
      AllowNodes = np.r_[AllowNodes,i]
    else: 
      RemoveNodes = np.r_[RemoveNodes,i]
  AllowNodes,RemoveNodes = np.array(AllowNodes,int),np.array(RemoveNodes,int)
  # remove possible repeated indices from "AllowNodes"
  #remove=find_repeats(AllowNodes)[0]
  if doTri==1:
    print "Update allowable triangles"
    for i in range(Tconnect.shape[0]):
      if find_repeats(np.r_[AllowNodes,Tconnect[i,]])[0].size>0:
	AllowTri.append(i)
    AllowTri = np.array(AllowTri,int)
  return AllowNodes,RemoveNodes,AllowTri
  