# -*- coding: utf-8 -*-

import numpy as np
import scipy as sp
from scipy.stats.mstats import find_repeats
import pprocess
from ptrimeshalt import *
from scipy.spatial import KDTree

LIM=8

def crestlines(nodes,NC,Tri,NN=3,ls=100,neighSpat=0,eOrk=0):
  # neighSpat == 1: Use mesh connectivity and NN-ring neighbours
  # otherwise only use spatial search with ls closest neighbours
  # eOrk == 1: use zero-crossing to determine crest nodes
  # default neighSpat=0,eOrk=0
  NCQ = NC.shape[0]
  N1 = nodes.size
  NPP1 = N1/LIM
  T1 = Tri.shape[0]
  TPP1 = T1/LIM
  # Determine Triangle unit normal
  TBC = np.c_[np.sum(np.c_[NC[Tri[:,0],0],NC[Tri[:,1],0],NC[Tri[:,2],0]],1)/3,
    np.sum(np.c_[NC[Tri[:,0],1],NC[Tri[:,1],1],NC[Tri[:,2],1]],1)/3,np.sum(np.c_[NC[Tri[:,0],2],NC[Tri[:,1],2],NC[Tri[:,2],2]],1)/3]
  TNORM = np.cross(NC[Tri[:,1],:]-NC[Tri[:,0],:],NC[Tri[:,2],:]-NC[Tri[:,0],:])
  TNORM = (TNORM.T/(np.ones((3,1))*np.sqrt(np.sum(np.array([TNORM*TNORM]),2)))).T
  VNORM = vrtxnormal(NC,Tri,TBC,TNORM)
  
  # get k-d tree for NC
  KDTNC = KDTree(NC,5)
  
  print 'Get principal curvature, -direction and -derivatives'
  # Determine the principal curvature directions at each node
  Kmax,Kmin,eMax,eMin,PdMax,PdMin = np.zeros((NCQ,1)),np.zeros((NCQ,1)),np.zeros((NCQ,1)),np.zeros((NCQ,1)),np.zeros((NCQ,3)),np.zeros((NCQ,3))
  #KmaxMCDval,KMCD =np.zeros((N1,10)),np.zeros((N1,10))
  results = pprocess.Map(limit=LIM)
  calc = results.manage(pprocess.MakeParallel(MLSmethod))
  for j in range(0,LIM):
      calc(nodes[np.array(range(0,NPP1))+j*NPP1],NC,Tri,KDTNC,ls,VNORM,Kmax,Kmin,eMax,eMin,PdMax,PdMin,NN,neighSpat)
  for j in range(0,LIM):
      [Kmax[nodes[np.array(range(0,NPP1))+j*NPP1],:],Kmin[nodes[np.array(range(0,NPP1))+j*NPP1],:],
      eMax[nodes[np.array(range(0,NPP1))+j*NPP1],:],eMin[nodes[np.array(range(0,NPP1))+j*NPP1],:],
      PdMax[nodes[np.array(range(0,NPP1))+j*NPP1],:],PdMin[nodes[np.array(range(0,NPP1))+j*NPP1],:]] = results[j]
  [Kmax[nodes[range(LIM*NPP1,N1)],:],Kmin[nodes[range(LIM*NPP1,N1)],:],
  eMax[nodes[range(LIM*NPP1,N1)],:],eMin[nodes[range(LIM*NPP1,N1)],:],PdMax[nodes[range(LIM*NPP1,N1)],:],
  PdMin[nodes[range(LIM*NPP1,N1)],:]] = MLSmethod(nodes[range(LIM*NPP1,N1)],NC,Tri,KDTNC,ls,VNORM,Kmax,Kmin,eMax,eMin,PdMax,PdMin,NN,neighSpat)
  #MLSmethod(nodes[range(0,N1)],NC,Tri,KDTNC,ls,VNORM,Kmax,Kmin,eMax,eMin,PdMax,PdMin,NN)
  #PdMax = (PdMax.T/(np.ones((3,1))*np.sqrt(np.sum(np.array([PdMax*PdMax]),2)))).T
  #PdMin = (PdMin.T/(np.ones((3,1))*np.sqrt(np.sum(np.array([PdMin*PdMin]),2)))).T
  
  print 'Determine Ridge nodes'
  ridgeN = ridgenodes(NC,Tri,PdMax,PdMin,Kmax,Kmin,eMax,eMin,1,eOrk)
  print 'Determine Valley nodes'
  valleyN = ridgenodes(NC,Tri,PdMax,PdMin,Kmax,Kmin,eMax,eMin,0,eOrk)
  
  
  return ridgeN,valleyN,Kmax,Kmin,eMax,eMin,PdMax,PdMin


def MLSmethod(nodes,NC,Tri,KDTNC,ls,VNORM,Kmax,Kmin,eMax,eMin,PdMax,PdMin,NN,neighSpat):
  #Kmax,Kmin,eMax,eMin,PdMax,PdMin = KM[:,0],KM[:,1],KM[:,2],KM[:,3],KM[:,4:7],KM[:,7:10]
  for ni in nodes:
    # Determine principal curvature direction for node "ni" on unstructured surface mesh.
    if neighSpat==1:
      'Use connectivity'
      neighb1,neighb = Nneighbours(ni,NC,Tri,NN)
      DistV = NC[neighb,:]-np.ones((neighb.size,1))*NC[ni,:]
      DistS = np.sqrt(np.sum(DistV*DistV,1))
    else:
      'Use spatial distribution'
      DistS,neighb = KDTNC.query(NC[ni,],ls)
    #h=np.average(DistS)
    ######################
    ######################   REMOVE NEIGHBOURS WHOSE NORMALS MAKE OBTUSE ANGLES WITH CURRENT NODE ##############
    #neighb1 = neighb1[np.where(np.dot(VNORM[neighb1,],VNORM[ni,].reshape((3,)))>0)[0]]
    keep = np.where(np.dot(VNORM[neighb,],VNORM[ni,].reshape((3,)))>0.2)[0]
    neighb,DistS = neighb[keep],DistS[keep]
    if neighSpat<>1:
      neighb1 = neighb[np.where(DistS<2*DistS[1])[0]]
    #h1 = NC[neighb1,:]-np.ones((neighb1.size,1))*NC[ni,:]
    #h1 = np.sqrt(np.sum(h1*h1,1))
    #h = np.sum(h1)/(neighb1.size)
    #Gw1 = np.exp(-np.power(h1,2)/(h*h)) # Gaussian weight
    PlaneN = np.sum(VNORM[neighb1,:],0)
    PlaneN = PlaneN/(np.ones((3,))*np.sqrt(np.sum(PlaneN*PlaneN)))
    keep = np.where(np.dot(VNORM[neighb,],PlaneN)>0)[0]
    neighb,DistS = neighb[keep],DistS[keep]
    DistV = NC[neighb,:]-np.ones((neighb.size,1))*NC[ni,:]
    #DistS = np.sqrt(np.sum(DistV*DistV,1))
     # average length of the 1-neighbourhood edges
    #Gw = np.exp(-np.power(DistS,2)/(h*h)) # Gaussian weight
    fi = np.dot(DistV,PlaneN)#np.sum((np.ones((neighb.size,1))*PlaneN)*DistV,1) # vertical distance from plane
    #fiV = np.c_[fi,fi,fi]*PlaneN
    #fiW = fi*Gw	# Weighted function values for moving least squares fitting
    # rotate coordinate axis to get local frame where n = z':
    rev=1
    if PlaneN[2]<0:
      PlaneN = -PlaneN
      rev=-1
    if PlaneN[2]==1:
      u=np.array([1,0,0])
    else:
      ax = np.array([-PlaneN[1],PlaneN[0],0])
      ax = ax/np.sqrt(np.sum(ax*ax))
      x,y,xy=ax[0],ax[1],ax[0]*ax[1]
      Theta = np.arccos(PlaneN[2])
      ct,st,omct = np.cos(Theta),np.sin(Theta),1-np.cos(Theta)
      u = np.array([omct*x*x+ct,omct*x*y,-st*y])
    PlaneN,u = rev*PlaneN,rev*u
    v = np.cross(PlaneN,u)
    xi = np.dot(DistV,u)
    yi = np.dot(DistV,v)
    # set up matrix for least squares approximation to C^3 polynomial coefficients
    x2,xy,y2=xi*xi,xi*yi,yi*yi
    Aw = np.matrix(np.c_[x2,2*xy,y2,xi*x2/3,xi*xy,yi*xy,yi*y2/3]/2)#*(np.ones((7,1))*Gw).T/2)
    coef = np.linalg.pinv(Aw.T*Aw)*Aw.T*np.matrix(fi).T
    # determine principal curvatures, principal curvature directions and their derivatives
    [k1,k2],t=np.linalg.eig(np.array([[coef[0,0],coef[1,0]],[coef[1,0],coef[2,0]]]))
    #k1,k2=np.abs(k1),np.abs(k2)
    t1 = t[:,0]
    t2 = t[:,1]
    t1M,t2M,eM = np.c_[t1,t1],np.c_[t2,t2],np.array([[coef[3,0],coef[4,0]],[coef[5,0],coef[6,0]]])
    e1,e2=np.sum(t1M.T*t1M*t1M*eM),np.sum(t2M.T*t2M*t2M*eM)
    Pd1=t1[0]*u + t1[1]*v
    Pd2=t2[0]*u + t2[1]*v
    #reverse=0
    if k1>k2:
      Kmax[ni,] = k1
      Kmin[ni,] = k2
      PdMax[ni,] = Pd1
      ts = np.sum(np.cross(PlaneN,Pd1)-Pd2)<1e-4
      PdMin[ni,] = ts*Pd2+(ts-1)*Pd2
      eMax[ni,] = e1
      eMin[ni,] = ts*e2+(ts-1)*e2
    else:
      Kmax[ni,] = k2
      Kmin[ni,] = k1
      PdMax[ni,] = Pd2
      ts = np.sum(np.cross(PlaneN,Pd2)-Pd1)<1e-4
      PdMin[ni,] = ts*Pd1+(ts-1)*Pd1
      eMax[ni,] = e2
      eMin[ni,] = ts*e1+(ts-1)*e1
  return Kmax[nodes,:],Kmin[nodes,:],eMax[nodes,:],eMin[nodes,:],PdMax[nodes,:],PdMin[nodes,:]
  

def ridgenodes(NC,Tri,PdMax,PdMin,Kmax,Kmin,eMax,eMin,RidgeVal,eOrk=0):
  #eOrk = 0	# 1=use eVal zero crossing otherwise use Kval maximum
  if RidgeVal==1: # search ridge nodes
    kReq = (Kmax>np.abs(Kmin))
    eVal = eMax
    kVal = Kmax
    Pd1 = PdMin
    Pd2 = PdMax
  else: # search valley nodes
    kReq = (Kmin<-np.abs(Kmax))
    eVal = eMin
    kVal = -Kmin
    Pd1 = PdMax
    Pd2 = PdMin
  ridgeN = []
  ridgeN = np.array(ridgeN,int)
  pos = np.where(kReq)[0]
  N1 = pos.size
  NPP1 = N1/LIM
  results = pprocess.Map(limit=LIM)
  calc = results.manage(pprocess.MakeParallel(findn))
  for j in range(0,LIM):
      calc(pos[np.array(range(0,NPP1))+j*NPP1],pos,NC,Tri,kVal,eVal,Pd1,Pd2,eOrk)
  for j in range(0,LIM):
      ridgeN = np.r_[ridgeN,results[j]]
  ridgeN = np.r_[ridgeN,findn(pos[range(LIM*NPP1,N1)],pos,NC,Tri,kVal,eVal,Pd1,Pd2,eOrk)]
  return ridgeN
  
def findn(nodes,pos,NC,Tri,kVal,eVal,Pd1,Pd2,eOrk):
  ridgeN =[]
  for ri in nodes:
    neighb = Nneighbours(ri,NC,Tri,1)[1]
    neighb = np.array(find_repeats(np.r_[pos,neighb])[0],int)
    if neighb.size>0:
      Ndirect = NC[neighb,]-np.ones((neighb.size,1))*NC[ri,]
      neighbkReq = neighb[np.where(np.abs(np.dot(Ndirect,Pd1[ri,]))<0.7*np.abs(np.dot(Ndirect,Pd2[ri,])))[0]]
      if neighbkReq.size>0:
	if eOrk==1:
	  eN = eVal[neighbkReq,]
	  reverse = np.where(np.dot(Pd1[neighbkReq,],Pd1[ri,].reshape((3,)))<0)[0]
	  eN[reverse,] = -eN[reverse,]
	  crestp = np.where(eN*eVal[ri]<0)[0]
	  if (crestp.size>0)&(np.prod(kVal[neighbkReq[crestp,]]>kVal[ri])==0):
	    ridgeN = ridgeN + [ri]
	else:
	  if np.prod(kVal[neighbkReq]<kVal[ri])==1:
	    ridgeN = ridgeN + [ri]	# Current node is a ridge node if all above conditions hold
	  #elif np.where(np.abs(np.dot(Ndirect,Pd1[ri,]))<2*np.abs(np.dot(Ndirect,Pd2[ri,])))[0].size==0:
	    #ridgeN = ridgeN + [ri]
      else:
	ridgeN = ridgeN + [ri]
  ridgeN = np.array(ridgeN,int)
  return ridgeN
	

def ThreshLines(lines,Thresh,NC,Kval):
  Tline = []
  noL = 0
  for ln in range(1,lines[0]+1):
    lineT = lines[ln][0]
    #lineP,lineT = lineO,lineO
    #for place in range(0,lines[ln][0].size):
      #for other in range(1,lines[0]+1):
	#if (other<>ln)&(lineO.size<lines[other].size):
	  #pos = np.where(lines[other][0]==lines[ln][0][place])[0]
	  #if pos.size>1:
	    #pos = pos[0]
	  #if (pos.size>0)&(pos>=place):
	    #lineP = lines[ln][0][0:place]
	    #if lineP.size<=lineT.size:
	      #lineT = lineP
    if (lineThresh(np.array(lineT,int),NC,Kval)>Thresh):#&(lineT.size>3):
      Tline = Tline+[np.array(lineT,int)]
      noL = noL+1
  Tline = [noL]+Tline
  return Tline
    
def lineThresh(line,NC,Kval):
  n = line.size-1
  segL2 = NC[line[range(1,n+1)],]-NC[line[range(0,n)],]
  segL2 = np.sqrt(np.sum(segL2*segL2,1))
  kA = Kval[line[range(1,n+1)],]+Kval[line[range(0,n)],]
  Th = np.sum(segL2*kA)/(2*n)
  return Th
  

def LineDraw(ridgeN,NC,Tri,PdC,Pd,kVal):
  lines = []
  nol=0
  drawY=0
  while drawY<ridgeN.size:
    drawY = drawY+1
    node = ridgeN[np.where(kVal[ridgeN]==np.max(kVal[ridgeN]))[0],]
    linePrev = []	# for all nodes that precede node
    lineFol = [node]	# for node and all nodes that succeed node
    ridgeN = ridgeN[np.where(ridgeN<>node)[0]]	# remove "next" form list of available nodes
    neig,neigh2=Nneighbours(node,NC,Tri,2)
    # find 1 and 2 ring neighbours in "ridgeN" & vector to neighbours
    neigh1 = np.array(find_repeats(np.r_[ridgeN,neig])[0],int)
    neigh2 = np.array(find_repeats(np.r_[ridgeN,neigh2])[0],int)
    neigh1P,neigh2P = neigh1,neigh2
    neigN,neigP = neig,neig
    if neigh2.size>0:
      print '######### node ',node,' ##########'
      nol = nol+1
      drawY = 0
      linePrev = []	# for all nodes that precede node
      lineFol = [node[0]]	# for node and all nodes that succeed node
      NdP1P,NdN1N = NC[neigh1,]-np.ones((neigh1.size,1))*NC[node,],NC[neigh1,]-np.ones((neigh1.size,1))*NC[node,]
      NdP2,NdN2 = NC[neigh2,]-np.ones((neigh2.size,1))*NC[node,],NC[neigh2,]-np.ones((neigh2.size,1))*NC[node,]
      NdP1,NdN1 = NC[neig,]-np.ones((neig.size,1))*NC[node,],NC[neig,]-np.ones((neig.size,1))*NC[node,]
      # check allowable neighbours in direction of principal connecting curvature:
      neighN = neigh1[np.where(np.dot(NdN1N,PdC[node].reshape((3,)))>0)[0]]
      neighP = neigh1[np.where(np.dot(NdP1P,PdC[node].reshape((3,)))<0)[0]]
      neighN2 = neigh2[np.where(np.dot(NdN2,PdC[node].reshape((3,)))>0)[0]]
      neighP2 = neigh2[np.where(np.dot(NdP2,PdC[node].reshape((3,)))<0)[0]]
      searchD1 = PdC[node,].reshape((3,))
      searchD2 = -PdC[node,].reshape((3,))
      while neighN2.size>0:	# connect forward
	if neighN.size==1:
	  next = neighN
	#elif neighN2.size==1:
	  #next = neighN2
	elif neighN.size>1:
	  next = neigh1[np.where(np.dot(NdN1N,searchD1)==np.max(np.dot(NdN1N,searchD1)))[0]]
	else:
	  next = neigN[np.where(np.dot(NdN1,searchD1)==np.max(np.dot(NdN1,searchD1)))[0]]
	  lineFol = lineFol+[next[0]]
	  next = neigh2[np.where(np.dot(NdN2,searchD1)==np.max(np.dot(NdN2,searchD1)))[0]]
	for i in range(0,neigh1.size):
	  ridgeN = ridgeN[np.where(ridgeN<>neigh1[i])[0]]
	next = next[0]
	ridgeN = ridgeN[np.where(ridgeN<>next)[0]]
	dY = np.sum(PdC[next,]*searchD1)>0
	searchD1 = (dY*PdC[next,]-(1-dY)*PdC[next,]).reshape((3,))
	Pd2Next = Pd[next,].reshape((3,))
	lineFol = lineFol+[next]
	neigN,neigh2=Nneighbours(next,NC,Tri,2)
	neigh1 = np.array(find_repeats(np.r_[ridgeN,neigN])[0],int)
	neigh2 = np.array(find_repeats(np.r_[ridgeN,neigh2])[0],int)
	NdN1 = NC[neigN,]-np.ones((neigN.size,1))*NC[next,]
	NdN1N = NC[neigh1,]-np.ones((neigh1.size,1))*NC[next,]
	NdN2 = NC[neigh2,]-np.ones((neigh2.size,1))*NC[next,]
	neighN = neigh1[np.where(np.dot(NdN1N,searchD1)>2*np.abs(np.dot(NdN1N,Pd2Next)))[0]]
	neighN2 = neigh2[np.where(np.dot(NdN2,searchD1)>2*np.abs(np.dot(NdN2,Pd2Next)))[0]]
      while neighP2.size>0:	# connect backward
	if neighP.size==1:
	  prev = neighP
	#elif neighP2.size==1:
	  #prev = neighP2
	elif neighP.size>1:
	  prev = neigh1P[np.where(np.dot(NdP1P,searchD2)==np.max(np.dot(NdP1P,searchD2)))[0]]
	else:
	  prev = neigP[np.where(np.dot(NdP1,searchD2)==np.max(np.dot(NdP1,searchD2)))[0]]
	  linePrev = [prev[0]]+linePrev
	  prev = neigh2P[np.where(np.dot(NdP2,searchD2)==np.max(np.dot(NdP2,searchD2)))[0]]
	for i in range(0,neigh1P.size):
	  ridgeN = ridgeN[np.where(ridgeN<>neigh1P[i])[0]]
	prev = prev[0]
	ridgeN = ridgeN[np.where(ridgeN<>prev)[0]]
	dY = np.sum(PdC[prev,]*searchD2)>0
	searchD2 = (dY*PdC[prev,]-(1-dY)*PdC[prev,]).reshape((3,))
	Pd2Prev = Pd[prev,].reshape((3,1))
	linePrev = [prev]+linePrev
	neigP,neigh2=Nneighbours(prev,NC,Tri,2)
	neigh1P = np.array(find_repeats(np.r_[ridgeN,neigP])[0],int)
	neigh2P = np.array(find_repeats(np.r_[ridgeN,neigh2])[0],int)
	NdP1 = NC[neigP,]-np.ones((neigP.size,1))*NC[prev,]
	NdP1P = NC[neigh1P,]-np.ones((neigh1P.size,1))*NC[prev,]
	NdP2 = NC[neigh2P,]-np.ones((neigh2P.size,1))*NC[prev,]
	neighP = neigh1P[np.where(np.dot(NdP1P,searchD2)>2*np.abs(np.dot(NdP1P,Pd2Prev)))[0]]
	neighP2 = neigh2P[np.where(np.dot(NdP2,searchD2)>2*np.abs(np.dot(NdP2,Pd2Prev)))[0]]
      line =  np.array([linePrev + lineFol])
      lines = lines + [line]
      print line
      print '	crest nodes left = ',ridgeN.size
  lines = [nol] + lines
  return lines
  
def LaplacianSmooth(nodes,NCL,Tri):
  for i in nodes:
    adj=Nneighbours(i,NCL,Tri,1)[1]
    NCL[i,:] = np.sum(NCL[adj,],0)/(adj.size)
  return NCL[nodes,:]
  
def ExtrSmooth(NC,Tri,eMax,eMin,kMax,kMin,PdMax,PdMin,it):
  eMaxT,eMinT = eMax,eMin
  RP = np.where(kMax>np.abs(kMin))[0]
  VP = np.where(kMin<-np.abs(kMax))[0]
  for inc in range(0,it):
    print 'iteration ',inc+1
    for i in range(0,NC.shape[0]):
      neigh = Nneighbours(i,NC,Tri,1)[0]
      #neigh = neigh[np.where(np.abs(np.dot(PdMax[neigh,],PdMax[i,]))>np.abs(np.dot(PdMin[neigh,],PdMax[i,])))[0]]
      if np.where(RP==i)[0].size>0:
	neigh = find_repeats(np.r_[neigh,RP])[0]
      else:
	neigh = find_repeats(np.r_[neigh,VP])[0]
      neigh = np.array(neigh,int)
      eMaN,eMiN = eMax[neigh,],eMin[neigh,]
      reverse = np.where(np.dot(PdMin[neigh,],PdMin[i,].reshape((3,)))<0)[0]
      if reverse.size>0:
	eMaN[reverse,],eMiN[reverse,] = -eMaN[reverse,],-eMiN[reverse,]
      eMaxT[i,],eMinT[i,] = (np.sum(eMaN)+eMaxT[i,])/(neigh.size+1),(np.sum(eMiN)+eMinT[i,])/(neigh.size+1)
  return eMaxT,eMinT