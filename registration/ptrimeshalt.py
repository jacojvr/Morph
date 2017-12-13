# -*- coding: utf-8 -*-
# Triangular surface coarsening

#do:
# from ptrimeshalt import *

import numpy as np
import scipy as sp
from scipy.stats.mstats import find_repeats
import pprocess
from scipy.spatial import KDTree

LIM = 8 # max number of parallel processes

def NodeselectFP(NC,Tri,Neig,vv):
  #spheres&saddles:
  spsad=np.where((Neig[:,0]<1.5*Neig[:,1])&(Neig[:,1]<1.5*Neig[:,2])&(Neig[:,2]>0))[0]
  rest=np.where((Neig[:,0]>=1.5*Neig[:,1])|(Neig[:,1]>=1.5*Neig[:,2])|(Neig[:,2]<=0))[0]
  print "Nr of Spheres & Sadles: ",spsad.size
  #ridges and valleys
  ridgeval1=np.where((Neig[rest,0]<1.5*Neig[rest,1])&(Neig[rest,1]>10*Neig[rest,2])&(np.abs(Neig[rest,2])<0.001))[0]
  ridgeval=rest[ridgeval1,]
  rest1=np.where((Neig[rest,0]>=1.5*Neig[rest,1])|(Neig[rest,1]<=10*Neig[rest,2])|(np.abs(Neig[rest,2])>=0.001))[0]
  rest=rest[rest1,]
  print "Nr of Ridges & Valleys: ",ridgeval.size
  #planes:
  plane1=np.where((Neig[rest,0]>10*Neig[rest,1])&(np.abs(Neig[rest,1])<0.001)&(np.abs(Neig[rest,2])<0.001))[0]
  plane=rest[plane1,]
  rest1=np.where((Neig[rest,0]<=10*Neig[rest,1])|(np.abs(Neig[rest,1])>=0.001)|(np.abs(Neig[rest,2])>=0.001))[0]
  rest=rest[rest1,]
  print "Nr of Planes: ",plane.size
  #keep all for eig1 < 3*eig2
  keep1=np.where(Neig[rest,0]<5*Neig[rest,1])[0]
  keep = np.r_[spsad,ridgeval,keep1]	# take all features&sharp edges except planes
  rest=np.where(Neig[rest,0]>=5*Neig[rest,1])[0]
  track=0
  ##keep every 2nd for eig1<20*eig2 and every 3rd for eig<500*eig2:
  #Neig=np.ma.load('EigvalDolph1n.txt')
  ##vv=np.array([50,100,1000])
  for inc in vv:
    track=track+1
    print "Step ",track," of ",vv.size
    rows2=np.where(Neig[rest,0]<inc*Neig[rest,1])[0]
    rest1=np.where(Neig[rest,0]>=inc*Neig[rest,1])[0]
    rows2=rest[rows2,]
    rest=rest[rest1,]
    if rows2.size > 0:
      cont=1
      while cont>0:
	node = rows2[0,]
	keep=np.r_[keep,node]
	neigh=Nneighbours(node,NC,Tri,track)[1]
	torem=find_repeats(np.r_[rows2,neigh])[0]
	for i in np.r_[node,torem]:
	  rows2 = rows2[np.r_[np.where(rows2<i)[0],np.where(rows2>i)[0]],]
	#print 'possible selections in current category: ',rows2.size
	if rows2.size==0:
	  cont=0
      
      
    #for i in rows2:
      #print "node: ",i
      #neigh=Nneighbours(i,NC,Tri,track)
      #if find_repeats(np.r_[keep,neigh])[0].size==0:
	#keep=np.r_[keep,i]
	#print "	keep"
  keep=np.array(keep,int)
  keep.sort()
  return keep


def Neigenvalues(NC,Tri,layer,ls):
  KDTNC = KDTree(NC,ls)
  VNORM = vertexnormal(NC,Tri)
  Neig = np.zeros(NC.shape)
  N1 = Neig.size/3;
  NPP1 = N1/LIM
  results = pprocess.Map(limit=LIM)
  calc = results.manage(pprocess.MakeParallel(NodeEigVal))
  for j in range(0,LIM):
      calc(np.array(range(0,NPP1))+j*NPP1,NC,Tri,KDTNC,ls,VNORM,layer,Neig)
  for j in range(0,LIM):
      Neig[np.array(range(0,NPP1))+j*NPP1,] = results[j]
  Neig[range(LIM*NPP1,N1),]=NodeEigVal(range(LIM*NPP1,N1),NC,Tri,KDTNC,ls,VNORM,layer,Neig)
  #queue=pprocess.Queue(limit=LIM)
  #results = []
  #Nneig = queue.manage(pprocess.MakeParallel(NodeEigVal))
  #for i in range(0,N):
    #Nneig(i,NC,Tri,KDTNC,ls,VNORM,2)
  #for i,NE in queue:
    #Neig[i,] = NE
  rows=np.where(Neig[:,0]<Neig[:,1])[0]
  Neig[rows,:]=np.c_[Neig[rows,1],Neig[rows,0],Neig[rows,2]]
  rows=np.where(Neig[:,0]<Neig[:,2])[0]
  Neig[rows,:]=np.c_[Neig[rows,2],Neig[rows,1],Neig[rows,0]]
  rows=np.where(Neig[:,1]<Neig[:,2])[0]
  Neig[rows,:]=np.c_[Neig[rows,0],Neig[rows,2],Neig[rows,1]]
  return Neig
 
def NodeEigVal(nodes,NC,Tri,KDTNC,ls,VNORM,layer,Neig):
  for node in nodes:
    #neigh=Nneighbours(node,NC,Tri,layer)[1]
    neigh = KDTNC.query(NC[node,],ls)[1]
    #neigh = neigh[np.where(np.dot(VNORM[neigh,],VNORM[node,].reshape((3,)))>-0.5)[0]]
    Tv=localstruct(VNORM[neigh,])
    Neig[node,] = np.linalg.eig(Tv)[0]
  return Neig[nodes,]

def Nneighbours(node,NC,Tri,layer):
  # determine neighbouring nodes of a specific node up to a set of "layers"
  # inputs:
  # node: specific node to query neighbours, NC&Tri: nodal coord and connectivity
  # layer: up to how many degrees of separation a neighbour shoud be searched
  sz = Tri.shape[1]
  aa=np.where(Tri==node)[0]
  pn=Tri[aa,].reshape((aa.size*sz,))
  neigh=find_repeats(pn)[0]
  neighOne = np.array(neigh,int)
  neighOne = neighOne[np.where(neighOne<>node)[0]]
  if layer>1:
    for l in range(1,layer):
      for nn in neigh:
	aa = np.where(Tri==nn)[0]
	pn=Tri[aa,].reshape((aa.size*sz,))
	pn=find_repeats(pn)[0]
	for i in pn:
	  if find_repeats(np.r_[neigh,i])[0].size==0:
	      neigh = np.r_[neigh,i]
  neigh=np.array(neigh[np.where(neigh<>node)[0],],int)
  #if layer==1:	#if only actual neighbours are calculated, returmn then in right hand rule order for retriangulation purposes
    #neighC=np.empty((neigh.size,))
    #neighC[0,]=neigh[0,]
    #rp,cp=np.where(Tri[neighpos,]==neigh[0,])
    #for i in range(1,neigh.size):
      #r,c=np.where((Tri[neighpos[rp,]]<>neighC[i-1,])&(Tri[neighpos[rp,]]<>node))
      #neighC[i,]=max(((cp+1==c)|((cp==2)*(c==0)))*Tri[neighpos[rp,],c])
      #rp,cp=np.where(Tri[neighpos,]==neighC[i,])
    #neigh = np.array(neighC,int)
  return neighOne,neigh
  
def TriEdge(node,NC,Tri):
  # determine edge returned as neighbouring node and triangles on either side of edges
  triangles = []
  aa=np.where(Tri==node)[0]
  pn=Tri[aa,].reshape((aa.size*3,))
  neigh=find_repeats(pn)[0]
  neighOne = np.array(neigh,int)
  neighOne = neighOne[np.where(neighOne<>node)[0]]
  for i in neighOne:
    tr = np.where(Tri[aa,]==i)[0]
    triangles = triangles+[[aa[tr[0]],aa[tr[1]]]]
  triangles = np.array(triangles,int)
  return neighOne,triangles
  
  

def localstruct(Neighnorm):
  # determine the local coordinate structure for a specific vertex:
  # inputs:
  # Neighnorm = vertex normal vector based on incident triangle normals for neighbouring triangles
  nx = Neighnorm[:,0]
  ny = Neighnorm[:,1]
  nz = Neighnorm[:,2]
  nxx = np.sum(nx*nx)
  nxy = np.sum(nx*ny)
  nxz = np.sum(nx*nz)
  nyy = np.sum(ny*ny)
  nyz = np.sum(ny*nz)
  nzz = np.sum(nz*nz)
  #LST = np.array([[nxx,nxy,nxz],[nxy,nyy,nyz],[nxz,nyz,nzz]])
  return np.array([[nxx,nxy,nxz],[nxy,nyy,nyz],[nxz,nyz,nzz]])

def vertexnormal(NC,Tri):
  # Program takes in the nodal coordinates (NC=[x,y,z]) of all vertices as well as the triangle connectivity (Tri=[n1,n2,n3]),
  # to determine the normals of each vertex (to do curvature estimation)
  # first, determine triangle centroids:
  Centr = np.c_[np.sum(np.c_[NC[Tri[:,0],0],NC[Tri[:,1],0],NC[Tri[:,2],0]],1)/3,
      np.sum(np.c_[NC[Tri[:,0],1],NC[Tri[:,1],1],NC[Tri[:,2],1]],1)/3,
      np.sum(np.c_[NC[Tri[:,0],2],NC[Tri[:,1],2],NC[Tri[:,2],2]],1)/3]
  # determine triangle normals
  T_normvec = np.cross(NC[Tri[:,1],:]-NC[Tri[:,0],:],NC[Tri[:,2],:]-NC[Tri[:,0],:])
  T_normvec = (T_normvec.T/(np.ones((3,1))*np.sqrt(np.sum(np.array([T_normvec*T_normvec]),2)))).T
  
  # number of triangles and nodes:
  nrT = np.array(Tri.shape)[0,]
  nrN = np.array(NC.shape)[0,]
  # inintialize vertex normals:
  Nnorm = np.zeros(NC.shape)
  
  # determine direction of vertex normal through contribution of each triangle (weight=1/distance from centroid):
  for i in range(0,nrT):
    weight = NC[Tri[i,],]-Centr[i,]
    weight = 1/np.sqrt(np.sum(weight*weight,1))
    Nnorm[Tri[i,],]=Nnorm[Tri[i,],]+np.ones((3,1))*T_normvec[i,]*(np.ones((3,1))*weight).T

  Nnorm = (Nnorm.T/(np.ones((3,1))*np.sqrt(np.sum(np.array([Nnorm*Nnorm]),2)))).T
  return Nnorm
  
def vrtxnormal(NC,Tri,Centr,T_normvec):
  # Program takes in the nodal coordinates (NC=[x,y,z]) of all vertices as well as the triangle connectivity (Tri=[n1,n2,n3]),
  # to determine the normals of each vertex (to do curvature estimation)  
  # number of triangles and nodes:
  nrT = np.array(Tri.shape)[0,]
  nrN = np.array(NC.shape)[0,]
  # inintialize vertex normals:
  Nnorm = np.zeros(NC.shape)
  
  # determine direction of vertex normal through contribution of each triangle (weight=1/distance from centroid):
  for i in range(0,nrT):
    weight = NC[Tri[i,],]-Centr[i,]
    weight = 1/np.sqrt(np.sum(weight*weight,1))
    Nnorm[Tri[i,],]=Nnorm[Tri[i,],]+np.ones((3,1))*T_normvec[i,]*(np.ones((3,1))*weight).T

  Nnorm = (Nnorm.T/(np.ones((3,1))*np.sqrt(np.sum(np.array([Nnorm*Nnorm]),2)))).T
  return Nnorm