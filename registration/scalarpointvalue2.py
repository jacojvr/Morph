# -*- coding: utf-8 -*-
from scipy.spatial.distance import cdist
import numpy as np
import pprocess
LIM = 8

def TetCellToVertex(NC,Connect,ScalarValue):
  N1 = NC.shape[0]
  NPP1 = N1/LIM
  ConnectNC = np.c_[np.sum(np.c_[NC[Connect[:,0],0],NC[Connect[:,1],0],NC[Connect[:,2],0],NC[Connect[:,3],0]],1)/4,
    np.sum(np.c_[NC[Connect[:,0],1],NC[Connect[:,1],1],NC[Connect[:,2],1],NC[Connect[:,3],1]],1)/4,
    np.sum(np.c_[NC[Connect[:,0],2],NC[Connect[:,1],2],NC[Connect[:,2],2],NC[Connect[:,3],2]],1)/4]
  Scalars = np.zeros(NC.shape[0])
  results = pprocess.Map(limit=LIM)
  calc = results.manage(pprocess.MakeParallel(weightSCVAL))
  for j in range(LIM):
    calc(np.array(range(0,NPP1))+j*NPP1,NC,Connect,ConnectNC,ScalarValue)
  for j in range(LIM):
    Scalars[np.array(range(0,NPP1))+j*NPP1] = results[j]
  #for j in range(LIM):
  if N1>LIM*NPP1:
    Scalars[range(LIM*NPP1,N1)] = weightSCVAL(np.array(range(LIM*NPP1,N1)),NC,Connect,ConnectNC,ScalarValue)
  return Scalars
  
def weightSCVAL(nodes,NC,Connect,ConnectNC,ScalarValue):
  SCVAL = np.zeros((nodes.size,))
  for i in range(nodes.size):
    tets = np.where(Connect==nodes[i])[0]
    Weights = cdist(ConnectNC[tets,],NC[nodes[i],].reshape((1,3)),'euclidean')
    SCVAL[i]= np.sum(Weights*ScalarValue[tets])/np.sum(Weights)
  return SCVAL

def TetListScalars(NodeTetList,ScalarVal,AoM = 0):
  # if AoM ==0, use average, if AoM<>0, use the minimum associated with that point
  N1 = len(NodeTetList)
  NPP1 = N1/LIM
  NodeScal = np.zeros((N1,))
  results = pprocess.Map(limit=LIM)
  calc = results.manage(pprocess.MakeParallel(TetListScalarsSec))
  for j in range(LIM):
    calc(np.array(range(0,NPP1))+j*NPP1,NodeTetList,ScalarVal,AoM)
  for j in range(LIM):
    NodeScal[np.array(range(0,NPP1))+j*NPP1] = results[j]
  #for j in range(LIM):
  if N1>LIM*NPP1:
    NodeScal[range(LIM*NPP1,N1)] = TetListScalarsSec(np.array(range(LIM*NPP1,N1)),NodeTetList,ScalarVal,AoM)
  return NodeScal
  
def TetListScalarsSec(NodeSEC,NodeTetList,ScalarVal,AoM):
  NN1 = NodeSEC.shape[0]
  ScalarSec = np.zeros((NN1,))
  for i in range(NN1):
    j = NodeSEC[i]
    if AoM==0:
      ScalarSec[i] = np.average(ScalarVal[NodeTetList[j]])
    else:
      ScalarSec[i] = np.min(ScalarVal[NodeTetList[j]])
  return ScalarSec
    
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  