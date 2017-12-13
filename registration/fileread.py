# -*- coding: utf-8 -*-

import numpy as np
def readmsh(filename):
  fil = open(filename)
  
  # Nodal coordinate array:
  temp = fil.readline()	
  temp = fil.readline()	
  temp = fil.readline()	
  temp = fil.readline()	
  temp = fil.readline()	
  temp = int(temp.rstrip('\r\n'))
  nodalC = np.zeros((temp,3))
  for i in range(0,temp):
    data = fil.readline()
    data = data.rstrip('\r\n')
    data = data.split(' ')
    nodalC[i,:] = np.array(data[1:4],float)
  
  temp = fil.readline()#Get triangle connectivity
  temp = fil.readline()	
  temp = fil.readline()
  temp = int(temp.rstrip('\r\n'))
  TriConnect = np.zeros((temp,3))
  for i in range(0,temp):
    data = fil.readline()
    data = data.rstrip('\r\n')
    data = data.split(' ')
    TriConnect[i,:] = np.array(data[6:9],int)
  TriConnect=np.array(TriConnect-1,int)

  #temp = fil.readline()#Get landmarks
  #temp = fil.readline()	
  #temp = fil.readline()
  #temp = int(temp.rstrip('\r\n'))
  #Land = np.zeros((temp,1))
  #for i in range(0,temp):
    #data = fil.readline()
    #data = data.rstrip('\r\n')
    #Land[i,:] = np.array(data,int)

  #Land = np.array(Land-1,int)
  fil.close
  return nodalC,TriConnect

  
def readEleMsh(filename):
  fil = open(filename)
  NCs = []
  #Nodal coordinate array:
  temp = fil.readline()	
  temp = fil.readline()
  data = fil.readline()
  while data[0:3]<>'end':
    data = data.rstrip('\n')
    data = data.split(' ')
    NCs = NCs + [data[-3:]]
    data = fil.readline()
  data = fil.readline()
  data = fil.readline()
  data = data.rstrip('\n')
  data = data.split(' ')
  data = np.array(data)
  Tets = np.c_[np.array(data[data<>''][-4:],int)]
  data = fil.readline()
  while data[0:3]<>'end':
    data = data.rstrip('\n')
    data = data.split(' ')
    data = np.array(data)
    Tets = np.c_[Tets,np.array(data[data<>''][-4:],int)]
    data = fil.readline()
  fil.close
  return np.array(NCs,float),np.array(Tets.T,int)-1