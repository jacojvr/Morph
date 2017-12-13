# -*- coding: utf-8 -*-
AA = ['000','030','035','040','045','050','055','060','065','070','100']

import numpy as np
import pyvtk as pv
#strSplit = 'Data = element stress xy, yz and xz\n'
#for inc in range(11):
  #print inc
  #fid = open('SkullPO_'+AA[inc]+'stress.log')
  #SS = fid.read()
  #SS = SS.split(strSplit)
  #SS1 = SS[0]
  #SS2 = SS[1]
  #SS1 = SS1.split('\n')
  #SS2 = SS2.split('\n')
  #SS1 = SS1[0:1687791]
  #SS2 = SS2[0:1687791]
  #for i in range(1687791):                                                                                                                                                        
    #SS1[i] = SS1[i].split(' ')
    #SS2[i] = SS2[i].split(' ')
  #ST1 = np.array(SS1,float)
  #ST2 = np.array(SS2,float)
  #Stress = np.c_[ST1[:,1:],ST2[:,1:]]
  #vonMisses = np.sqrt(Stress[:,0]*Stress[:,0]+Stress[:,1]*Stress[:,1]+Stress[:,2]*Stress[:,2]-Stress[:,1]*Stress[:,2]-Stress[:,0]*Stress[:,2]-Stress[:,0]*Stress[:,1]+3*(Stress[:,3]*Stress[:,3]+Stress[:,4]*Stress[:,4]+Stress[:,5]*Stress[:,5]))
  #Stress = np.c_[Stress,vonMisses]
  #np.ma.dump(Stress,'Skull'+AA[inc]+'Stress')

Tet = np.ma.load('SkullPO_TetUse')
NC = np.ma.load('SkullPO_NC050')
for inc in range(11):
  print 'Skull '+AA[inc]+' on Average Shape'
  Stress = np.ma.load('Skull'+AA[inc]+'Stress')
  print Stress[:,6]
  skvtk = pv.VtkData(pv.UnstructuredGrid(points = NC,tetra=Tet),'skull '+AA[inc]+' vM stress',pv.CellData(pv.Scalars(Stress[:,6],'von Misses')))
  skvtk.tofile('SkullAver_'+AA[inc]+'vM.vtk')