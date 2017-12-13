# -*- coding: utf-8 -*-
import numpy as np
import scipy as sp
from scipy.spatial import KDTree
from scipy.interpolate import Rbf
def RBFmorph(NCB,knownC,knownD,typeRBF):
  #Determine displacement of Base nodes NCB by displacing nodes at coordinates 
  #knownC by known landmark displacement
  rbfx = Rbf(knownC[:,0],knownC[:,1],knownC[:,2],knownD[:,0],function=typeRBF)
  rbfy = Rbf(knownC[:,0],knownC[:,1],knownC[:,2],knownD[:,1],function=typeRBF)
  rbfz = Rbf(knownC[:,0],knownC[:,1],knownC[:,2],knownD[:,2],function=typeRBF)
  DSx = rbfx(NCB[:,0],NCB[:,1],NCB[:,2])
  DSy = rbfy(NCB[:,0],NCB[:,1],NCB[:,2])
  DSz = rbfz(NCB[:,0],NCB[:,1],NCB[:,2])
  ##print np.c_[DSx,DSy,DSz]
  ##print NCB
  NCI = np.c_[DSx,DSy,DSz]+NCB
  return NCI