# -*- coding: utf-8 -*-
#

import numpy as np                                                                                   
from enthought.mayavi.mlab import *
from fileread import *
from pmeshdeform import *
import ptetmeshdeform as ptet
from scipy.io import *
from scipy.interpolate import Rbf
from scipy.spatial import KDTree
from scipy.stats.mstats import find_repeats
from ptrimeshalt import *
from ptrifeatures import *
from plinedeform import *
import pshapecontext as sh
import pmeshquality as qu
import time


### FEMUR MODELS:
### ICP
NC1 = np.ma.load('femur1NC')
NC2 = np.ma.load('femur2NC')
Tri1 = np.ma.load('femur1tri')
Tri2 = np.ma.load('femur2tri')


NC1 = np.ma.load('femur1RS1T')
NCrst,Conv = lineRST(NC2,[],[],NC1,[],[],0,1,0)	#(NCB,RlinesB,VlinesB,NCT,RlinesT,VlinesT,UseFeat,UseScale)
np.ma.dump(NCrst,'femur1RST')
#np.ma.dump(Conv,'femur1RT_Conv')

### Full registration:
gamm=2
sigm0=10
f=1.0715
#gammL = [1.2,1.5,1.8,2,2.2,2.5,2.8] #7
#sigm0L = [0.5,2,5,10,20,30,50] #7
#fL = [1.,1.1,1.3,1.5,1.7,1.9,2.] # 7
#NCT = np.ma.load('femur1RS1T')
#NCelas,Conv = elasticsurf(NC2,Tri2,np.array([]),np.array([]),np.array(range(Tri2.shape[0])),NCT,Tri1,np.array(range(Tri1.shape[0])),np.array(range(NC2.shape[0])),np.array(range(NCT.shape[0])),75,np.array([]),gamm,50,f)
#np.ma.dump(NCelas,'femurReg_gamm7_i75')
#NCelas,Conv = elasticsurf(NC2,Tri2,np.array([]),np.array([]),np.array(range(Tri2.shape[0])),NCT,Tri1,np.array(range(Tri1.shape[0])),np.array(range(NC2.shape[0])),np.array(range(NCT.shape[0])),40,np.array([]),gamm,2,1.1)
#np.ma.dump(NCelas,'femurReg_f1_i40')
#NCelas,Conv = elasticsurf(NC2,Tri2,np.array([]),np.array([]),np.array(range(Tri2.shape[0])),NCT,Tri1,np.array(range(Tri1.shape[0])),np.array(range(NC2.shape[0])),np.array(range(NCT.shape[0])),20,np.array([]),gamm,2,1.3)
#np.ma.dump(NCelas,'femurReg_f2_i20')
#np.ma.dump(Conv,'femurConv_f'+str(i))
#for i in range(7):
  #NCelas,Conv = elasticsurf(NC2,Tri2,np.array([]),np.array([]),np.array(range(Tri2.shape[0])),NCT,Tri1,np.array(range(Tri1.shape[0])),np.array(range(NC2.shape[0])),np.array(range(NCT.shape[0])),100,np.array([]),gammL[i],sigm0,f)
  #np.ma.dump(NCelas,'femurReg_gamm'+str(i))
  #np.ma.dump(Conv,'femurConv_gamm'+str(i))
#for i in range(7):
  #NCelas,Conv = elasticsurf(NC2,Tri2,np.array([]),np.array([]),np.array(range(Tri2.shape[0])),NCT,Tri1,np.array(range(Tri1.shape[0])),np.array(range(NC2.shape[0])),np.array(range(NCT.shape[0])),100,np.array([]),gamm,sigm0L[i],f)
  #np.ma.dump(NCelas,'femurReg_sig'+str(i))
  #np.ma.dump(Conv,'femurConv_sig'+str(i))
#for i in range(7):
  #NCelas,Conv = elasticsurf(NC2,Tri2,np.array([]),np.array([]),np.array(range(Tri2.shape[0])),NCT,Tri1,np.array(range(Tri1.shape[0])),np.array(range(NC2.shape[0])),np.array(range(NCT.shape[0])),100,np.array([]),gamm,sigm0,fL[i])
  #np.ma.dump(NCelas,'femurReg_f'+str(i))
  #np.ma.dump(Conv,'femurConv_f'+str(i))


#NCR = np.ma.load('femur1RS1T')
#figure(fgcolor=(0,0,0),bgcolor=(1,1,1))
##triangular_mesh(NC1[:,0],NC1[:,2],NC1[:,1],Tri1,color=(0.8,0.6,0.6),opacity=1)
#triangular_mesh(NC2[:,0],NC2[:,2],NC2[:,1],Tri2,color=(0,0,0),opacity=0.4,representation='wireframe',line_width=1)
#triangular_mesh(NCR[:,0],NCR[:,2],NCR[:,1],Tri1,color=(0.8,0.6,0.6),opacity=1)
#for i in range(47):
  #NC1 = np.ma.load('Femur1NC_'+str(i+1))
  #AA = triangular_mesh(NC1[:,0],NC1[:,2],NC1[:,1],Tri1,color=(0.8,0.6,0.6),opacity=1)
  #AA.visible= 1==2
#for i in range(99):
  #NC2 = np.ma.load('Femur2NC_'+str(i+2))
  #AA = triangular_mesh(NC2[:,0],NC2[:,2],NC2[:,1],Tri2,color=(0,0,0),opacity=0.4,representation='wireframe',line_width=1)
  #AA.visible= 1==2