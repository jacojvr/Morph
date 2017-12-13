# -*- coding: utf-8 -*-
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
import writeFEBio as wf
import time


## DONE_TO_TEST_UNIQUENESS

#Tri = 
NCTnc = np.ma.load('SkTetAverNC')
TetT = np.ma.load('SkTetAverTET')
TriTet = np.ma.load('BSkullFTtri')
outer = np.ma.load('BSkullFTn_outer')
inner = np.ma.load('BSkullFTn_inner')
neighbList = np.ma.load('BSkullFTneighbTri')
neighbListTet = np.ma.load('BSkullFTneighbTet')

#NCTt = np.ma.load('Skull5TetSymm20T')
NCTt = np.ma.load('Skull4Symm20T')
NCT = NCTt[outer,]
NCB = NCTnc[outer,]
TriB = readmsh('BSkull2.msh')[1]

LandmB = np.ma.load('BSkullto4_keepRI2_2')
LandmB_NC = NCT[LandmB,]
usen = np.array(range(outer.size))
useTri = np.array(range(TriB.shape[0]))
USENORMALS = np.array([1,2,3,11,12,21,22])	# Iterations where normal information is used
gamm = [2, 2]
sigm0 = [10, 20]
i=0
#NCdef,Conv = elasticsurf(NCB,TriB,LandmB,LandmB_NC,useTri,NCT,TriB,useTri,usen,usen,30*(i+2),USENORMALS,gamm[i],sigm0[i],f=1.0715)
for i in [0,1]:
  #NCB = np.ma.load('TempElasNodes_Iter0_TimeWedMar09_2011_20')
  NCB = np.ma.load('TempElasNodes_Iter0_TimeWedMar14_2011_20')
  NCdef,Conv = elasticsurf(NCB,TriB,LandmB,LandmB_NC,useTri,NCT,TriB,useTri,usen,usen,30*(i+2),USENORMALS,gamm[i],sigm0[i],f=1.0715)
  np.ma.dump(Conv,'SkullUnique2_Conv_'+str(i))
  Disp = NCdef - NCTnc[outer]
  Disp = Disp/10
  wf.NodeDisp(NCTnc,TetT,np.array([]),outer,Disp,'SkullOS2_sig'+str(sigm0[i])+'_iter'+str(30*(i+2))+'_')
 
import os  
os.system('~/Software/MRL/FEBio/febio.lnx -i ./SkullOS2_sig20_iter90_.feb')
os.system('~/Software/MRL/FEBio/febio.lnx -i ./SkullOS2_sig10_iter60_.feb')