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
import pprocess
import time
import pyvtk as pv
import os
LIM=8
Steps = 10
Deg = 0.15

fname = 'SkullPrognSymmSurf'
NCTnc = np.ma.load('SkTetAverNC')
TetT = np.ma.load('SkTetAverTET')
inner = np.ma.load('SkAinnerNdlist')
outer = np.ma.load('BSkullFTn_outer')
#neighbListTet = np.ma.load('SkAverNeighbListTet')
N1 = NCTnc.shape[0]
NPP1 = N1/LIM
#neighbListTet = [[0]]*N1
#results = pprocess.Map(limit=LIM)
#calc = results.manage(pprocess.MakeParallel(ptet.Get1neigh))
#for j in range(0,LIM):
    #calc(np.array(range(0,NPP1))+j*NPP1,NCTnc,TetT)
#for j in range(0,LIM):
    #neighbListTet[j*NPP1:(1+j)*NPP1] = results[j]
#neighbListTet[LIM*NPP1:N1]=ptet.Get1neigh(np.array(range(LIM*NPP1,N1)),NCTnc,TetT)
#np.ma.dump(neighbListTet,'SkAverNeighbListTet')


#NCTnc = np.ma.load('BSkullFTnc')
#TetT = np.ma.load('BSkullFTtet')
#TriTet = np.ma.load('BSkullFTtri')
#outer = np.ma.load('BSkullFTn_outer')
#inner = np.ma.load('BSkullFTn_inner')
#neighbList = np.ma.load('BSkullFTneighbTri')
#neighbListTet = np.ma.load('BSkullFTneighbTet')


NCouter = np.ma.load(fname)
NCS = np.r_[NCTnc]
#NCS[outer,] = NCDef[outer,]
#np.ma.dump(NCS,fname[0:6]+'NC_Ainit')
#NCS = ptet.LaplacMesh(inner,NCS,neighbListTet,100,1)
#np.ma.dump(NCS,fname[0:6]+'NC_AinitLapl')

NCSprev = np.r_[NCS]
DispOuter = (NCouter - NCTnc[outer,])/Steps
for inc in range(Steps):
  #NCS[outer,] = NCS[outer,]+DispOuter	#update boundary displacement
  # Deform mesh using Gaussian smoothing as suggested in paper by R.Bryan et al.
  NCp = NCS[outer,]
  DD = DispOuter
  DS = np.zeros(NCS.shape)
  print 'Do Gaussian Smooth on internal nodes'
  sigma_k2 = np.power(np.power(1.0715,-(inc+1))*10,2)
  results = pprocess.Map(limit=LIM)
  calc = results.manage(pprocess.MakeParallel(ptet.GaussianSmooth))
  for j in range(0,LIM):
    calc(np.array(range(0,NPP1))+j*NPP1,np.r_[NCS],NCp,DD,sigma_k2,2)
  for j in range(0,LIM):
    DS[np.array(range(0,NPP1))+j*NPP1,:] = results[j]
  DS[range(LIM*NPP1,N1),:]=ptet.GaussianSmooth(np.array(range(LIM*NPP1,N1)),np.r_[NCS],NCp,DD,sigma_k2,2)
  NCS[outer,] = NCS[outer,]+DispOuter
  NCS[inner,] = NCS[inner,]+DS[inner,]
  
  
  
  EQ,delt,Sn2,Sig = qu.elemQual_mu(np.array(range(TetT.shape[0])),NCS,TetT)
  print '					Average Element Quality: 	', np.average(EQ)
  print '					Degenerate (q<0.15): 		', np.where(EQ<0.15)[0].size
  print '					Inverted Elements: 		', np.where(Sig<0)[0].size 
  TetDeg = np.where(EQ<0.15)[0]
  DegNd = TetT[TetDeg,]
  DegNd = DegNd.reshape((DegNd.size,))
  DegRep = np.array(find_repeats(DegNd)[0],int)
  for i in DegRep:
    DegNd = DegNd[DegNd<>i]
  DegNd = np.r_[DegNd,DegRep]
  NCS[DegNd,]=NCSprev[DegNd,]+DS[DegNd,]
  #DegNd = np.array(find_repeats(np.r_[DegNd,outer])[0],int)
  DegNd.sort()
  PointConst = np.zeros((NCS.shape[0],))
  PointConst[outer,] = 1
  PointConst[DegNd,] = 0
  PointConst = np.array(PointConst,int)


  print 'Construct VTK object for optimization'
  ##ADD CONTRAINTS AS POINT SCALARS
  skvtk = pv.VtkData(pv.UnstructuredGrid(points = NCS,tetra=TetT),'skull 4 symm',pv.PointData(pv.Scalars(PointConst,'fixed')))
  skvtk.tofile(fname[0:6]+'015constr.vtk')
  os.system('~/Software/Trilinos/Mesquite/msqshape -s ./'+fname[0:6]+'015constr.vtk ./'+fname[0:6]+'015constrOUT.vtk')


  fil = open(fname[0:6]+'015constrOUT.vtk','r+')
  fil.writelines('# vtk DataFile Version 2.0\n')
  fil.close()
  skvtk = pv.VtkData(fname[0:6]+'015constrOUT.vtk')
  NCS = np.array(skvtk.structure.points)
  NCSprev = np.r_[NCS]
  EQ,delt,Sn2,Sig = qu.elemQual_mu(np.array(range(TetT.shape[0])),NCS,TetT)
  print '					Average Element Quality: 	', np.average(EQ)
  print '					Degenerate (q<0.15): 		', np.where(EQ<0.15)[0].size
  print '					Inverted Elements: 		', np.where(Sig<0)[0].size 


#NCTnc = np.ma.load('BSkullFTnc')
#TetT = np.ma.load('BSkullFTtet')
#TriTet = np.ma.load('BSkullFTtri')
#outer = np.ma.load('BSkullFTn_outer')
#inner = np.ma.load('BSkullFTn_inner')
#neighbList = np.ma.load('BSkullFTneighbTri')
#neighbListTet = np.ma.load('BSkullFTneighbTet')
#NCS = np.ma.load(fname)
#skvtk = pv.VtkData(fname[0:6]+'015constrOUT1.vtk')
#NCS = np.array(skvtk.structure.points)

#EQ,delt,Sn2,Sig = qu.elemQual_mu(np.array(range(TetT.shape[0])),NCS,TetT)
#print '					Average Element Quality: 	', np.average(EQ)
#print '					Degenerate (q<0.15): 		', np.where(EQ<0.15)[0].size
#print '					Inverted Elements: 		', np.where(Sig<0)[0].size 
#print 'Remove nodes and triangles on the surface associated with degenerate elements'

#TetDeg = np.where(Sig<0)[0]
##TetDeg = np.where(EQ<0.15)[0]
#DegNd = TetT[TetDeg,]
#DegNd = DegNd.reshape((DegNd.size,))
#DegRep = np.array(find_repeats(DegNd)[0],int)
#for i in DegRep:
  #DegNd = DegNd[DegNd<>i]
#DegNd = np.r_[DegNd,DegRep]
##DegNd = np.array(find_repeats(np.r_[DegNd,outer])[0],int)
#DegNd.sort()

#NNS = np.r_[DegNd]
#for i in DegNd:
  #NNS = np.r_[NNS,neighbListTet[i]]
#DegRep = np.array(find_repeats(NNS)[0],int)
#for i in DegRep:
  #NNS = NNS[NNS<>i]
#NNS = np.r_[NNS,DegRep]


#np.ma.dump(DegNd,'SK4REMNODES')

##NCSsmooth = ptet.MeshSmooth(NNS,np.r_[NCS],neighbListTet,100,1)
##np.ma.dump(NCSsmooth,'SK4NC_SmDeg')
#NCS = np.ma.load('SK4NC_SmDeg')

##for i in DegNd:
  ##if i==DegNd[0]:
    ##TD = np.where(TriTet==i)[0]
  ##else:
    ##TD = np.r_[TD,np.where(TriTet==i)[0]]
    
##TDrep = np.array(find_repeats(TD)[0],int)
##for i in TDrep:
  ##TD = TD[TD<>i]
##TD = np.r_[TD,TDrep]
##keepTri = np.array(range(TriTet.shape[0]),int)
##for i in TD:
  ##keepTri=keepTri[keepTri<>i]
##TriU = TriTet[keepTri,]
##np.ma.dump(TriU,'SK4TriCONSTR')

#keepT = np.array(range(TetT.shape[0]))
#for i in TetDeg:
  #keepT = keepT[keepT<>i]



#PointConst = np.zeros((NCS.shape[0],))
#PointConst[outer,] = 1
#PointConst[NNS,] = 0
#PointConst = np.array(PointConst,int)


#print 'Construct VTK object for optimization'
###ADD CONTRAINTS AS POINT SCALARS
#skvtk = pv.VtkData(pv.UnstructuredGrid(points = NCS,tetra=TetT[keepT]),'skull 4 symm initial',pv.PointData(pv.Scalars(PointConst,'fixed')))
#skvtk.tofile(fname[0:6]+'015constr1.vtk')
#os.system('~/Software/Trilinos/Mesquite/msqshape -s ./'+fname[0:6]+'015constr1.vtk ./'+fname[0:6]+'015constrOUT1.vtk')


#fil = open(fname[0:6]+'015constrOUT1.vtk','r+')
#fil.writelines('# vtk DataFile Version 2.0\n')
#fil.close()
#skvtk = pv.VtkData(fname[0:6]+'015constrOUT1.vtk')
#NCS2 = np.array(skvtk.structure.points)
#EQ2,delt2,Sn22,Sig2 = qu.elemQual_mu(np.array(range(TetT.shape[0])),NCS2,TetT)
#print '					Average Element Quality: 	', np.average(EQ2)
#print '					Degenerate (q<0.15): 		', np.where(EQ2<0.15)[0].size
#print '					Inverted Elements: 		', np.where(Sig2<0)[0].size 
#print
#print 'Construct VTK object for optimization'
###ADD CONTRAINTS AS POINT SCALARS
#skvtk = pv.VtkData(pv.UnstructuredGrid(points = NCS2,tetra=TetT),'skull 4 symm initial',pv.PointData(pv.Scalars(PointConst,'fixed')))
#skvtk.tofile(fname[0:6]+'015constr2.vtk')
#os.system('~/Software/Trilinos/Mesquite/msqshape -s ./'+fname[0:6]+'015constr2.vtk ./'+fname[0:6]+'015constrOUT2.vtk')


#fil = open(fname[0:6]+'015constrOUT2.vtk','r+')
#fil.writelines('# vtk DataFile Version 2.0\n')
#fil.close()
#skvtk = pv.VtkData(fname[0:6]+'015constrOUT2.vtk')
#NCS3 = np.array(skvtk.structure.points)
#EQ3,delt3,Sn23,Sig3 = qu.elemQual_mu(np.array(range(TetT.shape[0])),NCS3,TetT)
#print '					Average Element Quality: 	', np.average(EQ3)
#print '					Degenerate (q<0.15): 		', np.where(EQ3<0.15)[0].size
#print '					Inverted Elements: 		', np.where(Sig3<0)[0].size 
#print
#print 'Construct VTK object for optimization'
###ADD CONTRAINTS AS POINT SCALARS
#skvtk = pv.VtkData(pv.UnstructuredGrid(points = NCS3,tetra=TetT),'skull 4 symm initial',pv.PointData(pv.Scalars(PointConst,'fixed')))
#skvtk.tofile(fname[0:6]+'015constr3.vtk')
#os.system('~/Software/Trilinos/Mesquite/msqshape -s ./'+fname[0:6]+'015constr3.vtk ./'+fname[0:6]+'015constrOUT3.vtk')


#fil = open(fname[0:6]+'015constrOUT3.vtk','r+')
#fil.writelines('# vtk DataFile Version 2.0\n')
#fil.close()
#skvtk = pv.VtkData(fname[0:6]+'015constrOUT3.vtk')
#NCS4 = np.array(skvtk.structure.points)
#EQ4,delt4,Sn24,Sig4 = qu.elemQual_mu(np.array(range(TetT.shape[0])),NCS4,TetT)
#print '					Average Element Quality: 	', np.average(EQ4)
#print '					Degenerate (q<0.15): 		', np.where(EQ4<0.15)[0].size
#print '					Inverted Elements: 		', np.where(Sig4<0)[0].size 