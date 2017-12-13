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
#import pmeshquality as qu
import time

NCB,TriB = readmsh('BSkull2.msh')
NCB = np.ma.load('BSKull_Symm')
NCT,TriT = readmsh('Skull4F.msh')

keepFA = np.ma.load('BSkullRemoveInternalNC2')

#ridgeN,valleyN,Kmax,Kmin,eMax,eMin,PdMax,PdMin = crestlines(np.array(range(NCT.shape[0])),NCT,TriT)	#(nodes,NC,Tri,NN,ls,neighSpat,eOrk)
###Skull4:
#np.ma.dump(ridgeN,'Skull4RNodes100')
#np.ma.dump(valleyN,'Skull4VNodes100')
#np.ma.dump(Kmax,'Skull4Kmax100')
#np.ma.dump(Kmin,'Skull4Kmin100')
#np.ma.dump(eMax,'Skull4Emax100')
#np.ma.dump(eMin,'Skull4Emin100')
#np.ma.dump(PdMax,'Skull4PDmax100')
#np.ma.dump(PdMin,'Skull4PDmin100')

Kmax = np.ma.load('Skull4Kmax100')
Kmin = np.ma.load('Skull4Kmin100')
eMax = np.ma.load('Skull4Emax100')
eMin = np.ma.load('Skull4Emin100')
PdMax = np.ma.load('Skull4PDmax100')
PdMin = np.ma.load('Skull4PDmin100')
ridgeN = np.ma.load('Skull4RNodes100')
valleyN = np.ma.load('Skull4VNodes100')



valleyN = valleyN[np.where(Kmin[valleyN]<1.8)[0]]
ridgeN = ridgeN[np.where(Kmax[ridgeN]>1.8)[0]]
#print '############ Valley Lines ##################'
#valleyLN = LineDraw(valleyN,NCT,TriT,PdMax,PdMin,-Kmin)
#np.ma.dump(valleyLN,'Skull4VLines100K')
#print '############ Ridge Lines ##################'
#ridgeLN = LineDraw(ridgeN,NCT,TriT,PdMin,PdMax,Kmax)
#np.ma.dump(ridgeLN,'Skull4RLines100K')

#ridgeLN = np.ma.load('Skull4RLines100K')
#valleyLN = np.ma.load('Skull4VLines100K')

#print '############# Threshold ##############'
#valleyLNT = ThreshLines(valleyLN,10,NCT,-Kmin)
#np.ma.dump(valleyLNT,'Skull4VLines100T1')
#ridgeLNT = ThreshLines(ridgeLN,1,NCT,Kmax)
#np.ma.dump(ridgeLNT,'Skull4RLines100T1')
Kmax = np.ma.load('BSkullSymmKmax3NN')
Kmin = np.ma.load('BSkullSymmKmin3NN')
#eMax = np.ma.load('BSkullSymmEmax3NN')
#eMin = np.ma.load('BSkullSymmEmin3NN')
#PdMax = np.ma.load('BSkullSymmPDmax3NN')
#PdMin = np.ma.load('BSkullSymmPDmin3NN')
#ridgeN = np.ma.load('BSkullSymmRNodes200')
#valleyN = np.ma.load('BSkullSymmVNodes200')
KmaxT = np.ma.load('Skull4Kmax100')
KminT = np.ma.load('Skull4Kmin100')

#keep = np.where((np.abs(Kmin)>0.1)|(Kmax>0.1))[0]

ridgeLN1,valleyLN1 = np.ma.load('BSkullRLineskeepFA2'),np.ma.load('BSkullVLineskeepFA2')
#ridgeLN = np.ma.load('Skull4RLines100K')
#valleyLN = np.ma.load('Skull4VLines100K')
##print ridgeLN[0],valleyLN[0]

ridgeLN2 = np.ma.load('Skull4RLines100T1')
valleyLN2 = np.ma.load('Skull4VLines100T1')
ridgeLNT,valleyLNT = ridgeLN2,valleyLN2
##print ridgeLNT[0],valleyLNT[0]
##np.ma.dump(ridgeLNT,'BSkullRLines100KT50')
#ridgeLN1 = np.ma.load('BSkullSymmRLines200T1')
#valleyLN1 = np.ma.load('BSkullSymmVLines200T1')
#NCT = NCT*10
#transl = np.sum(NCT,0)/NCT.shape[0]
#NCT = NCT-transl
#NCT = np.ma.load('Skull4Translate')
##Rotate scale and translate
#NCrst = lineRST(NCB,ridgeLN1,valleyLN1,NCT,ridgeLN2,valleyLN2,0,1)	#(NCB,RlinesB,VlinesB,NCT,RlinesT,VlinesT,UseFeat,UseScale)
#np.ma.dump(NCrst,'Skull4_NCmoveRST')

NCT = np.ma.load('Skull4_NCmoveRST')
NCS = NCT
### Only keep registered nodes in outer area:
##RegB = np.ma.load('BSkullRemoveInternalNC')

#NCS,keepB,keepT = LineReg(NCB,TriB,ridgeLN1,valleyLN1,NCT,TriT,ridgeLN2,valleyLN2,50,5) #(NCB,TriB,RlinesB,VlinesB,NCT,TriT,RlinesT,VlinesT,percentReg,DistMax,RegB)
#np.ma.dump(NCS,'BSkullto4_deformRI2_2')
#np.ma.dump(keepB,'BSkullto4_keepRI2_2')
#np.ma.dump(keepT,'Skull4fromB_keepRI2_2')
NCS = np.ma.load('BSkullto4_deformRI2_2')
keepT = np.ma.load('Skull4fromB_keepRI2_2')
keepB = np.ma.load('BSkullto4_keepRI2_2')
##keepB = np.array([])
##for i in keepBI:
  ##if np.where(keepFA==i)[0].size>0:
    ##keepB=np.r_[keepB,i]
##keepB = np.array(keepB,int)
##print keepBI.shape, find_repeats(keepBI)
##print keepB.shape, find_repeats(keepB)
FeatNodesB = keepFA[np.where((Kmax[keepFA]>0.08)|(np.abs(Kmin[keepFA])>0.08))[0]]
FeatLessNodesB = np.where((Kmax<0.08)&(np.abs(Kmin)<0.08))[0]
FeatNodesT = np.where((KmaxT>2)|(np.abs(KminT)>2))[0]
FeatLessNodesT = np.where((KmaxT<2)&(np.abs(KminT)<2))[0]
##print 'Base nodes and triangles allowed'
#UseN_B,RemoveB,AllowableB = AllowableNaT(NCB,TriB,FeatNodesB,FeatLessNodesB,ridgeLN1,valleyLN1,keepB)
#np.ma.dump(UseN_B,'BSkullto4_useNRI2_2')
#np.ma.dump(AllowableB,'BSkullto4_AlloweTriRI2_2')
##UseN_B = np.array(range(NCB.shape[0]))
UseN_B = np.ma.load('BSkullto4_useNRI2_2')
##AllowableB = np.ma.load('BSkullto5_AlloweTriRI')
AllowableB = np.ma.load('BSkull_AlloweTriRI')
#FeatNodesB = UseN_B[np.where((Kmax[UseN_B]>0.1)|(np.abs(Kmin[UseN_B])>0.1))[0]]
##print 'Target nodes and triangles allowed'
#UseN_T,RemoveT,AllowableT = AllowableNaT(NCT,TriT,FeatNodesT,FeatLessNodesT,ridgeLN2,valleyLN2,keepT)
#np.ma.dump(UseN_T,'Skull4fromB_useNRI2_2')
#np.ma.dump(RemoveT,'Skull4fromB_removeNRI2_2')
#np.ma.dump(AllowableT,'Skull4fromB_AlloweTriRI2_2')
UseN_T = np.ma.load('Skull4fromB_useNRI2_2')
##RemoveT = np.ma.load('Skull5fromB_removeNRI2')
AllowableT = np.ma.load('Skull4fromB_AlloweTriRI2_2')
##FeatNodesT5 = UseN_T[np.where((KmaxT[UseN_T]>0.18)|(np.abs(KminT[UseN_T])>0.18))[0]]
##FeatNodesT =np.where((KmaxT>0.18)|(np.abs(KminT)>0.18))[0]

####### Load tetrahedral mesh information
NCTnc = np.ma.load('BSkullFTnc')
TetT = np.ma.load('BSkullFTtet')
TriTet = np.ma.load('BSkullFTtri')
outer = np.ma.load('BSkullFTn_outer')
inner = np.ma.load('BSkullFTn_inner')
neighbList = np.ma.load('BSkullFTneighbTri')
neighbListTet = np.ma.load('BSkullFTneighbTet')
LandmB,LandmB_NC = keepB,NCS[keepB,]
##USENORMALS = np.array([1,2,3,11,12,21,22])	# Iterations where normal information is used
#NCdef,Conv = ptet.elasticsurf(NCTnc,TetT,NCB,TriB,LandmB,LandmB_NC,AllowableB,NCT,TriT,AllowableT,UseN_B,UseN_T)
#np.ma.dump(NCdef,'BSkullto4_FulldeformRI2_TetSmooth')
##NCS1 = RBFmorph(NCB,NCB[LandmB],LandmB_NC-NCB[LandmB])
##np.ma.dump(NCS1,'BSkullto5_TPSm10')
##NCS = np.ma.load('TempElasNodes_Iter0_TimeMonOct4')
#NCS = np.ma.load('BSkullto4_FulldeformRI2_GW')
#NCB,TriB = readmsh('BSkullSymm.msh')
##NCS2 = np.ma.load('Bto5RegTemp40_21Octp20TaubinA')
NCS1 = np.ma.load('TempElasNodes_Iter60_TimeThuNov4')
##NCS2 = np.ma.load('Skull5TetSymm20T')
##NCS1 = np.c_[-NCS2[:,0],NCS2[:,1:]]

##EQ = np.ma.load('ElQual1')
##QMAT = np.ma.load('QMAT1')
##NCsmooth = ptet.LaplacMesh(inner,NCS2,neighbListTet,100,1)
##np.ma.dump(NCsmooth,'Skull5TetSymm100intLapl')
###NCopt = qu.OptEQ(NCS2,TetT,inner)
##NCopt = qu.OptINT(inner,outer,NCsmooth,TetT,neighbListTet,30)
##np.ma.dump(NCopt,'Skull5TetSymm10Topt')
##Def = NCS - NCB
##Def = np.sqrt(np.sum(Def*Def,1))
##SI = 0.5 - np.arctan((Kmax+Kmin)/(Kmax-Kmin))/np.pi
##figure(fgcolor=(0,0,0),bgcolor=(1,1,1))
##triangular_mesh(NCB[:,0],NCB[:,2],NCB[:,1],TriB,scalars = SI[:,0],colormap='RdGy')
##points3d(NCB[FP,0],NCB[FP,2],NCB[FP,1],color=(0,0,1),scale_factor=2)

figure(fgcolor=(0,0,0),bgcolor=(1,1,1))
triangular_mesh(NCT[:,0],NCT[:,2],NCT[:,1],TriT,color=(0.8,0.8,0.6),opacity=1)
triangular_mesh(NCB[:,0],NCB[:,2],NCB[:,1],TriB,color=(0.6,0.6,0.8),opacity=1)
#triangular_mesh(NCS[:,0],NCS[:,2],NCS[:,1],TriB,color=(0.6,0.8,0.8),opacity=1)
##triangular_mesh(NCS[:,0],NCS[:,2],NCS[:,1],TriB,scalars = Def)
triangular_mesh(NCS1[:,0],NCS1[:,2],NCS1[:,1],TriTet,color=(0.6,0.8,0.6),opacity=1)
#triangular_mesh(NCS2[:,0],NCS2[:,2],NCS2[:,1],TriTet,color=(0.6,0.8,0.6),opacity=1)

###figure(fgcolor=(0,0,0),bgcolor=(1,1,1))
####points3d(NCB[keepFA,0],NCB[keepFA,2],NCB[keepFA,1],color=(0,0,0),scale_factor=0.5)
#points3d(NCT[FeatNodesT,0],NCT[FeatNodesT,2],NCT[FeatNodesT,1],color=(0,0,0),scale_factor=0.5)
##points3d(NCT[FeatNodesT5,0],NCT[FeatNodesT5,2],NCT[FeatNodesT5,1],color=(1,0,0),scale_factor=0.6)
###points3d(NCT[T5Fkeep,0],NCT[T5Fkeep,2],NCT[T5Fkeep,1],color=(1,0,0),scale_factor=0.5)
#points3d(NCB[FeatNodesB,0],NCB[FeatNodesB,2],NCB[FeatNodesB,1],color=(0,0,1),scale_factor=0.5)
###points3d(NCB[AllowableB,0],NCB[AllowableB,2],NCB[AllowableB,1],color=(0,0,1),scale_factor=0.8)
###points3d(NCB[BFkeep,0],NCB[BFkeep,2],NCB[BFkeep,1],color=(1,0,0),scale_factor=0.5)

figure(fgcolor=(0,0,0),bgcolor=(1,1,1))
##points3d(NCT[keepT,0],NCT[keepT,2],NCT[keepT,1],color=(1,0,0),scale_factor=1)
##points3d(NCB[keepB,0],NCB[keepB,2],NCB[keepB,1],color=(0,0,1),scale_factor=1)
##points3d(NCS[keepB,0],NCS[keepB,2],NCS[keepB,1],color=(0,0,1),scale_factor=1)
###NCT=NCB
###figure(fgcolor=(0,0,0),bgcolor=(1,1,1))
#points3d(NCT[ridgeN,0],NCT[ridgeN,2],NCT[ridgeN,1],color=(1,0,0),scale_factor=1)
#points3d(NCT[valleyN,0],NCT[valleyN,2],NCT[valleyN,1],color=(0,0,1),scale_factor=1)
###points3d(NCS[rn,0],NCS[rn,2],NCS[rn,1],color=(1,0,0),scale_factor=0.1,opacity=0.2)
##points3d(NCS[vn,0],NCS[vn,2],NCS[vn,1],color=(0,0,1),scale_factor=0.1,opacity=0.2)
##points3d(NCS[r2,0],NCS[r2,2],NCS[r2,1],color=(1,0,0),scale_factor=0.4)#,opacity=0.2)
##points3d(NCS[v2,0],NCS[v2,2],NCS[v2,1],color=(0,0,1),scale_factor=0.4)#,opacity=0.2)
##points3d(NCS[keep1,0],NCS[keep1,2],NCS[keep1,1],color=(0,0,0),scale_factor=0.1,opacity=0.2)
##points3d(NCS[keep2,0],NCS[keep2,2],NCS[keep2,1],color=(0,0,0),scale_factor=0.1,opacity=0.2)
##points3d(NCS[keep3,0],NCS[keep3,2],NCS[keep3,1],color=(0,0,0),scale_factor=0.1,opacity=0.2)
##points3d(NCS[keep4,0],NCS[keep4,2],NCS[keep4,1],color=(0,0,0),scale_factor=0.1,opacity=0.2)
##points3d(NCS[keep5,0],NCS[keep5,2],NCS[keep5,1],color=(0,0,0),scale_factor=0.1,opacity=0.2)
##points3d(NCS[keep6,0],NCS[keep6,2],NCS[keep6,1],color=(0,0,0),scale_factor=0.1,opacity=0.2)
##KA = 1000*np.average(Kmax[ridgeN])
##KB = 10000*np.average(-Kmin[valleyN])
##count1 = 0
###count2 = 0
#for i in range(1,ridgeLN[0]+1):
  #####Thick = np.max(Kmax[ridgeLN[i][0]])
  ##Th = lineThresh(ridgeLN[i][0],NCT,Kmax)
  ##print Th
  #######Thick = Th/ridgeLN[i][0].size
  #######if Th>10:
    #######count1=count1+1
  ##if (ridgeLN[i][0].size>3)&(Th>30)&(Th<100):
    ##plot3d(NCB[ridgeLN[i][0],0],NCB[ridgeLN[i][0],2],NCB[ridgeLN[i][0],1],color=(1,0,0),representation='wireframe',line_width=1)
  ##if (ridgeLN[i][0].size>3)&(Th>100)&(Th<200):
    ##plot3d(NCB[ridgeLN[i][0],0],NCB[ridgeLN[i][0],2],NCB[ridgeLN[i][0],1],color=(1,0,0),representation='wireframe',line_width=2)#Thick/KA)
  ##if (ridgeLN[i][0].size>3)&(Th>200):
    ##plot3d(NCB[ridgeLN[i][0],0],NCB[ridgeLN[i][0],2],NCB[ridgeLN[i][0],1],color=(1,0,0),representation='wireframe',line_width=4)#Thick/KA)

##for i in range(1,valleyLN[0]+1):
  #####Thick = np.max(-Kmin[valleyLN[i][0]])
  ##Th = lineThresh(valleyLN[i][0],NCT,-Kmin)*5
  ###print Th
  #######if Th>10:
    #########count2=count2+1
  #########Thick = Th/valleyLN[i][0].size
  ##if (valleyLN[i][0].size>3):
    ##plot3d(NCB[valleyLN[i],0],NCB[valleyLN[i],2],NCB[valleyLN[i],1],color=(0,0,1),representation='wireframe',line_width=2)
  ##if (valleyLN[i][0].size>3)&(Th>30)&(Th<100):
    ##plot3d(NCB[valleyLN[i][0],0],NCB[valleyLN[i][0],2],NCB[valleyLN[i][0],1],color=(0,0,1),representation='wireframe',line_width=1)
  ##if (valleyLN[i][0].size>3)&(Th>100)&(Th<200):
    ##plot3d(NCB[valleyLN[i][0],0],NCB[valleyLN[i][0],2],NCB[valleyLN[i][0],1],color=(0,0,1),representation='wireframe',line_width=2)
  ##if (valleyLN[i][0].size>3)&(Th>200):
    ##plot3d(NCB[valleyLN[i][0],0],NCB[valleyLN[i][0],2],NCB[valleyLN[i][0],1],color=(0,0,1),representation='wireframe',line_width=4)#Thick/KB)
####print count1,count2
## Thresholded lines
##NCS=NCB
#for i in range(1,ridgeLNT[0]+1):
  ##Th = lineThresh(ridgeLNT[i],NCS,Kmax)
  #plot3d(NCS[ridgeLNT[i],0],NCS[ridgeLNT[i],2],NCS[ridgeLNT[i],1],color=(1,0,0),representation='wireframe',line_width=2)
#for i in range(1,valleyLNT[0]+1):
  ##Th = lineThresh(valleyLNT[i],NCS,-Kmin)
  #plot3d(NCS[valleyLNT[i],0],NCS[valleyLNT[i],2],NCS[valleyLNT[i],1],color=(0,0,1),representation='wireframe',line_width=2)


#### To call FEBio:
#### import os
####os.system('~/Software/MRL/FEBio/febio.lnx -i ~/Documents/MASTERS/Skulls/Prognathic/PrognI.feb') for exemple

#NC0 = NCS
NC1 = NCS1
#NC2 = NCS2
##NC0 = np.ma.load('BSkullto5_Fulldeform')
##NC1 = np.ma.load('BSkullto5_FulldeformRI2')
##NC0 = np.ma.load('BSkull_RBFdeform')
##NC1 = np.ma.load('TempBto5_Iter10')*10
##NC2 = np.ma.load('TempBto5_Iter20')*10
##NC3 = np.ma.load('TempBto5_Iter30')*10
##NC4 = np.ma.load('TempBto5_Iter40')*10
##NC5 = np.ma.load('TempBto5_Iter50')*10
##NC6 = np.ma.load('TempBto5_Iter60')*10
##NC7 = np.ma.load('TempBto5_Iter70')*10
##NC8 = np.ma.load('TempBto5_Iter80')*10
##NC9 = np.ma.load('TempBto5_Iter90')*10
##NC10 = np.ma.load('BSkullto5_Fulldeform')
##NC11 = np.ma.load('TempElasNodes_Iter11_TimeWedSep15')*10
##NC12 = np.ma.load('TempElasNodes_Iter12_TimeWedSep15')*10
##NC13 = np.ma.load('TempElasNodes_Iter13_TimeWedSep15')*10
targeT = triangular_mesh(NCT[:,0],NCT[:,2],NCT[:,1],TriT,color=(0.9,0.5,0.5),opacity=0)
basE = triangular_mesh(NCB[:,0],NCB[:,2],NCB[:,1],TriB,color=(0.5,0.9,0.5),opacity=0)
#plt0 = triangular_mesh(NC0[:,0],NC0[:,2],NC0[:,1],TriB,color=(0.5,0.9,0.5),opacity=0)
plt1 = triangular_mesh(NC1[:,0],NC1[:,2],NC1[:,1],TriTet,color=(0.5,0.9,0.5),opacity=0)
#plt2 = triangular_mesh(NC2[:,0],NC2[:,2],NC2[:,1],TriTet,color=(0.5,0.9,0.5),opacity=0)
##plt3 = triangular_mesh(NC3[:,0],NC3[:,2],NC3[:,1],TriB,color=(0.5,0.9,0.5),opacity=0)
##plt4 = triangular_mesh(NC4[:,0],NC4[:,2],NC4[:,1],TriB,color=(0.5,0.9,0.5),opacity=0)
##plt5 = triangular_mesh(NC5[:,0],NC5[:,2],NC5[:,1],TriB,color=(0.5,0.9,0.5),opacity=0)
##plt6 = triangular_mesh(NC6[:,0],NC6[:,2],NC6[:,1],TriB,color=(0.5,0.9,0.5),opacity=0)
##plt7 = triangular_mesh(NC7[:,0],NC7[:,2],NC7[:,1],TriB,color=(0.5,0.9,0.5),opacity=0)
##plt8 = triangular_mesh(NC8[:,0],NC8[:,2],NC8[:,1],TriB,color=(0.5,0.9,0.5),opacity=0)
##plt9 = triangular_mesh(NC9[:,0],NC9[:,2],NC9[:,1],TriB,color=(0.5,0.9,0.5),opacity=0)
##plt10 = triangular_mesh(NC10[:,0],NC10[:,2],NC10[:,1],TriB,color=(0.5,0.9,0.5),opacity=0)
###plt11 = triangular_mesh(NC11[:,0],NC11[:,2],NC11[:,1],TriB,color=(0.5,0.9,0.5),opacity=0)
###plt12 = triangular_mesh(NC12[:,0],NC12[:,2],NC12[:,1],TriB,color=(0.5,0.9,0.5),opacity=0)
###plt13 = triangular_mesh(NC13[:,0],NC13[:,2],NC13[:,1],TriB,color=(0.5,0.9,0.5),opacity=0)

cutpT = pipeline.scalar_cut_plane(targeT,color=(1,0,0),line_width=5,plane_orientation='z_axes',name='Target')
cutpT.implicit_plane.origin = (0,0,-40)
cutpT.implicit_plane.widget.enabled = False
cutpB = pipeline.scalar_cut_plane(basE,color=(0,0,0),line_width=3,plane_orientation='z_axes',name='Base')
cutpB.implicit_plane.origin = (0,0,-40)
cutpB.implicit_plane.widget.enabled = False
#cutp0 = pipeline.scalar_cut_plane(plt0,color=(0,0,1),line_width=3,plane_orientation='z_axes',name='iter_00')
#cutp0.implicit_plane.origin = (0,0,-40)#(20,0,0)
#cutp0.implicit_plane.widget.enabled = False
cutp1 = pipeline.scalar_cut_plane(plt1,color=(0,1,1),line_width=3,plane_orientation='z_axes',name='iter_01')
cutp1.implicit_plane.origin = (0,0,-40)#(20,0,0)
cutp1.implicit_plane.widget.enabled = False
#cutp2 = pipeline.scalar_cut_plane(plt2,color=(0,1,0),line_width=3,plane_orientation='z_axes',name='iter_01')
#cutp2.implicit_plane.origin = (0,0,-40)#(20,0,0)
#cutp2.implicit_plane.widget.enabled = False
##cutp3 = pipeline.scalar_cut_plane(plt3,color=(0,0,0),line_width=3,plane_orientation='z_axes',name='iter_01')
##cutp3.implicit_plane.origin = (0,0,-40)#(20,0,0)
##cutp3.implicit_plane.widget.enabled = False
##cutp4 = pipeline.scalar_cut_plane(plt4,color=(0,0,0),line_width=3,plane_orientation='z_axes',name='iter_01')
##cutp4.implicit_plane.origin = (0,0,-40)#(20,0,0)
##cutp4.implicit_plane.widget.enabled = False
##cutp5 = pipeline.scalar_cut_plane(plt5,color=(0,0,0),line_width=3,plane_orientation='z_axes',name='iter_01')
##cutp5.implicit_plane.origin = (0,0,-40)#(20,0,0)
##cutp5.implicit_plane.widget.enabled = False
##cutp6 = pipeline.scalar_cut_plane(plt6,color=(0,0,0),line_width=3,plane_orientation='z_axes',name='iter_01')
##cutp6.implicit_plane.origin = (0,0,-40)#(20,0,0)
##cutp6.implicit_plane.widget.enabled = False
##cutp7 = pipeline.scalar_cut_plane(plt7,color=(0,0,0),line_width=3,plane_orientation='z_axes',name='iter_01')
##cutp7.implicit_plane.origin = (0,0,-40)#(20,0,0)
##cutp7.implicit_plane.widget.enabled = False
##cutp8 = pipeline.scalar_cut_plane(plt8,color=(0,0,0),line_width=3,plane_orientation='z_axes',name='iter_01')
##cutp8.implicit_plane.origin = (0,0,-40)#(20,0,0)
##cutp8.implicit_plane.widget.enabled = False
##cutp9 = pipeline.scalar_cut_plane(plt9,color=(0,0,0),line_width=3,plane_orientation='z_axes',name='iter_01')
##cutp9.implicit_plane.origin = (0,0,-40)#(20,0,0)
##cutp9.implicit_plane.widget.enabled = False
##cutp10 = pipeline.scalar_cut_plane(plt10,color=(0,0,0),line_width=3,plane_orientation='z_axes',name='iter_01')
##cutp10.implicit_plane.origin = (0,0,-40)#(20,0,0)
##cutp10.implicit_plane.widget.enabled = False
##cutp11 = pipeline.scalar_cut_plane(plt11,color=(0,0,0),line_width=3)#,plane_orientation='z_axes',name='iter_01')
##cutp11.implicit_plane.origin = (20,0,0)
##cutp11.implicit_plane.widget.enabled = False
##cutp12 = pipeline.scalar_cut_plane(plt12,color=(0,0,0),line_width=3)#,plane_orientation='z_axes',name='iter_01')
##cutp12.implicit_plane.origin = (20,0,0)
##cutp12.implicit_plane.widget.enabled = False
##cutp13 = pipeline.scalar_cut_plane(plt13,color=(0,0,0),line_width=3)#,plane_orientation='z_axes',name='iter_01')
##cutp13.implicit_plane.origin = (20,0,0)
###cutp13.implicit_plane.widget.enabled = False