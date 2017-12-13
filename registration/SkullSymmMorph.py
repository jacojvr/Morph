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
import time

NCB,TriB = readmsh('BSkull2.msh')
NCB = np.ma.load('BSKull_Symm')
#NotRegB = readmsh('BSkullNotReg.msh')[0]
#NCB2,TriB2 = readmsh('BSkullNotReg3.msh')[0:2]
#TriC = np.c_[np.sum(np.c_[NCB[TriB[:,0],0],NCB[TriB[:,1],0],NCB[TriB[:,2],0]],1)/3,
    #np.sum(np.c_[NCB[TriB[:,0],1],NCB[TriB[:,1],1],NCB[TriB[:,2],1]],1)/3,np.sum(np.c_[NCB[TriB[:,0],2],NCB[TriB[:,1],2],NCB[TriB[:,2],2]],1)/3]
#TriC2 = np.c_[np.sum(np.c_[NCB2[TriB2[:,0],0],NCB2[TriB2[:,1],0],NCB2[TriB2[:,2],0]],1)/3,
    #np.sum(np.c_[NCB2[TriB2[:,0],1],NCB2[TriB2[:,1],1],NCB2[TriB2[:,2],1]],1)/3,np.sum(np.c_[NCB2[TriB2[:,0],2],NCB2[TriB2[:,1],2],NCB2[TriB2[:,2],2]],1)/3]
#AllowTriB = np.r_[TriB]
#KDTriB = KDTree(TriC,5)
#RegTB = np.array(range(TriB.shape[0]))
#print "Find Allowable Triangles From User specified Removed Surface"
#for i in range(TriC2.shape[0]):
  #neigh =KDTriB.query(TriC2[i,])[1]
  #RegTB = RegTB[np.where(RegTB<>neigh)[0]]
#print "Dump it to file"
#np.ma.dump(RegTB,'BSkull_AlloweTriRI')

#KDTncB = KDTree(NCB,5)
#RegTB = np.array(range(NCB.shape[0]))
#print "Find Allowable Nodes From User specified Removed Surface"
#for i in range(NCB2.shape[0]):
  #neigh =KDTncB.query(NCB2[i,])[1]
  #RegTB = RegTB[np.where(RegTB<>neigh)[0]]
#print "Dump it to file"
#np.ma.dump(RegTB,'BSkullRemoveInternalNC2')

keepFA = np.ma.load('BSkullRemoveInternalNC2')
#NCB = np.ma.load('BSkull_NCmove')  
#NCB = np.ma.load('BSKull_Symm')
#NeigB = Neigenvalues(NCB,TriB,3,50)
#np.ma.dump(NeigB,'BskullNeig50')
#NeigB = np.ma.load('BskullFNeig3')
#keep = np.where(NeigB[:,0]<50*NeigB[:,1])[0]
#keep = np.array(range(0,NCB.shape[0]))
#NCB = 10*NCB
##NCT = 10*NCT
#NCB[:,0] = NCB[:,0]-np.average(NCB[:,0])
#NCB[:,1] = NCB[:,1]-np.average(NCB[:,1])
#NCB[:,2] = NCB[:,2]-np.average(NCB[:,2])
#np.ma.dump(NCB,'BSkull_NCmove')

#TriT=np.c_[TriB[:,0],TriB[:,2],TriB[:,1]]
#NCT = np.c_[-NCB[:,0],NCB[:,1:3]]

#ridgeN,valleyN,Kmax,Kmin,eMax,eMin,PdMax,PdMin = crestlines(np.array(range(NCB.shape[0])),NCB,TriB,3,200,1)	#(nodes,NC,Tri,NN,ls,neighSpat,eOrk)
##np.ma.dump(ridgeN,'BSkullSymmRNodes3NN')
##np.ma.dump(valleyN,'BSkullSymmVNodes3NN')
#np.ma.dump(Kmax,'BSkullSymmKmax3NN')
#np.ma.dump(Kmin,'BSkullSymmKmin3NN')
#np.ma.dump(eMax,'BSkullSymmEmax3NN')
#np.ma.dump(eMin,'BSkullSymmEmin3NN')
#np.ma.dump(PdMax,'BSkullSymmPDmax3NN')
#np.ma.dump(PdMin,'BSkullSymmPDmin3NN')
##Skull5:
#np.ma.dump(ridgeN,'Skull5RNodes100')
#np.ma.dump(valleyN,'Skull5VNodes100')
#np.ma.dump(Kmax,'Skull5Kmax100')
#np.ma.dump(Kmin,'Skull5Kmin100')
#np.ma.dump(eMax,'Skull5Emax100')
#np.ma.dump(eMin,'Skull5Emin100')
#np.ma.dump(PdMax,'Skull5PDmax100')
#np.ma.dump(PdMin,'Skull5PDmin100')

Kmax = np.ma.load('BSkullSymmKmax3NN')
Kmin = np.ma.load('BSkullSymmKmin3NN')
eMax = np.ma.load('BSkullSymmEmax3NN')
eMin = np.ma.load('BSkullSymmEmin3NN')
PdMax = np.ma.load('BSkullSymmPDmax3NN')
PdMin = np.ma.load('BSkullSymmPDmin3NN')
#ridgeN = np.ma.load('BSkullSymmRNodes200')
#valleyN = np.ma.load('BSkullSymmVNodes200')
KmaxT = np.ma.load('Skull5Kmax500')
KminT = np.ma.load('Skull5Kmin500')
#eMax = np.ma.load('Skull5Emax100')
#eMin = np.ma.load('Skull5Emin100')
#PdMax = np.ma.load('Skull5PDmax100')
#PdMin = np.ma.load('Skull5PDmin100')
#ridgeN = np.ma.load('Skull5RNodes100')
#valleyN = np.ma.load('Skull5VNodes100')

#FP = sh.FeatPoints(Kmax,Kmin,NCB,TriB,10,0.1,0.1)
#print NCB.shape,FP.shape
#np.ma.dump(FP,'FeatPointsBase10rad0101NN')
FP = np.ma.load('FeatPointsBase10rad0101NN')
FP2 = np.ma.load('FeatPointsBase10rad0202')
#print 'Valley Nodes'
#valleyN = ridgenodes(NCB,TriB,PdMax,PdMin,Kmax,Kmin,eMax,eMin,0,1)	#(NC,Tri,PdMax,PdMin,Kmax,Kmin,eMax,eMin,RidgeVal,eOrk)
#np.ma.dump(valleyN,'BSkullSymmVNodes200E')
#print 'Ridge Nodes'
#ridgeN = ridgenodes(NCB,TriB,PdMax,PdMin,Kmax,Kmin,eMax,eMin,1,1)
#np.ma.dump(ridgeN,'BSkullSymmRNodes200E')
#keep = np.where((np.abs(Kmin)>0.1)|(Kmax>0.1))[0]

#valleyN = valleyN[np.where(np.abs(Kmin[valleyN])>0.18)[0]]
#ridgeN = ridgeN[np.where(Kmax[ridgeN]>0.18)[0]]
#print '############ Valley Lines ##################'
#valleyLN = LineDraw(valleyN,NCB,TriB,PdMax,PdMin,-Kmin)
#np.ma.dump(valleyLN,'BSkullVLines100K')
#np.ma.dump(valleyLN,'BSkullSymmVLines200')
##print
#print '############ Ridge Lines ##################'
#ridgeLN = LineDraw(ridgeN,NCB,TriB,PdMin,PdMax,Kmax)
#np.ma.dump(ridgeLN,'BSkullRLines100K')
#np.ma.dump(ridgeLN,'BSkullSymmRLines200')
#ridgeLN = np.ma.load('BSkullRLines100K')
#valleyLN = np.ma.load('BSkullVLines100K')
ridgeLN = np.ma.load('BSkullSymmRLines200')
valleyLN = np.ma.load('BSkullSymmVLines200')
#print ridgeLN[0],valleyLN[0]

#print '############# Threshold ##############'
#valleyLNT = ThreshLines(valleyLN,3,NCB,-Kmin)
#np.ma.dump(valleyLNT,'BSkullSymmVLines200T1')
##np.ma.dump(valleyLNT,'BSkullVLines100KT50')
#ridgeLNT = ThreshLines(ridgeLN,3,NCB,Kmax)
#np.ma.dump(ridgeLNT,'BSkullSymmRLines200T1')
#print ridgeLNT[0],valleyLNT[0]
#np.ma.dump(ridgeLNT,'BSkullRLines100KT50')
ridgeLN1 = np.ma.load('BSkullSymmRLines200T1')
valleyLN1 = np.ma.load('BSkullSymmVLines200T1')
#RL_keepFA,VL_keepFA = [],[]
#noR,noV = 0,0
#print 'find ridge and valley lines that fall within user specified area'
#for i in range(1,valleyLN1[0]+1):
  #if find_repeats(np.r_[valleyLN1[i],keepFA])[0].size>0:
    #noV=noV+1
    #VL_keepFA = VL_keepFA+[valleyLN1[i]]
#VL_keepFA = [noV]+VL_keepFA
#print 'valleys now ',noV,' in stead of ',valleyLN1[0]
#for i in range(1,ridgeLN1[0]+1):
  #if find_repeats(np.r_[ridgeLN1[i],keepFA])[0].size>0:
    #noR=noR+1
    #RL_keepFA = RL_keepFA+[ridgeLN1[i]]
#RL_keepFA = [noR]+RL_keepFA
#print 'ridges now ',noR,' in stead of ',ridgeLN1[0]
#np.ma.dump(VL_keepFA,'BSkullVLineskeepFA2')
#np.ma.dump(RL_keepFA,'BSkullRLineskeepFA2')
RL_keepFA,VL_keepFA = np.ma.load('BSkullRLineskeepFA2'),np.ma.load('BSkullVLineskeepFA2')
ridgeLNT , valleyLNT = RL_keepFA,VL_keepFA
#ridgeLN2 , valleyLN2 = ridgeLN1 , valleyLN1
#print ridgeLN[0],valleyLN[0]
## Rotate scale and translate
#NCrst = lineRST(NCB,ridgeLN1,valleyLN1,NCT,ridgeLN2,valleyLN2,1,1)	#(NCB,RlinesB,VlinesB,NCT,RlinesT,VlinesT,UseFeat,UseScale)
#NCrst = (NCrst+NCB)/2
#np.ma.dump(NCrst,'Skull5_NCmoveRST')
NCB2 = np.ma.load('BSkull_NCmoveRST')
#TriT=np.c_[TriB[:,0],TriB[:,2],TriB[:,1]]
NCTsymm = np.c_[-NCB[:,0],NCB[:,1:3]]

#NCS,keepB,keepT = LineReg(NCB,TriB,ridgeLN1,valleyLN1,NCT,TriT,ridgeLN1,valleyLN1,80,10) #(NCB,TriB,RlinesB,VlinesB,NCT,TriT,RlinesT,VlinesT,percentReg)
#np.ma.dump(NCS,'BSkullto5_deform')
#np.ma.dump(keepB,'BSkullto5_keep')
#np.ma.dump(keepT,'Skull5fromB_keep')
#NCS = np.ma.load('BSkull_deform')
#keepB = np.ma.load('BSkull_keep')
#keepT = np.ma.load('BSkull_Skeep')

#FeatNodesB = np.where((Kmax>2)|(np.abs(Kmin)>2))[0]
#FeatLessNodesB = np.where((Kmax<2)|(np.abs(Kmin)<2))[0]
#FeatNodesT,FeatLessNodesT = FeatNodesB,FeatLessNodesB

#print 'Base nodes and triangles allowed'
#UseN_B,AllowableB = AllowableNaT(NCB,TriB,FeatNodesB,FeatLessNodesB,ridgeLN1,valleyLN1,keepB)
#print 'Target nodes and triangles allowed'
#UseN_T,AllowableT = AllowableNaT(NCT,TriT,FeatNodesT,FeatLessNodesT,ridgeLN1,valleyLN1,keepT)
#BFkeep = UseN_B[np.where((Kmax[UseN_B]>2)|(np.abs(Kmin[UseN_B])>2))[0]]
#np.ma.dump(UseN_B,'BSkull_useN')
#np.ma.dump(AllowableB,'BSkull_AlloweTri')
#UseN_T = np.ma.load('BSkull_useN')
#AllowableT = np.ma.load('BSkull_AlloweTri')
#UseN_B,AllowableB = UseN_T,AllowableT

####LandmB = ridge and valley nodes retained and LandmB_NC the associated final coordinates of these landmarks
#LandmB,LandmB_NC = keepB,NCS[keepB,]
#NCrbf = RBFmorph(NCB,NCB[LandmB,],LandmB_NC-NCB[LandmB,])
#np.ma.dump(NCrbf,'BSkull_RBFdeform')
#NCS = np.ma.load('BSkull_RBFdeform')
#AllowableB,AllowableT = np.array(range(0,TriB.shape[0])),np.array(range(0,TriT.shape[0]))
#UseN_B,UseN_T = np.array(range(0,NCB.shape[0])),np.array(range(0,NCT.shape[0]))
#USENORMALS = np.array([1,2,6,7,11,12,21,22,31,32])	# Iterations where normal information is used
#NCdef,Conv = elasticsurf(NCB/10,TriB,LandmB,LandmB_NC/10,AllowableB,NCT/10,TriT,AllowableT,UseN_B,UseN_T,100,USENORMALS)
#np.ma.dump(NCdef*10,'BSkull_Fulldeform')
NCS = np.ma.load('BSkull_Fulldeform')
#NCB=(NCB+NCS)/2
#np.ma.dump(NCB,'BSKull_Symm')
NCB = np.ma.load('BSKull_Symm')

#######################################
#######################################
#######################################
NCT,TriT = readmsh('Skull5.msh')
#DoReg5 = readmsh('Skull5keepRF.msh')[0]
#KDTNCB = KDTree(NCT,5)
#Reg5 = []
#for i in range(0,DoReg5.shape[0]):
  #neigh =KDTNCB.query(DoReg5[i,])[1]
  #Reg5 = Reg5 + [neigh]
#Reg5.sort()
#np.ma.dump(Reg5,'Skull5_NCallow')
Reg5 = np.ma.load('Skull5_NCallow')

NCT = np.ma.load('Skull5_NCmoveRST')
#ridgeLN1 = np.ma.load('BSkullRLines100KT50')
#valleyLN1 = np.ma.load('BSkullVLines100KT50')
ridgeLN2 = np.ma.load('Skull5RLines100T1')
valleyLN2 = np.ma.load('Skull5VLines100T1')



## Only keep registered nodes in outer area:
#RegB = np.ma.load('BSkullRemoveInternalNC')

#NCS,keepB,keepT = LineReg(NCB,TriB,ridgeLNT,valleyLNT,NCT,TriT,ridgeLN2,valleyLN2,50,5) #(NCB,TriB,RlinesB,VlinesB,NCT,TriT,RlinesT,VlinesT,percentReg,DistMax,RegB)
#np.ma.dump(NCS,'BSkullto5_deformRI2_2')
#np.ma.dump(keepB,'BSkullto5_keepRI2_2')
#np.ma.dump(keepT,'Skull5fromB_keepRI2_2')
NCS = np.ma.load('BSkullto5_deformRI2_2')
keepT = np.ma.load('Skull5fromB_keepRI2_2')
keepB = np.ma.load('BSkullto5_keepRI2_2')
#keepB = np.array([])
#for i in keepBI:
  #if np.where(keepFA==i)[0].size>0:
    #keepB=np.r_[keepB,i]
#keepB = np.array(keepB,int)
#print keepBI.shape, find_repeats(keepBI)
#print keepB.shape, find_repeats(keepB)
FeatNodesB = keepFA[np.where((Kmax[keepFA]>0.2)|(np.abs(Kmin[keepFA])>0.2))[0]]
#FeatLessNodesB = np.where((Kmax<0.08)&(np.abs(Kmin)<0.08))[0]
FeatNodesT = np.where((KmaxT>0.18)|(np.abs(KminT)>0.18))[0]
#FeatNodesT5 = Reg5[np.where((KmaxT[Reg5]>0.08)|(np.abs(KminT[Reg5])>0.08))[0]]
#FeatLessNodesT = np.where((KmaxT<0.08)&(np.abs(KminT)<0.08))[0]
#print 'Base nodes and triangles allowed'
#UseN_B,RemoveB,AllowableB = AllowableNaT(NCB,TriB,FeatNodesB,FeatLessNodesB,ridgeLN1,valleyLN1,keepB)
#np.ma.dump(UseN_B,'BSkullto5_useNRI2_2')
#np.ma.dump(AllowableB,'BSkullto5_AlloweTriRI2_2')
#UseN_B = np.array(range(NCB.shape[0]))
UseN_B = np.ma.load('BSkullto5_useNRI2_2')
#AllowableB = np.ma.load('BSkullto5_AlloweTriRI')
AllowableB = np.ma.load('BSkull_AlloweTriRI')
#FeatNodesB = UseN_B[np.where((Kmax[UseN_B]>0.1)|(np.abs(Kmin[UseN_B])>0.1))[0]]
#print 'Target nodes and triangles allowed'
#UseN_T,RemoveT,AllowableT = AllowableNaT(NCT,TriT,FeatNodesT5,FeatLessNodesT,ridgeLN2,valleyLN2,keepT)
#np.ma.dump(UseN_T,'Skull5fromB_useNRI2_2')
#np.ma.dump(RemoveT,'Skull5fromB_removeNRI2_2')
#np.ma.dump(AllowableT,'Skull5fromB_AlloweTriRI2_2')
UseN_T = np.ma.load('Skull5fromB_useNRI2_2')
#RemoveT = np.ma.load('Skull5fromB_removeNRI2')
AllowableT = np.ma.load('Skull5fromB_AlloweTriRI2_2')
FeatNodesT5 = UseN_T[np.where((KmaxT[UseN_T]>0.18)|(np.abs(KminT[UseN_T])>0.18))[0]]
#FeatNodesT =np.where((KmaxT>0.18)|(np.abs(KminT)>0.18))[0]

#UseN_B =  UseN_B[range(UseN_B.size/10)*10]
##UseN_T =  UseN_T[range(UseN_T.size/100)*100]
#NCB = np.ma.load('BSkullto5_FulldeformRI2_GWpND')
##NCS = np.ma.load('TempElasNodes_Iter10_TimeMonSep20')*10
###### Load tetrahedral mesh information
NCTnc = np.ma.load('BSkullFTnc')
TetT = np.ma.load('BSkullFTtet')
TriTet = np.ma.load('BSkullFTtri')
outer = np.ma.load('BSkullFTn_outer')
inner = np.ma.load('BSkullFTn_inner')
neighbList = np.ma.load('BSkullFTneighbTri')
neighbListTet = np.ma.load('BSkullFTneighbTet')
#LandmB,LandmB_NC = keepB,NCS[keepB,]
#USENORMALS = np.array([1,2,3,11,12,21,22])	# Iterations where normal information is used
#NCdef,Conv = ptet.elasticsurf(NCTnc,TetT,NCB,TriB,LandmB,LandmB_NC,AllowableB,NCT,TriT,AllowableT,UseN_B,UseN_T)
#np.ma.dump(NCdef,'BSkullto5_FulldeformRI2_TetSmooth')
#NCS1 = RBFmorph(NCB,NCB[LandmB],LandmB_NC-NCB[LandmB])
#np.ma.dump(NCS1,'BSkullto5_TPSm10')
#NCS = np.ma.load('TempElasNodes_Iter0_TimeMonOct4')
#NCS = np.ma.load('BSkullto5_FulldeformRI2_GW')
NCB,TriB = readmsh('BSkullSymm.msh')
NCS2 = np.ma.load('Bto5RegTemp40_21Octp20TaubinA')
NCS1 = np.ma.load('TempElasNodes_Iter40_TimeThuOct21')
#NCS2 = np.ma.load('Skull5TetSymm20T')
#NCS1 = np.c_[-NCS2[:,0],NCS2[:,1:]]

#EQ = np.ma.load('ElQual1')
#QMAT = np.ma.load('QMAT1')
#NCsmooth = ptet.LaplacMesh(inner,NCS2,neighbListTet,100,1)
#np.ma.dump(NCsmooth,'Skull5TetSymm100intLapl')
##NCopt = qu.OptEQ(NCS2,TetT,inner)
#NCopt = qu.OptINT(inner,outer,NCsmooth,TetT,neighbListTet,30)
#np.ma.dump(NCopt,'Skull5TetSymm10Topt')
#Def = NCS - NCB
#Def = np.sqrt(np.sum(Def*Def,1))
#SI = 0.5 - np.arctan((Kmax+Kmin)/(Kmax-Kmin))/np.pi
#figure(fgcolor=(0,0,0),bgcolor=(1,1,1))
#triangular_mesh(NCB[:,0],NCB[:,2],NCB[:,1],TriB,scalars = SI[:,0],colormap='RdGy')
#points3d(NCB[FP,0],NCB[FP,2],NCB[FP,1],color=(0,0,1),scale_factor=2)


#####ridgeLN1,valleyLN1 = np.ma.load('BSkullRLineskeepFA2'),np.ma.load('BSkullVLineskeepFA2')
#####ridgeLN2 = np.ma.load('Skull5RLines100T1')
#####valleyLN2 = np.ma.load('Skull5VLines100T1')
#####NCS = np.ma.load('BSkullto5_deformRI2_2')
#####NC1,Tri1 = readmsh('BSkullSymm.msh')
#####NC2,Tri2 = readmsh('Skull5.msh')
#####NC2 = np.ma.load('Skull5_NCmoveRST')
#####keepT = np.ma.load('Skull5fromB_keepRI2_2')
#####keepB = np.ma.load('BSkullto5_keepRI2_2')

#####figure(fgcolor=(0,0,0),bgcolor=(1,1,1))
#####triangular_mesh(NC1[:,0],NC1[:,2],NC1[:,1],Tri1,color=(0.6,0.6,0.8),opacity=1,name='Base')
#####for i in range(1,ridgeLN1[0]+1):
  #####plot3d(NC1[ridgeLN1[i],0],NC1[ridgeLN1[i],2],NC1[ridgeLN1[i],1],color=(0,0,1),representation='wireframe',line_width=2)
#####for i in range(1,valleyLN1[0]+1):
  #####plot3d(NC1[valleyLN1[i],0],NC1[valleyLN1[i],2],NC1[valleyLN1[i],1],color=(0,0,1),representation='wireframe',line_width=2)

#####figure(fgcolor=(0,0,0),bgcolor=(1,1,1))
#####triangular_mesh(NC2[:,0],NC2[:,2],NC2[:,1],Tri2,color=(0.8,0.6,0.6),opacity=1,name='Target')
#####for i in range(1,ridgeLN1[0]+1):
  #####plot3d(NC1[ridgeLN1[i],0],NC1[ridgeLN1[i],2],NC1[ridgeLN1[i],1],color=(0,0,1),representation='wireframe',line_width=2)
#####for i in range(1,valleyLN1[0]+1):
  #####plot3d(NC1[valleyLN1[i],0],NC1[valleyLN1[i],2],NC1[valleyLN1[i],1],color=(0,0,1),representation='wireframe',line_width=2)

#####figure(fgcolor=(0,0,0),bgcolor=(1,1,1))
#####triangular_mesh(NC2[:,0],NC2[:,2],NC2[:,1],Tri2,color=(0.8,0.6,0.6),opacity=1,name='Target')
#####for i in range(1,ridgeLN1[0]+1):
  #####plot3d(NCS[ridgeLN1[i],0],NCS[ridgeLN1[i],2],NCS[ridgeLN1[i],1],color=(0,0,1),representation='wireframe',line_width=2)
#####for i in range(1,valleyLN1[0]+1):
  #####plot3d(NCS[valleyLN1[i],0],NCS[valleyLN1[i],2],NCS[valleyLN1[i],1],color=(0,0,1),representation='wireframe',line_width=2)

#############################################################3
#####figure(fgcolor=(0,0,0),bgcolor=(1,1,1))
#####triangular_mesh(NC2[:,0],NC2[:,2],NC2[:,1],Tri2,color=(0.8,0.6,0.6),opacity=0.2,name='TargetRST') 

#####for i in range(1,ridgeLN1[0]+1):
  #####if find_repeats(np.r_[ridgeLN1[i],keepB])[0].size>0:
    #####plot3d(NCS[ridgeLN1[i],0],NCS[ridgeLN1[i],2],NCS[ridgeLN1[i],1],color=(0,0,1),representation='wireframe',line_width=2)
#####for i in range(1,valleyLN1[0]+1):
  #####if find_repeats(np.r_[valleyLN1[i],keepB])[0].size>0:
    #####plot3d(NCS[valleyLN1[i],0],NCS[valleyLN1[i],2],NCS[valleyLN1[i],1],color=(0,0,1),representation='wireframe',line_width=2)
  
#####for i in range(1,ridgeLN2[0]+1):
  #####if (ridgeLN2[i].size>4)&(find_repeats(np.r_[ridgeLN2[i],keepT])[0].size>4):
    #####plot3d(NC2[ridgeLN2[i],0],NC2[ridgeLN2[i],2],NC2[ridgeLN2[i],1],color=(1,0,0),representation='wireframe',line_width=2)
#####for i in range(1,valleyLN2[0]+1):
  #####if (valleyLN2[i].size>4)&(find_repeats(np.r_[valleyLN2[i],keepT])[0].size>4):
    #####plot3d(NC2[valleyLN2[i],0],NC2[valleyLN2[i],2],NC2[valleyLN2[i],1],color=(1,0,0),representation='wireframe',line_width=2)

  
###############################################################3
######figure(fgcolor=(0,0,0),bgcolor=(1,1,1))
######triangular_mesh(NC2[:,0],NC2[:,2],NC2[:,1],Tri2,color=(0.8,0.6,0.6),opacity=0.2,name='TargetOrig')

######for i in range(1,ridgeLN2[0]+1):
  ######if ridgeLN2.size>4:
    ######plot3d(NC2[ridgeLN2[i],0],NC2[ridgeLN2[i],2],NC2[ridgeLN2[i],1],color=(1,0,0),representation='wireframe',line_width=2)
######for i in range(1,valleyLN2[0]+1):
  ######if valleyLN2.size>4:
    ######plot3d(NC2[valleyLN2[i],0],NC2[valleyLN2[i],2],NC2[valleyLN2[i],1],color=(1,0,0),representation='wireframe',line_width=2)

######for i in range(1,ridgeLN1[0]+1):
  ######plot3d(NC1[ridgeLN1[i],0],NC1[ridgeLN1[i],2],NC1[ridgeLN1[i],1],color=(0,0,1),representation='wireframe',line_width=2)
######for i in range(1,valleyLN1[0]+1):
  ######plot3d(NC1[valleyLN1[i],0],NC1[valleyLN1[i],2],NC1[valleyLN1[i],1],color=(0,0,1),representation='wireframe',line_width=2)





















#figure(fgcolor=(0,0,0),bgcolor=(1,1,1))
#triangular_mesh(NCT[:,0],NCT[:,2],NCT[:,1],TriT,color=(0.8,0.6,0.6),opacity=1)
#triangular_mesh(NCB[:,0],NCB[:,2],NCB[:,1],TriB,color=(0.6,0.6,0.8),opacity=1)
##triangular_mesh(NCS[:,0],NCS[:,2],NCS[:,1],TriB,color=(0.6,0.8,0.8),opacity=1)
###triangular_mesh(NCS[:,0],NCS[:,2],NCS[:,1],TriB,scalars = Def)
#triangular_mesh(NCS1[:,0],NCS1[:,2],NCS1[:,1],TriTet,color=(0.6,0.6,0.8),opacity=1)
##triangular_mesh(NCS2[:,0],NCS2[:,2],NCS2[:,1],TriTet,color=(0.6,0.6,0.8),opacity=1)

#figure(fgcolor=(0,0,0),bgcolor=(1,1,1))
#####points3d(NCB[keepFA,0],NCB[keepFA,2],NCB[keepFA,1],color=(0,0,0),scale_factor=0.5)
#triangular_mesh(NCT[:,0],NCT[:,2],NCT[:,1],TriT,color=(0.8,0.8,0.7),opacity=1)
#triangular_mesh(NCB[:,0],NCB[:,2],NCB[:,1],TriB,color=(0.8,0.8,0.7),opacity=1)
#points3d(NCT[FeatNodesT,0],NCT[FeatNodesT,2],NCT[FeatNodesT,1],color=(0,0,0),opacity=0.2,scale_factor=0.5)
#points3d(NCT[FeatNodesT5,0],NCT[FeatNodesT5,2],NCT[FeatNodesT5,1],color=(0,0,0),opacity=0.2,scale_factor=0.6)
####points3d(NCT[T5Fkeep,0],NCT[T5Fkeep,2],NCT[T5Fkeep,1],color=(1,0,0),scale_factor=0.5)
#points3d(NCB[FeatNodesB,0],NCB[FeatNodesB,2],NCB[FeatNodesB,1],color=(0,0,0),opacity=0.2,scale_factor=0.5)
####points3d(NCB[AllowableB,0],NCB[AllowableB,2],NCB[AllowableB,1],color=(0,0,1),scale_factor=0.8)
####points3d(NCB[BFkeep,0],NCB[BFkeep,2],NCB[BFkeep,1],color=(1,0,0),scale_factor=0.5)

#figure(fgcolor=(0,0,0),bgcolor=(1,1,1))
##points3d(NCT[keepT,0],NCT[keepT,2],NCT[keepT,1],color=(1,0,0),scale_factor=1)
##points3d(NCB[keepB,0],NCB[keepB,2],NCB[keepB,1],color=(0,0,1),scale_factor=1)
##points3d(NCS[keepB,0],NCS[keepB,2],NCS[keepB,1],color=(0,0,1),scale_factor=1)
###NCT=NCB
###figure(fgcolor=(0,0,0),bgcolor=(1,1,1))
##points3d(NCB[ridgeN,0],NCB[ridgeN,2],NCB[ridgeN,1],color=(1,0,0),scale_factor=1)
##points3d(NCB[valleyN,0],NCB[valleyN,2],NCB[valleyN,1],color=(0,0,1),scale_factor=1)
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
##for i in range(1,ridgeLN[0]+1):
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
##for i in range(1,ridgeLNT[0]+1):
  ##Th = lineThresh(ridgeLNT[i],NCS,Kmax)
  ##plot3d(NCS[ridgeLNT[i],0],NCS[ridgeLNT[i],2],NCS[ridgeLNT[i],1],color=(1,0,0),representation='wireframe',line_width=2)
##for i in range(1,valleyLNT[0]+1):
  ##Th = lineThresh(valleyLNT[i],NCS,-Kmin)
  ##plot3d(NCS[valleyLNT[i],0],NCS[valleyLNT[i],2],NCS[valleyLNT[i],1],color=(0,0,1),representation='wireframe',line_width=2)


#### To call FEBio:
#### import os
####os.system('~/Software/MRL/FEBio/febio.lnx -i ~/Documents/MASTERS/Skulls/Prognathic/PrognI.feb') for exemple

##NC0 = NCS
#NC1 = NCS1
#NC2 = NCS2
###NC0 = np.ma.load('BSkullto5_Fulldeform')
###NC1 = np.ma.load('BSkullto5_FulldeformRI2')
###NC0 = np.ma.load('BSkull_RBFdeform')
###NC1 = np.ma.load('TempBto5_Iter10')*10
###NC2 = np.ma.load('TempBto5_Iter20')*10
###NC3 = np.ma.load('TempBto5_Iter30')*10
###NC4 = np.ma.load('TempBto5_Iter40')*10
###NC5 = np.ma.load('TempBto5_Iter50')*10
###NC6 = np.ma.load('TempBto5_Iter60')*10
###NC7 = np.ma.load('TempBto5_Iter70')*10
###NC8 = np.ma.load('TempBto5_Iter80')*10
###NC9 = np.ma.load('TempBto5_Iter90')*10
###NC10 = np.ma.load('BSkullto5_Fulldeform')
###NC11 = np.ma.load('TempElasNodes_Iter11_TimeWedSep15')*10
###NC12 = np.ma.load('TempElasNodes_Iter12_TimeWedSep15')*10
###NC13 = np.ma.load('TempElasNodes_Iter13_TimeWedSep15')*10
#targeTF = triangular_mesh(NCT[:,0],NCT[:,2],NCT[:,1],TriT,color=(0.9,0.5,0.5),opacity=0)
#targeTA = triangular_mesh(NCT[:,0],NCT[:,2],NCT[:,1],TriT[AllowableT,],color=(0.9,0.5,0.5),opacity=0)
#basE = triangular_mesh(NCB[:,0],NCB[:,2],NCB[:,1],TriB,color=(0.5,0.9,0.5),opacity=0)
##plt0 = triangular_mesh(NC0[:,0],NC0[:,2],NC0[:,1],TriB,color=(0.5,0.9,0.5),opacity=0)
#plt1 = triangular_mesh(NC1[:,0],NC1[:,2],NC1[:,1],TriTet,color=(0.5,0.9,0.5),opacity=0)
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

cutpTF = pipeline.scalar_cut_plane(targeTF,color=(1,1,0),line_width=5,plane_orientation='z_axes',name='TargetF')
cutpTF.implicit_plane.origin = (0,0,-40)
cutpTF.implicit_plane.widget.enabled = False
cutpTA = pipeline.scalar_cut_plane(targeTA,color=(1,0,0),line_width=5,plane_orientation='z_axes',name='TargetA')
cutpTA.implicit_plane.origin = (0,0,-40)
cutpTA.implicit_plane.widget.enabled = False
cutpB = pipeline.scalar_cut_plane(basE,color=(0,0,0),line_width=3,plane_orientation='z_axes',name='Base')
cutpB.implicit_plane.origin = (0,0,-40)
cutpB.implicit_plane.widget.enabled = False
#cutp0 = pipeline.scalar_cut_plane(plt0,color=(0,0,1),line_width=3,plane_orientation='z_axes',name='iter_00')
#cutp0.implicit_plane.origin = (0,0,-40)#(20,0,0)
#cutp0.implicit_plane.widget.enabled = False
cutp1 = pipeline.scalar_cut_plane(plt1,color=(0,0,0),line_width=3,plane_orientation='z_axes',name='iter_01')
cutp1.implicit_plane.origin = (0,0,-40)#(20,0,0)
cutp1.implicit_plane.widget.enabled = False
cutp2 = pipeline.scalar_cut_plane(plt2,color=(0,0,0),line_width=3,plane_orientation='z_axes',name='iter_01')
cutp2.implicit_plane.origin = (0,0,-40)#(20,0,0)
cutp2.implicit_plane.widget.enabled = False
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