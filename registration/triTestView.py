# -*- coding: utf-8 -*-
import numpy as np                                                                                   
from enthought.mayavi.mlab import *
from fileread import *
from pmeshdeform import *
from scipy.io import *
from scipy.interpolate import Rbf
from scipy.spatial import KDTree
from scipy.stats.mstats import find_repeats
from ptrimeshalt import *
from ptrifeatures import *
from plinedeform import *
import time

#NCS,Tri,LMB = readmsh('Bishop.msh')
#NCS,Tri,LMB = readmsh('hand1.msh')
#NCS,Tri,LMB = readmsh('dolph1n.msh')
#NCS,Tri,LMB = readmsh('StarGeo.msh')
NC2,Tri2 = readmsh('dolph3n.msh')
NC1,Tri1 = readmsh('dolph1n.msh')
#NCS,Tri,LM = readmsh('skrewdriver.msh')
#NCS,Tri,LM1 = readmsh('Skull5.msh')
#NCS,Tri,LM1 = readmsh('BSkullFine.msh')
#NCS,Tri,LM1 = readmsh('tstar.msh')

LM1 = [5995,2365,9593,1448,8542,6271,13555,4554,11849,14787,6483,6705,12554,5420]
LM2 = [5257,883,6203,4187,6022,5332,7122,4972,3561,104,7777,1925,3645,6923]
knownC = NC1[LM1,]
knownD = NC2[LM2,]-knownC
NCrbf = RBFmorph(NC1,knownC,knownD)
np.ma.dump(NCrbf,'dolph1nRBF_TPS')
#LMB=np.array(LMB[:,0],int)
#LM2=np.array(LM2[:,0],int)
#triB = np.array(triB,int)
#tri2 = np.array(tri2,int)

#NeigB = Neigenvalues(NCS,Tri,3,50)
#np.ma.dump(NeigB,'dolph3nFNeig50NN')
#Neig = np.ma.load('dolph3nFNeig50NN')
#keep = np.where(Neig[:,0]<50*Neig[:,1])[0]
#print NCS.size/3, keep.size
#keep=np.array(range(0,NCS.size/3))

#for i in range(0,5):	# smooth before extraction of surface information
  #print 'smoothing iteration ',i+1
  #NC = LaplacianSmooth(np.array(range(0,NC.size/3)),NC,Tri)

#NCS = np.ma.load('BskullNC3LaplacianS')
#NCS = np.ma.load('dolph3nNC3LaplacianS')
#ridgeN,valleyN,Kmax,Kmin,eMax,eMin,PdMax,PdMin = crestlines(keep,NCS,Tri,3,100,0)

#np.ma.dump(ridgeN,'dolph3nFRNodes100T50')
#np.ma.dump(valleyN,'dolph3nFVNodes100T50')
#np.ma.dump(Kmax,'dolph3nFKmax100T50')
#np.ma.dump(Kmin,'dolph3nFKmin100T50')
#np.ma.dump(eMax,'dolph3nFEmax100T50')
#np.ma.dump(eMin,'dolph3nFEmin100T50')
#np.ma.dump(PdMax,'dolph3nFPDmax100T50')
#np.ma.dump(PdMin,'dolph3nFPDmin100T50')

#keep1 = np.where(Neig[:,0]>2000*Neig[:,1])[0]
#keep2 = np.where(Neig[:,0]>1000*Neig[:,1])[0]
#keep2 = np.where(Neig[:,0]>20000*Neig[:,1])[0]
#keep3 = np.where(Neig[:,0]<5*Neig[:,1])[0]
#keep5 = np.where(Neig[:,0]<20*Neig[:,1])[0]
#keep4 = np.where(Neig[:,0]<100*Neig[:,1])[0]

#Kmax = np.ma.load('dolph3nFKmax100T50')
#Kmin = np.ma.load('dolph3nFKmin100T50')
#eMax = np.ma.load('dolph3nFEmax100T50')
#eMin = np.ma.load('dolph3nFEmin100T50')
#PdMax = np.ma.load('dolph3nFPDmax100T50')
#PdMin = np.ma.load('dolph3nFPDmin100T50')
#ridgeN = np.ma.load('dolph3nFRNodes100T50')
#valleyN = np.ma.load('dolph3nFVNodes100T50')
#ridgeLN = np.ma.load('dolph3nFRLines100KT')
#valleyLN = np.ma.load('dolph3nFVLines100KT')

ridgeLN1 = np.ma.load('dolph1nFRLines100KT10T')
valleyLN1 = np.ma.load('dolph1nFVLines100KT10T')
ridgeLN2 = np.ma.load('dolph3nFRLines100KT10T')
valleyLN2 = np.ma.load('dolph3nFVLines100KT10T')

#NC2 = np.ma.load('dolph3nFRNC_FICP')

#NCS = lineRST(NC1,ridgeLN1,valleyLN1,NC2,ridgeLN2,valleyLN2,0,0,1)
#np.ma.dump(NCS,'dolph3nFRS1T')
#NC3 = np.ma.load('dolph3nFRS1T')
#NCS=NC3
NC3=NCrbf

#NCS,keepB,keepT = LineReg(NC1,Tri1,ridgeLN1,valleyLN1,NC3,Tri2,ridgeLN2,valleyLN2,50,20)
#np.ma.dump(NCS,'dolph1nFRS1T_deform')
#np.ma.dump(keepB,'dolph1nFRS1T_keep')
#np.ma.dump(keepT,'dolph3nFRS1T_keep')
NCS = np.ma.load('dolph1nFRS1T_deform')
keepB = np.ma.load('dolph1nFRS1T_keep')
keepT = np.ma.load('dolph3nFRS1T_keep')
#valleyN = ridgenodes(NCS,Tri,PdMax,PdMin,Kmax,Kmin,eMax,eMin,0,0)
#np.ma.dump(valleyN,'dolph3nFVNodes100T50')
#ridgeN = ridgenodes(NCS,Tri,PdMax,PdMin,Kmax,Kmin,eMax,eMin,1,0)
#np.ma.dump(ridgeN,'dolph3nFRNodes100T50')

#r2 = ridgeN[np.where(Neig[ridgeN,0]<50*Neig[ridgeN,1])[0]]
#v2 = valleyN[np.where(Neig[valleyN,0]<50*Neig[valleyN,1])[0]]

#print '############ Valley Lines ##################'
#valleyLN = LineDraw(v2,NCS,Tri,PdMax,PdMin,-Kmin)
#np.ma.dump(valleyLN,'dolph3nFVLines100KT')
#print
#print '############ Ridge Lines ##################'
#ridgeLN = LineDraw(r2,NCS,Tri,PdMin,PdMax,Kmax)
#np.ma.dump(ridgeLN,'dolph3nFRLines100KT')

#print '############# Threshold ##############'
#valleyLNT = ThreshLines(valleyLN,0.005,NCS,-Kmin)
#np.ma.dump(valleyLNT,'skrewdVLines5T0.005')
#ridgeLNT = ThreshLines(ridgeLN,0.005,NCS,Kmax)
#np.ma.dump(ridgeLNT,'skrewdRLines5T0.005')

#rn = np.where(Kmax>np.abs(Kmin))[0]
#vn = np.where(Kmin<-np.abs(Kmax))[0]

#NCS = 10*NCS
NCB,NCT,TriB,TriT = NC1,NC3,Tri1,Tri2
#LandmB,LandmB_NC = keepB,NCS[keepB,]
#AllowableB,AllowableT = np.array(range(0,TriB.shape[0])),np.array(range(0,TriT.shape[0]))
#UseN_B,UseN_T = np.array(range(0,NCB.shape[0]/500)*500),np.array(range(0,NCT.shape[0]/700)*700)
#USENORMALS = np.array([1])#1,2,3,4,5])	# Iterations where normal information is used
#NCdef,Conv = elasticsurf(NCB/10,TriB,LandmB,LandmB_NC/10,AllowableB,NCT/10,TriT,AllowableT,UseN_B,UseN_T,5,USENORMALS)
#np.ma.dump(NCdef,'dolph1n_FulldeformP')
#NCS = np.ma.load('dolph1n_FulldeformP')



figure(fgcolor=(0,0,0),bgcolor=(1,1,1))
triangular_mesh(NC2[:,0],NC2[:,2],NC2[:,1],Tri2,color=(0.8,0.6,0.6),opacity=1)
triangular_mesh(NC1[:,0],NC1[:,2],NC1[:,1],Tri1,color=(0.6,0.6,0.8),opacity=1)
triangular_mesh(NC3[:,0],NC3[:,2],NC3[:,1],Tri1,color=(0.8,0.6,0.6),opacity=1)

#figure(fgcolor=(0,0,0),bgcolor=(1,1,1))
#points3d(NC2[keepT,0],NC2[keepT,2],NC2[keepT,1],color=(1,0,0),scale_factor=0.4)
#points3d(NC1[keepB,0],NC1[keepB,2],NC1[keepB,1],color=(0,0,1),scale_factor=0.4)
#points3d(NC3[keepT,0],NC3[keepT,2],NC3[keepT,1],color=(1,0,0),scale_factor=0.4)
#points3d(NCS[keepB,0],NCS[keepB,2],NCS[keepB,1],color=(0,0,1),scale_factor=0.4)
##points3d(NCS[keep,0],NCS[keep,2],NCS[keep,1],color=(0,0,0),scale_factor=0.2)
#points3d(NCS[ridgeN,0],NCS[ridgeN,2],NCS[ridgeN,1],color=(1,0,0),scale_factor=0.5)
#points3d(NCS[valleyN,0],NCS[valleyN,2],NCS[valleyN,1],color=(0,0,1),scale_factor=0.5)
#points3d(NCS[rn,0],NCS[rn,2],NCS[rn,1],color=(1,0,0),scale_factor=0.1,opacity=0.2)
#points3d(NCS[vn,0],NCS[vn,2],NCS[vn,1],color=(0,0,1),scale_factor=0.1,opacity=0.2)
#points3d(NCS[r2,0],NCS[r2,2],NCS[r2,1],color=(1,0,0),scale_factor=0.4)#,opacity=0.2)
#points3d(NCS[v2,0],NCS[v2,2],NCS[v2,1],color=(0,0,1),scale_factor=0.4)#,opacity=0.2)
#points3d(NCS[keep1,0],NCS[keep1,2],NCS[keep1,1],color=(0,0,0),scale_factor=0.1,opacity=0.2)
#points3d(NCS[keep2,0],NCS[keep2,2],NCS[keep2,1],color=(0,0,0),scale_factor=0.1,opacity=0.2)
#points3d(NCS[keep3,0],NCS[keep3,2],NCS[keep3,1],color=(0,0,0),scale_factor=0.1,opacity=0.2)
#points3d(NCS[keep4,0],NCS[keep4,2],NCS[keep4,1],color=(0,0,0),scale_factor=0.1,opacity=0.2)
#points3d(NCS[keep5,0],NCS[keep5,2],NCS[keep5,1],color=(0,0,0),scale_factor=0.1,opacity=0.2)
#points3d(NCS[keep6,0],NCS[keep6,2],NCS[keep6,1],color=(0,0,0),scale_factor=0.1,opacity=0.2)
#KA = 1000*np.average(Kmax[ridgeN])
#KB = 10000*np.average(-Kmin[valleyN])
#count1 = 0
#count2 = 0
#figure(fgcolor=(0,0,0),bgcolor=(1,1,1))

#for i in range(1,ridgeLN2[0]+1):
  #plot3d(NC2[ridgeLN2[i],0],NC2[ridgeLN2[i],2],NC2[ridgeLN2[i],1],color=(1,0,0),representation='wireframe',line_width=2)
#for i in range(1,valleyLN2[0]+1):
  #plot3d(NC2[valleyLN2[i],0],NC2[valleyLN2[i],2],NC2[valleyLN2[i],1],color=(1,0,0),representation='wireframe',line_width=2)
#triangular_mesh(NC2[:,0],NC2[:,2],NC2[:,1],Tri2,color=(0.8,0.6,0.6),opacity=0.2,name='TargetOrig')

#for i in range(1,ridgeLN1[0]+1):
  #plot3d(NC1[ridgeLN1[i],0],NC1[ridgeLN1[i],2],NC1[ridgeLN1[i],1],color=(0,0,1),representation='wireframe',line_width=2)
#for i in range(1,valleyLN1[0]+1):
  #plot3d(NC1[valleyLN1[i],0],NC1[valleyLN1[i],2],NC1[valleyLN1[i],1],color=(0,0,1),representation='wireframe',line_width=2)

####################################################3
#figure(fgcolor=(0,0,0),bgcolor=(1,1,1))

#for i in range(1,ridgeLN2[0]+1):
  #plot3d(NC3[ridgeLN2[i],0],NC3[ridgeLN2[i],2],NC3[ridgeLN2[i],1],color=(1,0,0),representation='wireframe',line_width=2)
#for i in range(1,valleyLN2[0]+1):
  #plot3d(NC3[valleyLN2[i],0],NC3[valleyLN2[i],2],NC3[valleyLN2[i],1],color=(1,0,0),representation='wireframe',line_width=2)
#triangular_mesh(NC3[:,0],NC3[:,2],NC3[:,1],Tri2,color=(0.8,0.6,0.6),opacity=0.2,name='TargetRST')  

#for i in range(1,ridgeLN1[0]+1):
  #plot3d(NC1[ridgeLN1[i],0],NC1[ridgeLN1[i],2],NC1[ridgeLN1[i],1],color=(0,0,1),representation='wireframe',line_width=2)
#for i in range(1,valleyLN1[0]+1):
  #plot3d(NC1[valleyLN1[i],0],NC1[valleyLN1[i],2],NC1[valleyLN1[i],1],color=(0,0,1),representation='wireframe',line_width=2)


########################################################3
#figure(fgcolor=(0,0,0),bgcolor=(1,1,1))

#for i in range(1,ridgeLN1[0]+1):
  #if find_repeats(np.r_[ridgeLN1[i],keepB])[0].size>0:
    #plot3d(NCS[ridgeLN1[i],0],NCS[ridgeLN1[i],2],NCS[ridgeLN1[i],1],color=(0,0,1),representation='wireframe',line_width=2)
#for i in range(1,valleyLN1[0]+1):
  #if find_repeats(np.r_[valleyLN1[i],keepB])[0].size>0:
    #plot3d(NCS[valleyLN1[i],0],NCS[valleyLN1[i],2],NCS[valleyLN1[i],1],color=(0,0,1),representation='wireframe',line_width=2)
    
#for i in range(1,ridgeLN2[0]+1):
  #if find_repeats(np.r_[ridgeLN2[i],keepT])[0].size>0:
    #plot3d(NC3[ridgeLN2[i],0],NC3[ridgeLN2[i],2],NC3[ridgeLN2[i],1],color=(1,0,0),representation='wireframe',line_width=2)
#for i in range(1,valleyLN2[0]+1):
  #if find_repeats(np.r_[valleyLN2[i],keepT])[0].size>0:
    #plot3d(NC3[valleyLN2[i],0],NC3[valleyLN2[i],2],NC3[valleyLN2[i],1],color=(1,0,0),representation='wireframe',line_width=2)
#triangular_mesh(NC3[:,0],NC3[:,2],NC3[:,1],Tri2,color=(0.8,0.6,0.6),opacity=0.2,name='TargetRST') 
  
#for i in range(1,ridgeLN[0]+1):
  ##Thick = np.max(Kmax[ridgeLN[i][0]])
  #Th = lineThresh(ridgeLN[i][0],NCS,Kmax)
  ##print Th
  ###Thick = Th/ridgeLN[i][0].size
  ###if Th>10:
    ###count1=count1+1
  #if (ridgeLN[i][0].size>3)&(Th>10)&(Th<100):
    #plot3d(NCS[ridgeLN[i][0],0],NCS[ridgeLN[i][0],2],NCS[ridgeLN[i][0],1],color=(1,0,0),representation='wireframe',line_width=1)
  #if (ridgeLN[i][0].size>3)&(Th>100)&(Th<800):
    #plot3d(NCS[ridgeLN[i][0],0],NCS[ridgeLN[i][0],2],NCS[ridgeLN[i][0],1],color=(1,0,0),representation='wireframe',line_width=2)#Thick/KA)
  #if (ridgeLN[i][0].size>3)&(Th>800):
    #plot3d(NCS[ridgeLN[i][0],0],NCS[ridgeLN[i][0],2],NCS[ridgeLN[i][0],1],color=(1,0,0),representation='wireframe',line_width=4)#Thick/KA)
#for i in range(1,valleyLN[0]+1):
  ###Thick = np.max(-Kmin[valleyLN[i][0]])
  #Th = lineThresh(valleyLN[i][0],NCS,-Kmin)
  ##print Th
  ##if Th>10:
    ##count2=count2+1
  ##Thick = Th/valleyLN[i][0].size
  #if (valleyLN[i][0].size>3)&(Th>10)&(Th<100):
    #plot3d(NCS[valleyLN[i][0],0],NCS[valleyLN[i][0],2],NCS[valleyLN[i][0],1],color=(1,0,0),representation='wireframe',line_width=1)
  #if (valleyLN[i][0].size>3)&(Th>100)&(Th<800):
    #plot3d(NCS[valleyLN[i][0],0],NCS[valleyLN[i][0],2],NCS[valleyLN[i][0],1],color=(1,0,0),representation='wireframe',line_width=2)
  #if (valleyLN[i][0].size>3)&(Th>800):
    #plot3d(NCS[valleyLN[i][0],0],NCS[valleyLN[i][0],2],NCS[valleyLN[i][0],1],color=(1,0,0),representation='wireframe',line_width=4)#Thick/KB)
#print count1,count2
# Thresholded lines
#for i in range(1,ridgeLNT[0]+1):
  #Th = lineThresh(ridgeLNT[i],NCS,Kmax)
  #plot3d(NCS[ridgeLNT[i],0],NCS[ridgeLNT[i],2],NCS[ridgeLNT[i],1],color=(1,0,0),representation='wireframe',line_width=Th/20)
#for i in range(1,valleyLNT[0]+1):
  #Th = lineThresh(valleyLNT[i],NCS,-Kmin)
  #plot3d(NCS[valleyLNT[i],0],NCS[valleyLNT[i],2],NCS[valleyLNT[i],1],color=(0,0,1),representation='wireframe',line_width=Th/20)

#keep1 = np.ma.load('KeepN_5_20_100.dolph1n.txt')
#keep2 = np.ma.load('KeepN_5_20_1000_dolph3n.txt')
#NCED = np.ma.load('ElasticNodes7May2.txt')

#lmkeep=np.where((NC2[LM2[:,],0]<1)|(NC2[LM2[:,],2]>-50))[0]
#LMB=LMB[lmkeep,]
#LM2=LM2[lmkeep,]
#node4SFP = np.where((NC2[:,0]<0)|(NC2[:,2]>-50))[0]
#nodes = np.where((NC2[keep2,0]<0)|(NC2[keep2,0]>-50))[0]
#keep2=keep2[nodes,]
#k=0
#SFPTT=np.zeros((tri2.size/3,))
#for i in range(0,tri2.size/3):
  #if find_repeats(np.r_[tri2[i,].T,node4SFP])[0].size==3:
    #SFPTT[k,]=i
    #k=k+1
#SFPTT=SFPTT[range(0,k),]
#SFPTT=np.array(SFPTT,int)

#knownC = NCB[LMB,]
#knownD = NC2[LM2,]-knownC
#rbfx = Rbf(knownC[:,0],knownC[:,1],knownC[:,2],knownD[:,0])
#rbfy = Rbf(knownC[:,0],knownC[:,1],knownC[:,2],knownD[:,1])
#rbfz = Rbf(knownC[:,0],knownC[:,1],knownC[:,2],knownD[:,2])
#DSx = rbfx(NCB[:,0],NCB[:,1],NCB[:,2])
#DSy = rbfy(NCB[:,0],NCB[:,1],NCB[:,2])
#DSz = rbfz(NCB[:,0],NCB[:,1],NCB[:,2])
#NCI = NCB + np.c_[DSx,DSy,DSz]



#figure(1)
##triangular_mesh(NCB[:,0],NCB[:,1],NCB[:,2],triB,color=(1,1,1),representation='wireframe',opacity=0.2)
#triangular_mesh(NCB[:,0],NCB[:,2],NCB[:,1],triB,color=(0,0,1),representation='wireframe',opacity=0.2)
#triangular_mesh(NCI[:,0],NCI[:,2],NCI[:,1],triB,color=(0,1,1),representation='wireframe',opacity=0.2)
#triangular_mesh(NCED[:,0],NCED[:,2],NCED[:,1]-40,triB,color=(0,1,0),representation='wireframe',opacity=0.2)
#points3d(NCB[LMB,0],NCB[LMB,2],NCB[LMB,1]+40,color=(1,0,0),scale_factor=2)
#points3d(NCI[LMB,0],NCI[LMB,2],NCI[LMB,1],color=(1,0,0),scale_factor=2)
#points3d(NCED[LMB,0],NCED[LMB,2],NCED[LMB,1]-40,color=(1,0,0),scale_factor=2)

#figure(2)
#triangular_mesh(NCED[:,0],NCED[:,2],NCED[:,1],triB,color=(0,1,0),representation='wireframe',opacity=0.2)
#triangular_mesh(NC2[:,0],NC2[:,2],NC2[:,1],tri2,color=(0,1,1),opacity=0.3)
#triangular_mesh(NC2[:,0],NC2[:,2],NC2[:,1],tri2[SFPTT,],color=(0,0,0),opacity=0.3)
#points3d(NCED[LMB,0],NCED[LMB,2],NCED[LMB,1],color=(1,0,0),scale_factor=2)
#points3d(NC2[LM2,0],NC2[LM2,2],NC2[LM2,1],color=(0,0,0),scale_factor=2)
