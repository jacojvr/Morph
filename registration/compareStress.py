# -*- coding: utf-8 -*-
import numpy as np
import pyvtk as pv

### Von Mises of average stress tensor of P and O skulls 
fname = 'SkullPO_AveragePO'
S = (np.ma.load('Skull100Stress')[:,0:6] + np.ma.load('Skull000Stress')[:,0:6])/2
NCS = (np.ma.load('SkullPO_NC100')+np.ma.load('SkullPO_NC000'))/2
TetT = np.ma.load('SkullPO_TetUse')
VonM = np.sqrt((np.power(S[:,0]-S[:,1],2)+np.power(S[:,1]-S[:,2],2)+np.power(S[:,0]-S[:,2],2)+6*(S[:,3]*S[:,3]+S[:,4]*S[:,4]+S[:,5]*S[:,5]))/2)
skvtk = pv.VtkData(pv.UnstructuredGrid(points = NCS,tetra=TetT),'skull unique mesh',pv.CellData(pv.Scalars(VonM,'VonMises')))
skvtk.tofile(fname+'.vtk')
np.ma.dump(VonM,fname+'VonMises')

### Von Mises of Difference in stress tensor for Average and Average between P and O
fname = 'SkullPO_AvPOminAverage'
S = S - np.ma.load('Skull050Stress')[:,0:6]
NCS = np.ma.load('SkullPO_NC050')
TetT = np.ma.load('SkullPO_TetUse')
VonM = np.sqrt((np.power(S[:,0]-S[:,1],2)+np.power(S[:,1]-S[:,2],2)+np.power(S[:,0]-S[:,2],2)+6*(S[:,3]*S[:,3]+S[:,4]*S[:,4]+S[:,5]*S[:,5]))/2)
skvtk = pv.VtkData(pv.UnstructuredGrid(points = NCS,tetra=TetT),'skull unique mesh',pv.CellData(pv.Scalars(VonM,'VonMises')))
skvtk.tofile(fname+'.vtk')
np.ma.dump(VonM,fname+'VonMises')

###detS = S[:,0]*(S[:,1]*S[:,2]-S[:,4]*S[:,4]) - S[:,3]*(S[:,3]*S[:,2]-S[:,5]*S[:,4]) + S[:,5]*(S[:,3]*S[:,4]-S[:,5]*S[:,1])

###i=0
###Stensor = np.array([[S[i,0],S[i,3],S[i,5]],[S[i,3],S[i,1],S[i,4]],[S[i,5],S[i,2],S[i,4]]])
###Pstresses = np.linalg.svd(Stensor)[1]

###### Compare Von Mises of difference in stress tensor for 3 Prognathic and 3 orthognathic skulls
###fname = 'SkullPO_00_00'
###S = np.ma.load('Skull100Stress')[:,0:6] - np.ma.load('Skull000Stress')[:,0:6]
###NCS = (np.ma.load('SkullPO_NC100')+np.ma.load('SkullPO_NC000'))/2
###TetT = np.ma.load('SkullPO_TetUse')
###VonM = np.sqrt((np.power(S[:,0]-S[:,1],2)+np.power(S[:,1]-S[:,2],2)+np.power(S[:,0]-S[:,2],2)+6*(S[:,3]*S[:,3]+S[:,4]*S[:,4]+S[:,5]*S[:,5]))/2)
###skvtk = pv.VtkData(pv.UnstructuredGrid(points = NCS,tetra=TetT),'skull unique mesh',pv.CellData(pv.Scalars(VonM,'VonMises')))
###skvtk.tofile(fname+'.vtk')
###np.ma.dump(VonM,fname+'VonMises')

###fname = 'SkullPO_00_10'
###S = np.ma.load('Skull100Stress')[:,0:6] - np.ma.load('Skull000Stress10')[:,0:6]
###NCS = (np.ma.load('SkullPO_NC100')+np.ma.load('SkullPO_NC000_10'))/2
###TetT = np.ma.load('SkullPO_TetUse')
###VonM = np.sqrt((np.power(S[:,0]-S[:,1],2)+np.power(S[:,1]-S[:,2],2)+np.power(S[:,0]-S[:,2],2)+6*(S[:,3]*S[:,3]+S[:,4]*S[:,4]+S[:,5]*S[:,5]))/2)
###skvtk = pv.VtkData(pv.UnstructuredGrid(points = NCS,tetra=TetT),'skull unique mesh',pv.CellData(pv.Scalars(VonM,'VonMises')))
###skvtk.tofile(fname+'.vtk')
###np.ma.dump(VonM,fname+'VonMises')

###fname = 'SkullPO_00_20'
###S = np.ma.load('Skull100Stress')[:,0:6] - np.ma.load('Skull000Stress20')[:,0:6]
###NCS = (np.ma.load('SkullPO_NC100')+np.ma.load('SkullPO_NC000_20'))/2
###TetT = np.ma.load('SkullPO_TetUse')
###VonM = np.sqrt((np.power(S[:,0]-S[:,1],2)+np.power(S[:,1]-S[:,2],2)+np.power(S[:,0]-S[:,2],2)+6*(S[:,3]*S[:,3]+S[:,4]*S[:,4]+S[:,5]*S[:,5]))/2)
###skvtk = pv.VtkData(pv.UnstructuredGrid(points = NCS,tetra=TetT),'skull unique mesh',pv.CellData(pv.Scalars(VonM,'VonMises')))
###skvtk.tofile(fname+'.vtk')
###np.ma.dump(VonM,fname+'VonMises')

###fname = 'SkullPO_10_00'
###S = np.ma.load('Skull100Stress10')[:,0:6] - np.ma.load('Skull000Stress')[:,0:6]
###NCS = (np.ma.load('SkullPO_NC100_10')+np.ma.load('SkullPO_NC000'))/2
###TetT = np.ma.load('SkullPO_TetUse')
###VonM = np.sqrt((np.power(S[:,0]-S[:,1],2)+np.power(S[:,1]-S[:,2],2)+np.power(S[:,0]-S[:,2],2)+6*(S[:,3]*S[:,3]+S[:,4]*S[:,4]+S[:,5]*S[:,5]))/2)
###skvtk = pv.VtkData(pv.UnstructuredGrid(points = NCS,tetra=TetT),'skull unique mesh',pv.CellData(pv.Scalars(VonM,'VonMises')))
###skvtk.tofile(fname+'.vtk')
###np.ma.dump(VonM,fname+'VonMises')

###fname = 'SkullPO_10_10'
###S = np.ma.load('Skull100Stress10')[:,0:6] - np.ma.load('Skull000Stress10')[:,0:6]
###NCS = (np.ma.load('SkullPO_NC100_10')+np.ma.load('SkullPO_NC000_10'))/2
###TetT = np.ma.load('SkullPO_TetUse')
###VonM = np.sqrt((np.power(S[:,0]-S[:,1],2)+np.power(S[:,1]-S[:,2],2)+np.power(S[:,0]-S[:,2],2)+6*(S[:,3]*S[:,3]+S[:,4]*S[:,4]+S[:,5]*S[:,5]))/2)
###skvtk = pv.VtkData(pv.UnstructuredGrid(points = NCS,tetra=TetT),'skull unique mesh',pv.CellData(pv.Scalars(VonM,'VonMises')))
###skvtk.tofile(fname+'.vtk')
###np.ma.dump(VonM,fname+'VonMises')

###fname = 'SkullPO_10_20'
###S = np.ma.load('Skull100Stress10')[:,0:6] - np.ma.load('Skull000Stress20')[:,0:6]
###NCS = (np.ma.load('SkullPO_NC100_10')+np.ma.load('SkullPO_NC000_20'))/2
###TetT = np.ma.load('SkullPO_TetUse')
###VonM = np.sqrt((np.power(S[:,0]-S[:,1],2)+np.power(S[:,1]-S[:,2],2)+np.power(S[:,0]-S[:,2],2)+6*(S[:,3]*S[:,3]+S[:,4]*S[:,4]+S[:,5]*S[:,5]))/2)
###skvtk = pv.VtkData(pv.UnstructuredGrid(points = NCS,tetra=TetT),'skull unique mesh',pv.CellData(pv.Scalars(VonM,'VonMises')))
###skvtk.tofile(fname+'.vtk')
###np.ma.dump(VonM,fname+'VonMises')

###fname = 'SkullPO_20_00'
###S = np.ma.load('Skull100Stress20')[:,0:6] - np.ma.load('Skull000Stress')[:,0:6]
###NCS = (np.ma.load('SkullPO_NC100_20')+np.ma.load('SkullPO_NC000'))/2
###TetT = np.ma.load('SkullPO_TetUse')
###VonM = np.sqrt((np.power(S[:,0]-S[:,1],2)+np.power(S[:,1]-S[:,2],2)+np.power(S[:,0]-S[:,2],2)+6*(S[:,3]*S[:,3]+S[:,4]*S[:,4]+S[:,5]*S[:,5]))/2)
###skvtk = pv.VtkData(pv.UnstructuredGrid(points = NCS,tetra=TetT),'skull unique mesh',pv.CellData(pv.Scalars(VonM,'VonMises')))
###skvtk.tofile(fname+'.vtk')
###np.ma.dump(VonM,fname+'VonMises')

###fname = 'SkullPO_20_10'
###S = np.ma.load('Skull100Stress20')[:,0:6] - np.ma.load('Skull000Stress10')[:,0:6]
###NCS = (np.ma.load('SkullPO_NC100_20')+np.ma.load('SkullPO_NC000_10'))/2
###TetT = np.ma.load('SkullPO_TetUse')
###VonM = np.sqrt((np.power(S[:,0]-S[:,1],2)+np.power(S[:,1]-S[:,2],2)+np.power(S[:,0]-S[:,2],2)+6*(S[:,3]*S[:,3]+S[:,4]*S[:,4]+S[:,5]*S[:,5]))/2)
###skvtk = pv.VtkData(pv.UnstructuredGrid(points = NCS,tetra=TetT),'skull unique mesh',pv.CellData(pv.Scalars(VonM,'VonMises')))
###skvtk.tofile(fname+'.vtk')
###np.ma.dump(VonM,fname+'VonMises')

###fname = 'SkullPO_20_20'
###S = np.ma.load('Skull100Stress20')[:,0:6] - np.ma.load('Skull000Stress20')[:,0:6]
###NCS = (np.ma.load('SkullPO_NC100_20')+np.ma.load('SkullPO_NC000_20'))/2
###TetT = np.ma.load('SkullPO_TetUse')
###VonM = np.sqrt((np.power(S[:,0]-S[:,1],2)+np.power(S[:,1]-S[:,2],2)+np.power(S[:,0]-S[:,2],2)+6*(S[:,3]*S[:,3]+S[:,4]*S[:,4]+S[:,5]*S[:,5]))/2)
###skvtk = pv.VtkData(pv.UnstructuredGrid(points = NCS,tetra=TetT),'skull unique mesh',pv.CellData(pv.Scalars(VonM,'VonMises')))
###skvtk.tofile(fname+'.vtk')
###np.ma.dump(VonM,fname+'VonMises')