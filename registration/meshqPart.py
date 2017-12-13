# -*- coding: utf-8 -*-
#NCS = NC
#TetT = np.ma.load('SkTetAverTET')
######EQ,delt,Sn2,Sig = qu.elemQual_mu(np.array(range(TetT.shape[0])),NC,TetT)

#TetT = np.r_[Tet]
#whereInv = np.where(Sig<0)[0]
#Tempo = TetT[whereInv,0]
#TetT[whereInv,0]=TetT[whereInv,1]
#TetT[whereInv,1]=Tempo



#####EQ = EQ2
#####Sig = Sig2
fname = 'SkullUnique2_s20'





##TetDeg = np.where(EQ<0.15)[0]
##DegNd = TetT[TetDeg,]
##DegNd = DegNd.reshape((DegNd.size,))
##DegRep = np.array(find_repeats(DegNd)[0],int)
##for i in DegRep:
  ##DegNd = DegNd[DegNd<>i]
##DegNd = np.r_[DegNd,DegRep]
##DegNd.sort()
#PointConst = np.zeros((NCS.shape[0],))
#PointConst[outer,] = 1
##PointConst[DegNd,] = 0
#PointConst = np.array(PointConst,int)


#print 'Construct VTK object for optimization'
###ADD CONTRAINTS AS POINT SCALARS
#skvtk = pv.VtkData(pv.UnstructuredGrid(points = NCS,tetra=TetT),'skull unique mesh',pv.PointData(pv.Scalars(PointConst,'fixed')))
#skvtk.tofile(fname+'.vtk')


##os.system('~/Software/Trilinos/Mesquite/msqshape -s ./'+fname+'.vtk ./'+fname+'OUT.vtk')




fil = open(fname+'OUT.vtk','r+')
fil.writelines('# vtk DataFile Version 2.0\n')
fil.close()
skvtk = pv.VtkData(fname+'OUT.vtk')
NCS = np.array(skvtk.structure.points)