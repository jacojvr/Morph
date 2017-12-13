# -*- coding: utf-8 -*-

#execfile('SkullSymmMorph.py')
#print 'Do CG optimization iteration 1'
#import pysparse as ps
Tet = TetT
#NC = np.ma.load('NCopt_temp8')
#NCprev = NC
#Gmat = ps.spmatrix.ll_mat_from_mtx('Gmat9')
#GTG = ps.spmatrix.dot(Gmat,Gmat)
#gVec = ps.spmatrix.ll_mat(NC.shape[0]+outer.size,3)
#for j in outer:
  #rows = np.array([j,j,j])+NC.shape[0]
  #cols = np.array([0,1,2])
  #vals = NC[j,]
  #gVec.put(vals,rows,cols)
#GTg = ps.spmatrix.dot(Gmat,gVec)
#GTgNUMPY = np.zeros(NC.shape)
#GTgNUMPY[GTg.keys()] = GTg.values()
##mu = GTg.norm('inf')
##mu=mu*mu
##if mu<1:
#mu=1
#GTG.update_add_at(mu*np.ones((NC.shape[0],)),np.array(range(NC.shape[0])),np.array(range(NC.shape[0])))
#RHS = GTgNUMPY+mu*NCprev
#X,Y,Z = np.zeros((NC.shape[0],)),np.zeros((NC.shape[0],)),np.zeros((NC.shape[0],))
#print 'X'
#infx,itex,resx = ps.itsolvers.cgs(GTG,RHS[:,0],X,0.00000000001,100)
#print 'Y'
#infy,itey,resy = ps.itsolvers.cgs(GTG,RHS[:,1],Y,0.00000000001,100)
#print 'Z'
#infz,itez,resz = ps.itsolvers.cgs(GTG,RHS[:,2],Z,0.00000000001,100)
#NCcur = np.c_[X,Y,Z]
EQ,delt,Sn2,Sig = qu.elemQual_mu(np.array(range(Tet.shape[0])),NCcur,Tet)
print '			Max ',np.max(EQ),'; Min ',np.min(EQ),'; Average ',np.average(EQ)