# -*- coding: utf-8 -*-
import ptetmeshdeform as tet
NCTnc = np.ma.load('BSkullFTnc')
TetT = np.ma.load('BSkullFTtet')
TriTet = np.ma.load('BSkullFTtri')
outer = np.ma.load('BSkullFTn_outer')
inner = np.ma.load('BSkullFTn_inner')
neighbList = np.ma.load('BSkullFTneighbTri')
neighbListTet = np.ma.load('BSkullFTneighbTet')
#N1 = outer.size
#N1T = inner.size
#LIM = 24
#NPP1 = N1/LIM
#NPP1T = N1T/LIM

#print "Set up 1-ring neighbor list for all points on the generic mesh surface"
#neighbList = [[0]]*N1
#results = pprocess.Map(limit=LIM)
#calc = results.manage(pprocess.MakeParallel(tet.Get1neigh))
#for j in range(0,LIM):
  #calc(outer[np.array(range(0,NPP1))+j*NPP1],NCTnc,TriTet)
#for j in range(0,LIM):
  #neighbList[j*NPP1:(1+j)*NPP1] = results[j]
#neighbList[LIM*NPP1:N1]=tet.Get1neigh(outer[np.array(range(LIM*NPP1,N1))],NCTnc,TriTet)
#np.ma.dump(neighbList,'BSkullFTneighbTri')
#print "Set up 1-ring neighbor list for all points internal to the generic mesh"
#neighbListTet = [[0]]*N1T
#results = pprocess.Map(limit=LIM)
#calc = results.manage(pprocess.MakeParallel(tet.Get1neigh))
#for j in range(0,LIM):
  #calc(inner[np.array(range(0,NPP1T))+j*NPP1T],NCTnc,TetT)
#for j in range(0,LIM):
  #neighbListTet[j*NPP1T:(1+j)*NPP1T] = results[j]
#neighbListTet[LIM*NPP1T:N1T]=tet.Get1neigh(inner[np.array(range(LIM*NPP1T,N1T))],NCTnc,TetT)
#np.ma.dump(neighbListTet,'BSkullFTneighbTet')