# -*- coding: utf-8 -*-
print
print "SELECTIVE MESH MORPHING ALGORITHM USING ELASTIC SURFACE REGISTRATION"
print "	-G.J.J.v.Rensburg - 22/04/2010-"
t_start = time.clock()
ConnectB = np.array(ConnectB,int)
ConnectT = np.array(ConnectT,int)
#LandmB = np.array(LandmB[:,0],int)		# do -1 later to be consistent with python indexing, first need to do other "temporary landmarks"& check that they dont fall on actual landmark positions!

# Settings for elastic surface registration:
m=10 # nearest neighbour parameter
alph=0.5 # normilization factor
gamm=2 # smoothing parameter1
sigm0=10 # smoothing parameter2
f=1.0715 # smoothing parameter3
Tol=0.0001 # stopping criteria

# determine N1,N2,T1 and T2:
N1=NCB.shape[0]
N2=NCT.shape[0]
T1=ConnectB.shape[0]
T2=ConnectT.shape[0]
NL=LandmB.shape[0]

################################     INITIALIZE & NODES OF CONCERN:    #################################
########################################################################################################

print
print "INITIALIZE SURFACE DEFORMATION"
k = 1
CONV = []
print " 	enquire nodes where required displacement is checked"
###remove Landmarks from FDNB and SFPNT:
#for i in range(0,NL):
  #if find_repeats(np.r_[UseN_B,LandmB[i,]])[0].size>0:
    #r=np.where(UseN_B==LandmB[i,])[0]
    #UseN_B = np.r_[UseN_B[0:r,],UseN_B[r+1:UseN_B.size,]]
SamplingB=UseN_B.size
SamplingT=UseN_T.size
## Full list of nodes used in Surface registration:
LMB = np.r_[UseN_B]#,LandmB]	# Last NL entries are reserved for Landmarks that HAVE TO FIT points on the target mesh
LMT = np.r_[UseN_T]


# For parallel programming devide Nr of computations by number of parallel processes (LIM)
NPP1 = N1/LIM
NPP2 = N2/LIM
SBPP = SamplingB/LIM
STPP = SamplingT/LIM
FMorph = 0

print "	Compute known displacement for Base_Landmarks "
knownC = NCB[LandmB,]
knownD = LandmB_NC-knownC
print	
print "COARSE SURFACE REGISTRATION"
print "	using landmark displacements to deform using RBF"
W_km1 = RBFmorph(NCB,NCB[LandmB,],LandmB_NC-NCB[LandmB,])

################################    MAIN MESH DEFORMATION ALGORITHM:   #################################
########################################################################################################
print
print "ELASTIC SURFACE REGISTRATION"
print "determine vertex normals of target surface"
#Compute target-mesh triangle centroids:
print "determining centroids of target surface triangles"
S_2_centr = np.c_[np.sum(np.c_[NCT[ConnectT[:,0],0],NCT[ConnectT[:,1],0],NCT[ConnectT[:,2],0]],1)/3,
  np.sum(np.c_[NCT[ConnectT[:,0],1],NCT[ConnectT[:,1],1],NCT[ConnectT[:,2],1]],1)/3,np.sum(np.c_[NCT[ConnectT[:,0],2],NCT[ConnectT[:,1],2],NCT[ConnectT[:,2],2]],1)/3]
print "determine triangle and vertex normals of target surface"
TNORMT = np.cross(NCT[ConnectT[:,1],:]-NCT[ConnectT[:,0],:],NCT[ConnectT[:,2],:]-NCT[ConnectT[:,0],:])
TNORMT = (TNORMT.T/(np.ones((3,1))*np.sqrt(np.sum(np.array([TNORMT*TNORMT]),2)))).T
VNORMT = vrtxnormal(NCT,ConnectT,S_2_centr,TNORMT)

print "determining kd-trees of target surface centroids and nodal coordinates"
KDT_TC = KDTree(S_2_centr,m)
KDT_TN = KDTree(NCT,m)


D1 = np.zeros((SamplingB,3))
D2 = np.zeros((SamplingT,3))
DS = np.zeros((N1,3))
print
print "MESH DEFORMATION ITERATION",k
print "	determining known displacement of landmarks"
knownD = LandmB_NC-W_km1[LandmB,]
print "	determining centroids of deforming mesh"
W_km1_centr = np.c_[np.sum(np.c_[W_km1[ConnectB[:,0],0],W_km1[ConnectB[:,1],0],W_km1[ConnectB[:,2],0]],1)/3,
  np.sum(np.c_[W_km1[ConnectB[:,0],1],W_km1[ConnectB[:,1],1],W_km1[ConnectB[:,2],1]],1)/3,np.sum(np.c_[W_km1[ConnectB[:,0],2],W_km1[ConnectB[:,1],2],W_km1[ConnectB[:,2],2]],1)/3]
print "	determine triangle and vertex normals of deforming surface"
TNORMB = np.cross(W_km1[ConnectB[:,1],:]-W_km1[ConnectB[:,0],:],W_km1[ConnectB[:,2],:]-W_km1[ConnectB[:,0],:])
TNORMB = (TNORMB.T/(np.ones((3,1))*np.sqrt(np.sum(np.array([TNORMB*TNORMB]),2)))).T
VNORMB = vrtxnormal(W_km1,ConnectB,W_km1_centr,TNORMB)

print "	determining kd-tree of current deforming surface centroids and nodal coordinates"
KDT_KC = KDTree(W_km1_centr,m)
KDT_KN = KDTree(W_km1,m)
