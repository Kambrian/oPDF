'''pick additional information to add to the set of star particles within 500kpc'''
import h5py,os,sys
import numpy as np

halo=sys.argv[1]
starfile="/gpfs/data/wenting/starpop_output/APCBower/Aq_"+halo+"2_WMAP1_Cold/tags_han.hdf5"
#subfile="/gpfs/data/wenting/starpop_output/APCBower/Aq_"+halo+"2_WMAP1_Cold/subhalo_membership_han.hdf5"
outfile="/gpfs/data/jvbq85/DynDistr/data/"+halo+"2star.hdf5"
f1=h5py.File(starfile,"r")
#f2=h5py.File(subfile,"r")
f3=h5py.File(outfile,"a")

ids=f3['/Pid'][...]
f3.create_dataset('/PartMass', data=f1['/PartType0/Mass'][...][ids])
f3.create_dataset('/SubID', data=f1['/PartType0/SubhaloID'][...][ids], dtype='int32')

DMID=f1['/PartType0/DMID'][...][ids]
_,ind=np.unique(DMID, return_index=True)
flag=np.zeros(ids.shape,dtype='int32')
flag[ind]=1
f3.create_dataset('/HaloID', data=flag)#unique set of DM particles
f3['/HaloID'].attrs['Note']='flags to select a unique set of DM particles; not really haloID'

f1.close()
#f2.close()
f3.close()

#treeind=f2['/PartType0/TreeIndex'][...][ids]
#TreeFlag=np.zeros([treeind.max()+1,1],dtype='int32')
#subid=f1['/PartType0/SubhaloID'][...][ids]
#for i in xrange(ids.shape[0]):
  #if subid[i]>0:
	#TreeFlag[treeind[i]]=1
#f3.create_dataset('/HaloID', data=TreeFlag[treeind])#tree branch selection