import h5py,os,sys
import numpy as np

halo=sys.argv[1]
subfile="/gpfs/data/wenting/starpop_output/APCBower/Aq_"+halo+"2_WMAP1_Cold/subhalo_membership_han.hdf5"
newfile="/gpfs/data/jvbq85/DynDistr/data/"+halo+"2star.hdf5"
f2=h5py.File(subfile,"r")
f4=h5py.File(newfile,"a")

uids=f4['/Pid'][...]
treeind=f2['/PartType0/TreeIndex'][...]#all the branches
_,_,ind_recon=np.unique(treeind[uids], return_index=True, return_inverse=True)
f4.create_dataset('/BranchID', data=ind_recon, dtype='int32')
f4['/BranchID'].attrs['Comment']='branch index. relabelled from 1-n using np.unique()'

print 'All done'

f2.close()
f4.close()