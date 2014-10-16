import h5py,os,sys
import numpy as np

halo=sys.argv[1]
starfile="/gpfs/data/wenting/starpop_output/APCBower/Aq_"+halo+"2_WMAP1_Cold/tags_han.hdf5"
subfile="/gpfs/data/wenting/starpop_output/APCBower/Aq_"+halo+"2_WMAP1_Cold/subhalo_membership_han.hdf5"
oldfile="/gpfs/data/jvbq85/DynDistr/data/old/"+halo+"2star.hdf5"
newfile="/gpfs/data/jvbq85/DynDistr/data/"+halo+"2star.hdf5"
f1=h5py.File(starfile,"r")
f2=h5py.File(subfile,"r")
f3=h5py.File(oldfile,"r")
f4=h5py.File(newfile,"w")

ids=f3['/Pid'][...]
DMID=f1['/PartType0/DMID'][...][ids]
_,ind,ind_recon=np.unique(DMID, return_index=True, return_inverse=True)
uids=ids[ind] #unique set of ids

#copy x,v,pid,subid
f4.create_dataset('/x', data=f3['/x'][...][ind])
f4['/x'].attrs['x0']=f3['/x'].attrs['x0']
f4.create_dataset('/v', data=f3['/v'][...][ind])
f4['/v'].attrs['v0']=f3['/v'].attrs['v0']
f4.create_dataset('/Pid', data=uids)
f4.create_dataset('/SubID', data=f1['/PartType0/SubhaloID'][...][uids], dtype='int32')

#stack particle mass
mnew=np.bincount(ind_recon, weights=f1['/PartType0/Mass'][...][ids].flat).reshape(-1,1) 
f4.create_dataset('/PartMass', data=mnew, dtype='float32')
#tag branch ids
treeind=f2['/PartType0/TreeIndex'][...]#all the branches
subid=f1['/PartType0/SubhaloID'][...]#all subhaloes
TreeFlag=np.in1d(treeind[uids].flat, treeind[subid>0].flat).reshape(-1,1) #tag all the alive branches
f4.create_dataset('/HaloID', data=TreeFlag, dtype='int32')#tree alive flag. >0 means alive; =0 means dead
f4['/HaloID'].attrs['Comment']='FlagBranchAlive; whether the merger tree branch is still alive. not really haloID'

print 'All done'

f1.close()
f2.close()
f3.close()
f4.close()