#using nested_views rather than jointEL
from dynio import *
from myutils import *
import h5py,os,sys
from scipy.stats import chi2
plt.ion()

#estimator=8
#proxy='LE'
#nbins=[10,10]
npart=1000

lib.open()
FullSample=Tracer('AqA2starN')
Sample=FullSample.copy(0,npart)
Sample.FlagUseWeight=0
print Sample.freeze_and_like([1,1], 8)
Sample.FlagUseWeight=1
print Sample.freeze_and_like([1,1], 8)

FullSample.clean()
Sample.clean()
lib.close()