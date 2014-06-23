#using nested_views rather than jointEL
import matplotlib
matplotlib.use('Agg')
from dynio import *
from myutils import *
import h5py,os,sys
from scipy.stats import chi2


estimator=8
proxy='LE'
nbins=[10,10]
npart=1000

lib.open()
FullSample=Tracer('Mock')
#FullSample=Tracer('AqA4N',DynRMAX=100)
Sample=FullSample.copy(0,npart)
FullSample.clean()

from scipy.optimize import *
#x0=10**((np.random.rand(2)-0.5)*2)
x0=[2,1]
print 'Initial: ', x0
print '-------------raw est---------------'
like= lambda x: -Sample.freeze_and_like(x,estimator)
x=fmin(like, x0, xtol=0.001, ftol=1e-4)
print x,-Sample.freeze_and_like(x,estimator)

  
#Sample.clean()
#lib.close()