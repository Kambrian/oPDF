#using nested_views rather than jointEL
import matplotlib
matplotlib.use('Agg')
from dynio import *
from myutils import *
import h5py,os,sys
from scipy.stats import chi2


estimator=10
proxy='LE'
nbins=[10,10]
npart=1000

lib.open()
#FullSample=Tracer('Mock')
FullSample=Tracer('AqA4N',DynRMAX=100)
Sample=FullSample.copy(0,npart)
FullSample.clean()

from scipy.optimize import *
x0=10**((np.random.rand(2)-0.5)*2)
print 'Initial: ', x0
print '-------------raw est---------------'
like= lambda x: Sample.jointE_Flike(x, estimator, nbinE=1)
x=fmin(like, x0, xtol=0.001, ftol=1e-4)
print x,Sample.joint_Flike(x, estimator, proxy, nbins)
print '------------jointL --------------------'
like= lambda x:Sample.jointLE_Flike(x, estimator, nbinL=nbins[0]*nbins[1], nbinE=1)
x=fmin(like, x0, xtol=0.001, ftol=1e-4)
print x,Sample.joint_Flike(x, estimator, proxy, nbins)
print '------------jointLE --------------------'
like= lambda x:Sample.jointLE_Flike(x, estimator, nbinL=nbins[0], nbinE=nbins[1])
x=fmin(like, x0, xtol=0.001, ftol=1e-4)
print x, Sample.joint_Flike(x, estimator, proxy, nbins)
print '--------------iterative fit----------------------'
like= lambda x:Sample.nested_views_Flike(x, estimator) #note: nested_views_like() does not work here.
  
x=x0
x0=[x[0]+1,x[1]]
while abs(x0[0]-x[0])>0.01:
  x0=x
  Sample.create_nested_views(x0, proxy, nbins)
  x=fmin(like, x0, xtol=0.001, ftol=1e-4)
  print x, like(x)
  
Sample.clean()
lib.close()