#using nested_views rather than jointEL
from dynio import *
from myutils import *
import h5py,os,sys
from scipy.stats import chi2


#estimator=8
#proxy='LE'
#nbins=[10,10]
npart=1000

lib.open()
FullSample=Tracer('Mock')
#FullSample=Tracer('AqA4N',DynRMAX=100)
Sample=FullSample.copy(5000,npart)
FullSample.clean()
#plt.ion()
x=Sample.fmin_jointLE(8, 100, 1)
#like=lambda m,c, beta,cc, g1, g2: -lib.wenting_like(lib.ParType(m,c,beta,cc,g1,g2), Sample._pointer)
#print like(1.5401/1.873, 28.0619/16.3349, 0.7173/0.715, 84.3133/69.014, 3.0/2.301, 8.2016/7.467)
#m=Minuit(like, m=1,c=1,beta=1,cc=1,g1=1,g2=1, print_level=3, pedantic=False, errordef=1, frontend=ConsoleFrontend())
#m.set_strategy(0)
#m.tol=10   #default convergence edm<1e-4*tol*errordef, but we do not need that high accuracy
#m.migrad()
#Sample.wenting_like_marginal([1,1])
#from scipy.optimize import *
#x0=10**((np.random.rand(2)-0.5)*2)
#x0=[2,1]
#print 'Initial: ', x0
#print '-------------raw est---------------'
#like= lambda x: -Sample.freeze_and_like(x,estimator)
#x=fmin(like, x0, xtol=0.001, ftol=1e-4)
#print x,-Sample.freeze_and_like(x,estimator)

  
#Sample.clean()
#lib.close()