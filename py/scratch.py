#using nested_views rather than jointEL
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
Sample=FullSample.copy(5000,npart)
FullSample.clean()
#plt.ion()

#def plot_pot(pars, linestyle='-'):
  #print pars
  #r=np.logspace(0,2,20)
  #Halo=NFWHalo()
  #Halo.define_halo(pars)
  #y=-np.array(map(Halo.pot, r))
  #plt.loglog(r/Halo.halo.Rs,y, linestyle)
  
#plot_pot([1,1],'r')
#plot_pot([0.6013295 ,  2.22171445], 'g')
#plot_pot([0.82423562, 1.28484559], 'c')
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