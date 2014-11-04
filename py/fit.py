from dynio import *
from myutils import *
#plt.ion()

## fit halo
estimator=4 #4: rbin; 8: AD
nbin_r=50 #number of radial bins
Sample=Tracer('myhalo', DynDataFile='myhalo.hdf5') #load the data
Sample.radial_count(nbin_r) #prepare radial binning, only required if doing RBin fit
pars0=[1., 1.] #initial parameters, in units specified by the DynM0 and DynC0 in DataFiles.cfg (using DEFAULT if no section).
x1=Sample.gfmin_like(estimator, [1.,1.]) #fit
print 'x1', x1

## scan likelihood contour

##est error
sigm=P2Sig(chi2.sf(ts.min(0),1))
sigc=P2Sig(chi2.sf(ts.min(0),2))

## ML fit of density profile

## plot density profile


Sample.clean() #free the memory