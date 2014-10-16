from dynio import *
from myutils import *
#plt.ion()

estimator=4 #4: rbin; 8: AD
nbin_r=50
Sample=Tracer('myhalo', DynDataFile='myhalo.hdf5') #load the data
Sample.radial_count(nbin_r) #prepare radial binning
pars0=[1., 1.] #initial parameters, in units specified by the DynM0 and DynC0 in DataFiles.cfg (using DEFAULT if no section).
x1=Sample.gfmin_like(estimator, [1.,1.]) #fit
print 'x1', x1

Sample.clean() #free the memory