from dynio import *
from myutils import *
#plt.ion()

## fit halo
estimator=4 #4: rbin; 8: AD
nbin_r=50 #number of radial bins
with Tracer('Mock') as FullSample:#load the data
  Sample=FullSample.copy(0, 100) #pick a random subsample of 1000 particles.

##define the fittype and fit units
#Note: only assign to C_global_var.value, never assign to the var directly!!!
Halo=NFWHalo()
#set the meaning of pars[], 
#can be HaloFitTypes.(MC,RhosRs,PotsRs,CoreRhosRs,CorePotsRs). default is MC (mass, concentration).
Halo.FitType.value=HaloFitTypes.CoreRhosRs 
#set fit units. by default, the units of pars[] are determined from DynM0 and DynC0 in configuration file when loading Tracer(). 
#now we manually overide the units here. 
Halo.Rhos0.value=1e-4 #set scale of Rhos, in units of 10^10Msun/kpc^3
Halo.Rs0.value=10 #set scale of Rs, in units of kpc

Sample.radial_count(nbin_r) #prepare radial binning, only required if doing RBin fit
pars0=[1., 1.] #initial parameters, in units specified by the DynM0 and DynC0 in DataFiles.cfg (using DEFAULT if no section).
x1=Sample.gfmin_like(estimator, [1.,1.]) #fit
print 'x1', x1

# scan likelihood contour...

#Sample.clean() #free the memory