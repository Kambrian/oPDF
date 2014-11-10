from math import *
import numpy as np
import ctypes,os
from myutils import fmin_gsl,density_of_points,get_extent,NamedEnum,NamedValues
#from scipy.optimize import fmin, fmin_powell
from scipy.stats import norm,chi2
from iminuit import Minuit
from iminuit.ConsoleFrontend import ConsoleFrontend
#import copy

if os.uname()[1]=='Medivh':
    rootdir='/work/Projects/DynDistr/'
else:
    rootdir='/gpfs/data/jvbq85/DynDistr/'

#=======================load the library===============================
lib=ctypes.CDLL("../liboPDF.so")
#=======================prototype the library==========================
#==globals.h
class global_tol(ctypes.Structure):
  _fields_=[('bin', ctypes.c_double),
			('bin_abs',ctypes.c_double),
			('rel', ctypes.c_double)]
class global_cosm(ctypes.Structure):
  _fields_=[('OmegaM0',ctypes.c_double),
			('OmegaL0',ctypes.c_double),
			('Redshift',ctypes.c_double)]
class global_const(ctypes.Structure):
  _fields_=[('G',ctypes.c_double),
			('H0',ctypes.c_double)]
class global_units(ctypes.Structure):
  _fields_=[('MassInMsunh',ctypes.c_double),
			('LengthInKpch',ctypes.c_double),
			('VelInKms',ctypes.c_double),
			('Const',global_const)]
lib.set_units.restype=None
lib.set_units.argtypes=[ctypes.c_double, ctypes.c_double, ctypes.c_double]
class globals_t(ctypes.Structure):
  _fields_=[('tol',global_tol),
			('cosmology',global_cosm),
			('units',global_units)]
  def set_defaults(self):
	'''set default global parameters, 
	including precision, cosmology and units'''
	lib.default_global_pars()
  
  def set_units(self, MassInMsunh=0.73e10, LengthInKpch=0.73, VelInKms=1.):
	'''set system of units,
	specify Mass in Msun/h, Length in kpc/h, Velocity in km/s'''
	lib.set_units(MassInMsunh, LengthInKpch, VelInKms)
VirTypes=NamedEnum('TH C200 B200')	
Globals=globals_t.in_dll(lib, 'Globals')

#===halo.h
MaxNPar=10
Param_t=ctypes.c_double*MaxNPar
class Halo_t(ctypes.Structure):
  '''all quantities are physical'''
  _fields_=[('pars', Param_t), #parameter values
			('scales', Param_t), #parameter scales. real parameters are set by pars*scales
			('z', ctypes.c_double),
			('M', ctypes.c_double),
			('c', ctypes.c_double),
			('Rv',ctypes.c_double),
			('Pots',ctypes.c_double),#-4*pi*G*rhos*rs^2, the potential at r=0
			('Rs',ctypes.c_double),
			('Rhos',ctypes.c_double),
			('Ms',ctypes.c_double),#4*pi*rs^3*rhos
			('RScale',ctypes.c_double),#for TMP profile, Rs/Rs0
			('PotScale',ctypes.c_double),#for TMP profile, Pots/Pots0
			#('TMPid',ctypes.c_int),#for TMP profile
			('virtype',ctypes.c_int),
			('type',ctypes.c_int)
			]
Halo_p=ctypes.POINTER(Halo_t) 
lib.halo_set_type.restype=None
lib.halo_set_type.argtypes=[ctypes.c_int, ctypes.c_int, ctypes.c_double, Param_t, Halo_p, ctypes.c_int]
lib.halo_set_param.restype=None
lib.halo_set_param.argtypes=[Param_t, Halo_p]
lib.halo_mass.restype=ctypes.c_double
lib.halo_mass.argtypes=[ctypes.c_double, Halo_p]
lib.halo_pot.restype=ctypes.c_double
lib.halo_pot.argtypes=[ctypes.c_double, Halo_p]

HaloTypes=NamedEnum('NFWMC NFWPotsRs NFWRhosRs TMPMC TMPPotScaleRScale CoreRhosRs CorePotsRs')
class Halo(Halo_t):
  '''general halo'''
  def __init__(self, halotype=HaloTypes.NFWMC, virtype=VirTypes.C200, redshift=0., scales=None, TMPid=-1):
	'''define a halo by specifiying the parametrization, virial definition and redshift of halo
	halotype: halo parametrization, one of the HaloTypes members
	virtype: virial definition, one of the VirTypes members
	redshift: redshift of halo
	scales: scales of halo parameters, array-like, of the same shape as parameters. default to all-ones if None. physical parameters will be params*scales
	TMPid: template id. only required when halotype is of template type'''
	Halo_t.__init__(self)
	#print "initing halo"
	self.set_type(halotype, virtype, redshift, scales, TMPid)

  def set_type(self, halotype=HaloTypes.NFWMC, virtype=VirTypes.C200, redshift=0., scales=None, TMPid=-1):
	'''set the parametrization, virial definition and redshift of halo
	halotype: halo parametrization, one of the HaloTypes members
	virtype: virial definition, one of the VirTypes members
	redshift: redshift of halo
	scales: scales of halo parameters, array-like, of the same shape as parameters. default to ones if not specified. physical parameters will be params*scales'''	
	if scales is None:
	  scales=np.ones(MaxNPar)
	lib.halo_set_type(halotype.value, virtype.value, redshift, Param_t(*scales), ctypes.byref(self), TMPid)
	
  def set_param(self, pars=[1.,1.]):
	'''set the parameters of the halo
	pars: parameters describing the halo'''
	lib.halo_set_param(Param_t(*pars), ctypes.byref(self))
  def mass(self, r):
	'''cumulative mass profile'''
	return lib.halo_mass(r, ctypes.byref(self))
  def pot(self, r):
	'''potential'''
	return lib.halo_pot(r, ctypes.byref(self))
  def get_current_TMPid():
	'''get the id of the template currently loaded in the system.
	this func can be used to check whether the loaded template 
	is the template of the current halo, just in case the template does not match'''
	return lib.get_current_TMPid(ctypes.byref(self))

#==tracer.h
class Particle_t(ctypes.Structure):
  _fields_=[('haloid', ctypes.c_int),
		('subid', ctypes.c_int),
		#('strmid', ctypes.c_int),
		('flag', ctypes.c_int),
		('w', ctypes.c_double),
	    ('r', ctypes.c_double),
	    ('K', ctypes.c_double),
	    ('L2', ctypes.c_double),#L^2
	    ('x', ctypes.c_double*3),
	    ('v', ctypes.c_double*3),
	    ('E', ctypes.c_double),
	    ('T', ctypes.c_double),
	    ('vr', ctypes.c_double),
	    ('theta', ctypes.c_double),
	    ('rlim', ctypes.c_double*2)
	    ]  
Particle_p=ctypes.POINTER(Particle_t)  

class Tracer_t(ctypes.Structure):
  pass
Tracer_p=ctypes.POINTER(Tracer_t)
Tracer_t._fields_=[('lnL', ctypes.c_double),
				   ('nP', ctypes.c_int),
				   ('mP', ctypes.c_double),
				   ('P', Particle_p),
				   ('nbin_r', ctypes.c_int),
				   ('FlagRLogBin', ctypes.c_int),
				   ('RadialCount', ctypes.POINTER(ctypes.c_double)),
				   ('rmin', ctypes.c_double),
				   ('rmax', ctypes.c_double),
				   ('proxybin', ctypes.c_double*2),
				   ('halo', Halo),
				   ('nView', ctypes.c_int),
				   ('ViewType', ctypes.c_char),
				   ('Views', Tracer_p)
				  ]
lib.load_tracer_particles.restype=None
lib.load_tracer_particles.argtypes=[ctypes.c_char_p, Tracer_p]
lib.cut_tracer_particles.restype=None
lib.cut_tracer_particles.argtypes=[Tracer_p, ctypes.c_double, ctypes.c_double]
lib.shuffle_tracer_particles.restype=None
lib.shuffle_tracer_particles.argtypes=[ctypes.c_ulong, Tracer_p]
lib.squeeze_tracer_particles.restype=None
lib.squeeze_tracer_particles.argtypes=[Tracer_p]
lib.free_tracer_particles.restype=None
lib.free_tracer_particles.argtypes=[Tracer_p]
lib.print_tracer_particle.restype=None
lib.print_tracer_particle.argtypes=[Tracer_p, ctypes.c_int]
lib.resample_tracer_particles.restype=None
lib.resample_tracer_particles.argtypes=[ctypes.c_ulong, Tracer_p, Tracer_p]
lib.copy_tracer_particles.restype=None
lib.copy_tracer_particles.argtypes=[ctypes.c_int, ctypes.c_int, Tracer_p, Tracer_p]
lib.count_tracer_radial.restype=None
lib.count_tracer_radial.argtypes=[Tracer_p, ctypes.c_int, ctypes.c_int]
lib.free_tracer_rcounts.restype=None
lib.free_tracer_rcounts.argtypes=[Tracer_p]
lib.sort_part_flag.restype=None
lib.sort_part_flag.argtypes=[Particle_p, ctypes.c_int]
lib.sort_part_E.restype=None
lib.sort_part_E.argtypes=[Particle_p, ctypes.c_int]
lib.sort_part_L.restype=None
lib.sort_part_L.argtypes=[Particle_p, ctypes.c_int]
lib.sort_part_R.restype=None
lib.sort_part_R.argtypes=[Particle_p, ctypes.c_int]
lib.create_tracer_views.restype=None
lib.create_tracer_views.argtypes=[Tracer_p, ctypes.c_int, ctypes.c_char]
lib.free_tracer_views.restype=None
lib.free_tracer_views.argtypes=[Tracer_p]
lib.create_nested_views.restype=None
lib.create_nested_views.argtypes=[ctypes.POINTER(ctypes.c_int), ctypes.c_char_p, Tracer_p];
lib.free_tracer.restype=None
lib.free_tracer.argtypes=[Tracer_p]
lib.tracer_set_energy.restype=None
lib.tracer_set_energy.argtypes=[Tracer_p]
lib.tracer_set_orbits.restype=None
lib.tracer_set_orbits.argtypes=[Tracer_p, ctypes.c_int]

#==nfw.h
lib.NFW_like.restype=ctypes.c_double
lib.NFW_like.argtypes=[Param_t, Tracer_p]

#==template.h
lib.get_current_TMPid.restype=ctypes.c_int
lib.get_current_TMPid.argtypes=[]

#==models.h
class NamedValuesEst(NamedValues):
  def __init__(self, value, name):
	NamedValues.__init__(self,value, name)
	self.need_phase=True
class NamedEstimators(object):
  def __init__(self, names):
	for number, name in enumerate(names.split()):
	  setattr(self, name, NamedValuesEst(number, name))
Estimators=NamedEstimators('RBinLike AD MeanPhase MeanPhaseRaw')
Estimators.RBinLike.need_phase=False
Estimators.RBinLike.nbin_r=20
Estimators.RBinLike.logscale=True

lib.alloc_integration_space.restype=None
lib.alloc_integration_space.argtypes=[]
lib.free_integration_space.restype=None
lib.free_integration_space.argtypes=[]

lib.like_eval.restype=ctypes.c_double
lib.like_eval.argtypes=[ctypes.c_int, Tracer_p]
lib.nested_views_like.restype=ctypes.c_double
lib.nested_views_like.argtypes=[ctypes.c_int, Tracer_p, ctypes.c_int, ctypes.c_int]
lib.DynFit.restype=ctypes.c_int
lib.DynFit.argtypes=[ctypes.POINTER(ctypes.c_double), ctypes.c_int, ctypes.c_double, ctypes.c_double, ctypes.c_int, ctypes.c_int, ctypes.c_int, Tracer_p]
lib.predict_radial_count.restype=None
lib.predict_radial_count.argtypes=[ctypes.POINTER(ctypes.c_double), ctypes.c_int, ctypes.c_int, Tracer_p]

def lib_open():
  '''initialize the library'''
  Globals.set_defaults()
  lib.alloc_integration_space()
  
def lib_close():  
  '''finalize the library'''
  lib.free_integration_space()
  
class Tracer(Tracer_t):
  ''' Tracer: a population of tracer particles.'''
  def __init__(self, datafile=None, rmin=None, rmax=None, shuffle=True):
	'''load a tracer from the datafile. 
	optionally, can apply radial cut given by rmin and rmax
	by default, the tracer particles will be shuffled after loading, for easy creation of subsamples by copying later.
	to keep the original ordering of particles, set shuffle=False'''
	Tracer_t.__init__(self)
	self._pointer=ctypes.byref(self)
	self.nP=0
	self.nView=0
	self.halo=Halo()
	if datafile!=None:
	  self.load(datafile)
	  self.radial_cut(rmin,rmax)
	  if shuffle:
		self.shuffle()

  def __enter__(self):
	return self
  
  def __exit__(self, type, value, traceback):
	self.clean()

  def __update_array(self):
	'''this should be called whenever the pointer self.P has changed
	to point the numpy array to the new memory'''
	Parr=(Particle_t*self.nP).from_address(ctypes.addressof(self.P.contents)) #struct array
	self.data=np.frombuffer(Parr, np.dtype(Parr))[0] #the numpy array

  def load(self, datafile):
	'''load particles from datafile'''
	if self.nP>0:
	  lib.free_tracer(self._pointer)
	lib.load_tracer_particles(datafile, self._pointer)
	print  self.nP, self.P
	self.__update_array()
	self.rmin=self.data['r'].min()
	self.rmax=self.data['r'].max()

  #def fill(x,v):
	#'''load particles from array. To be implemented.'''
	#FIXME.
	#pass
	
  def clean(self):
	''' never define it as __del__ and rely on the garbage collector.
	it's dangerous. gc may never call your __del__.
	call it yourself.'''
	#print "cleaning Tracer ", id(self)
	lib.free_tracer(self._pointer)
	#print self.nView

  def copy(self, offset=0, n=0):
	'''create a subsample by copying n particles starting from offset.
	if n==0, then copy all the particles starting from offset.
	return the subsample'''
	newsample=Tracer()
	lib.copy_tracer_particles(int(offset), int(n), newsample._pointer, self._pointer)
	newsample.__update_array()
	return newsample
  
  def select(self, flags):
	'''select particles according to flags array, to create a subsample'''
	sample=self.copy(0,0)
	sample.data['flag']=flags
	sample.squeeze() # __update_array() is called automatically
	return sample

  def squeeze(self):
	'''remove P.flag==0 (data['flag']==0) particles'''
	lib.squeeze_tracer_particles(self._pointer)
	self.__update_array()
	
  def shuffle(self, seed=1024):
	'''shuffle particles randomly.
	seed: optional, seeds the random number generator for the shuffle'''
	lib.shuffle_tracer_particles(ctypes.c_ulong(seed), self._pointer)

  def sort(self, proxy, offset=0,n=0):
	'''sort the particles according to proxy
	proxy can be 'E','L2','r' or 'flag'.
	offset, n: optional, sort n particles starting from offset.
	n=0 means sort all particles starting from offset.'''
	if n==0:
	  n=self.nP-offset
	sort_func={'flag': lib.sort_part_flag, 'E': lib.sort_part_E, 'L2': lib.sort_part_L, 'r': lib.sort_part_R}
	sort_func[proxy](ctypes.byref(Particle_t.from_buffer(self.P[offset])), n)

  def resample(self, seed=1024):
	''' create a bootstrap sample (sampling with replacement) of the same size from tracer
	return the new sample'''
	newsample=Tracer()
	lib.resample_tracer_particles(seed, newsample._pointer, self._pointer)
	newsample.__update_array()
	return newsample

  def radial_cut(self, rmin=None, rmax=None):
	'''cut the tracer with bounds [rmin, rmax]. 
	if only rmin or rmax is given, the other bound is not changed.'''
	if not ((rmin is None) and (rmax is None)):
	  if rmin is None:
		rmin=self.rmin
	  if rmax is None:
		rmax=self.rmax
	  lib.cut_tracer_particles(self._pointer, rmin, rmax)
	  self.__update_array()
	  
  def radial_count(self, nbin=10, logscale=True):
	'''bin the particles radially, to be used for radial likelihood calculation
	The histogram will be recorded in Tracer.RadialCount[].'''
	lib.count_tracer_radial(self._pointer, nbin, logscale)

  def predict_radial_count(self, nbin=100, logscale=True):
	'''predict radial counts according to oPDF.
	Tracer must be prepared with a halo (set_energy,set_orbits) 
	before predicting.
	return predicted counts.'''
	n=np.empty(nbin,dtype='f8')
	lib.predict_radial_count(n.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), nbin, logscale, self._pointer)
	#lib.predict_radial_count(n.ctypes.data, nbin) //also work?
	return n
  		
  def gen_bin(self,bintype,nbin=30,logscale=True, equalcount=False):
	'''return bin edges. divide into nbin bins, with nbin+1 edges.'''
	if bintype=='L':
	  bintype='L2'
	proxy=self.data[bintype]
	n=nbin+1
	if equalcount:
		x=np.array(np.percentile(proxy,list(np.linspace(0,100,n))))
	else:  
		if logscale:  
			x=np.logspace(np.log10(proxy[proxy>0].min()),np.log10(proxy.max()),n)
		else:
			x=np.linspace(proxy.min(),proxy.max(),n)
	return x
  
  def create_views(self, n=10, proxy='L'):
	'''sort the particles according to proxy, 
	and divide into n equal-size subsamples sequentially.
	these subsamples does not copy the particle data,
	but only points to the corresponding segments of data in the parent sample,
	so they are called views, and can be accessed through Tracer.Views[i] from
	the parent sample
	the energy need to have been set before calling if proxy is E'''
	lib.create_tracer_views(self._pointer, n, proxy)
	
  def destroy_views(self):
	'''erase any views from a tracer'''
	lib.free_tracer_views(self._pointer)  
  
  def print_data(self, i=0):
	'''print information of particle i'''
	print '---from python---'
	print self.nP, '%0x'%ctypes.addressof(self.P.contents)
	print self.P[i].x[0], self.P[i].r
	print self.data[i]['x'][0], self.data[i]['r']
	print '---  from C   ---'
	lib.print_tracer_particle(self._pointer, i)
	print '================='

  def NFW_like(self, pars=[1,1]):
	'''NFW log-likelihood. a NFW halo should have been attached before calling this
	pars are the parameters to be passed to the halo.
	return log(likelihood)'''
	return lib.NFW_like(Param_t(*pars), self._pointer)
	
  def set_energy(self):
	'''determine particle energy in the attached halo'''
	lib.tracer_set_energy(self._pointer)

  def set_orbits(self, set_phase=True):
	'''prepare particle orbits inside the attached halo
	set_phase: whether to calculate the phase-angle of each particle. 
	           phase angle is needed by AD and MeanPhase estimators,
	           but not by RBinLike'''
	lib.tracer_set_orbits(self._pointer, set_phase)

  def like_eval(self, estimator):
	'''evaluate likelihood or fig of merit with the given estimator in the attached halo.
	one has to call set_energy() and set_orbits() before this.
	'''
	return lib.like_eval(estimator.value, self._pointer)
  
  def likelihood(self, pars, estimator):
	'''calculate likelihood. automatically attach halo, prepare orbits and eval like'''
	self.halo.set_param(pars)
	self.set_energy()
	self.set_orbits(estimator.need_phase)
	return self.like_eval(estimator)

  def create_nested_views(self, viewtypes='EL', nbins=[10,10]):
	'''create nested views, i.e., create views according to first proxy,
	then create sub-views for each view according to the second proxy and so on.
	viewtypes can be one, two or more proxies, e.g, 'E','EL','LEr'.
	len(nbins) must match len(viewtypes).
	the energy need to be set before calling if creating E views'''
	try:
	  nbins=list(nbins)
	except:
	  nbins=[nbins]
	lib.create_nested_views((ctypes.c_int*(len(nbins)+1))(*nbins), ctypes.c_char_p(viewtypes), self._pointer)
  
  def nested_views_like(self, estimator=Estimators.AD):
	'''evaluate likelihood in at the deepest views, and return the sum of them.
	The likelihood for each view is also availabel in Views[i].lnL'''
	nbin_r=0
	logscale=False
	if estimator==Estimators.RBinLike:
	  nbin_r=estimator.nbin
	  logscale=estimator.logscale
	return lib.nested_views_like(estimator.value, self._pointer, nbin_r, logscale)
  
  def scan_like(self, estimator, x, y):
	'''scan a likelihood surface to be used for contour plots as contour(x,y,z)'''
	z=[self.likelihood([m,c], estimator) for m in x for c in y]
	return x,y,np.array(z).reshape([len(y),len(x)], order="F")
    
  def TSprof(self, pars, estimator=Estimators.MeanPhase, proxy='L', nbin=100):
	'''calculate the likelihood inside equal-count bins of proxy.
	return the loglike or f.o.m. for the estimator in each bin, and the bin edges.
	proxy and nbin can also be of len>1; in that case, use 
	self.Views[i].Views[j].lnL and self.Views[i].Views[j].proxybin
	to get the likelihood and bins in each node'''
	self.halo.set_param(pars)
	self.set_energy()
	self.set_orbits()
	if viewtype=='L2':
	  viewtype='L'
	self.create_nested_views(viewtype, nbin)
	self.nested_views_like(estimator)
	ts=[self.Views[i].lnL for i in xrange(self.nView)] #these are like_eval, not chi2
	ts.append(np.nan)
	ts=np.array(ts)
	x=self.gen_bin(viewtypes, nbins, equalcount=True)
	return ts,x
  
  def TSprofCum(self, pars, estimator=Estimators.AD, proxy='r', bins=100):
	'''cumulative TS profile
	reuturn bin edges, ts, counts'''
	self.halo.set_param(pars)
	self.set_energy()
	self.sort(proxy)
	if proxy=='E':
	  self.data[...]=self.data[-1::-1] #reverse the order
	n,bins=np.histogram(self.data[proxy],bins)
	if proxy=='E':
	  n=n[-1::-1].cumsum()
	else:
	  n=n.cumsum()
	ts=[]
	if proxy!='r':
	  self.set_orbits(True)
	for i,x in enumerate(bins[1:]):
	  if n[i]==0:
		y=np.nan
	  else:
		with self.copy(0,n[i]) as s:
		  if proxy!='r':
			y=s.like_eval(pars, estimator)
		  else:
			s.rmax=x
			s.set_orbits()
			y=s.like_eval(pars,estimator)
	  ts.append(y)
	if proxy=='E':
	  ts=ts[-1::-1]
	return np.array(bins),np.array(ts),np.array(n)

  def NFW_fit(self, x0=[1,1], minuittol=1):
	''' to fit an NFW density PDF with maximum likelihood.
	This is only intended for fitting the Dark Matter density profile to get the NFW parameters.
	The tracer particle mass should have been properly assigned or adjusted, 
	so that mP*number_density=physical_density.
	If you have sampled n particles from the full sample of n0 particles, 
	remember to adjust the mP of the sample to be mP0*n0/n, so that total mass is conserved.
	results will be printed on screen.
	also return minuit result and the minuit minimizer'''
	like=lambda m,c: -self.NFW_like([m,c])
	#profilelikelihood ratio error-def: chi-square1
	m=Minuit(like, m=x0[0], c=x0[1], print_level=0, pedantic=False, error_m=0.1, error_c=0.1, errordef=0.5,frontend=ConsoleFrontend())
	m.tol=minuittol   #default convergence edm<1e-4*tol*errordef, but we do not need that high accuracy
	result=m.migrad()
	m.print_fmin()
	m.print_matrix()
	return result,m
  
  def dyn_fit(self, estimator, x0=[1,1], xtol=1e-3, ftol_abs=0.01, maxiter=500, verbose=0):
	''' dynamical fit with the given estimator
	input: 	estimator: estimator to use. select one from Estimators.
			x0: initial parameter values
			xtol: tolerance in x to consider convergence
			ftol_abs: tolerance in function values to consider convergence. 
			          convergence is reached when both dx<xtol and df<ftol_abs between subsequent steps in the search.
			maxiter: maximum number of iterations
			verbose: whether to print during each step.
	return: [x, fval, status_success], 
			x: the best fit parameter  
			fval: log-likelihood or fig of merit, depending on estimator
			status_success: whether the search converged successfully, 1 if yes, 0 if no.
			'''
	status_success=lib.DynFit(Param_t(*x0), len(x0), xtol, ftol_abs, maxiter, verbose, estimator.value, self._pointer)	
	x=np.array(self.halo.pars[:len(x0)])
	return x,self.lnL,status_success

  def phase_density(self, proxy='E', bins=100, method='hist', logscale=False, weight=False):
	'''estimate density in proxy-theta space'''
	if logscale:
	  f=self.data[proxy]>0
	  data=np.array((np.log10(self.data[proxy][f]), self.data['theta'][f]))
	  w=self.data['w'][f]
	else:
	  data=np.array((self.data[proxy], self.data['theta']))
	  w=self.data['w']
	if weight:
	  X,Y,Z=density_of_points(data, bins=bins, method=method, weights=w)
	else:
	  X,Y,Z=density_of_points(data, bins=bins, method=method)
	extent=get_extent(X,Y)
	return X,Y,Z,extent,data

  #def gfmin_like(self, estimator, x0=[1,1], xtol=1e-3, ftolabs=0.01):
	#like=lambda x: lib.like_to_chi2(self.freeze_and_like(x, estimator), estimator) #the real likelihood prob
	#return fmin_gsl(like, x0, xtol=xtol, ftolabs=ftolabs, maxiter=500, full_output=True)
  
  #def gfmin_dist(self, estimator, x0=[1,1], xtol=1e-3, ftolabs=0.01):
	#like=lambda x: -self.freeze_and_like(x, estimator) #distance
	#return fmin_gsl(like, x0, xtol=xtol, ftolabs=ftolabs, maxiter=500, full_output=True)

lib_open() #allocate integration space

if __name__=="__main__":
  datafile=rootdir+'/data/mockhalo_wenting.hdf5'
  FullSample=Tracer(datafile=datafile)
  Sample=FullSample.copy(0,1000)
  FullSample.print_data(10)
  Sample.print_data(10)
  #Sample.data[1]['r']=2
  Sample.radial_cut(1,1000)
  Sample.print_data(1)
  Sample.sort('L2')
  Sample.print_data(1)
  Sample.halo.set_type(HaloTypes.NFWMC, scales=[183.5017,16.1560])
  Sample.radial_count(10)
  result=Sample.dyn_fit(Estimators.RBinLike,verbose=1)
  print Sample.likelihood([1,1], Estimators.MeanPhase)
  #FullSample.clean()
  #Sample.clean()
  
  #good practice: use with!
  #with Tracer(datafile) as FullSample:
	#with FullSample.copy() as Sample:
	  #FullSample.print_data()
	  #Sample.print_data()
	  #Sample.data[1]['r']=2
	  #Sample.radial_cut(rmin=10)
	  #Sample.print_data()
  #lib.close()