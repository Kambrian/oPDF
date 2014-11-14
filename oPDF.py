""" Python interface to the C code for oPDF modelling.

It wraps the C functions into python classes with ctypes."""
from math import *
import numpy as np
import ctypes,os
from myutils import Chi2Sig,AD2Sig,density_of_points,get_extent,NamedEnum,NamedValues
#from myutils import fmin_gsl
#from scipy.optimize import fmin, fmin_powell
import matplotlib.pyplot as plt
from iminuit import Minuit
from iminuit.ConsoleFrontend import ConsoleFrontend
#import copy

rootdir=os.path.dirname(__file__)
if rootdir=='':
  rootdir='.'
#=======================load the library===============================
lib=ctypes.CDLL(rootdir+"/liboPDF.so")
#=======================prototype the library==========================
#==globals.h
class global_tol(ctypes.Structure):
  _fields_=[('bin', ctypes.c_double),
			('bin_abs',ctypes.c_double),
			('rel', ctypes.c_double)]
class global_cosm(ctypes.Structure):
  _fields_=[('OmegaM0',ctypes.c_double),
			('OmegaL0',ctypes.c_double)]
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
  ''' global variables of the module. It controls numerical precision, internal units, and cosmology. '''
  _fields_=[('tol',global_tol),
			('cosmology',global_cosm),
			('units',global_units)]
  def set_defaults(self):
	'''set default global parameters, 
	including precision, cosmology and units'''
	lib.default_global_pars()
  
  def set_units(self, MassInMsunh=0.73e10, LengthInKpch=0.73, VelInKms=1.):
	'''set system of units.
	specify Mass in Msun/h, Length in kpc/h, Velocity in km/s'''
	lib.set_units(MassInMsunh, LengthInKpch, VelInKms)

  def get_units(self):
	'''query the units'''
	print 'Mass  :', self.units.MassInMsunh, 'Msun/h'
	print 'Length:', self.units.LengthInKpch, 'kpc/h'
	print 'Vel   :', self.units.VelInKms, 'km/s'
	return (self.units.MassInMsunh, self.units.LengthInKpch, self.units.VelInKms)
VirTypes=NamedEnum('TH C200 B200')	
Globals=globals_t.in_dll(lib, 'Globals')

#===halo.h
MaxNPar=10
Param_t=ctypes.c_double*MaxNPar
class Halo_t(ctypes.Structure):
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
  '''a general halo describing the potential. It has the following properties
  
  :ivar pars: raw parameter values. do not change them manually, use :func:`set_param` to set them.
  :ivar scales: parameter scales. use :func:`set_type` to set them.
  :ivar virtype: virial definition
  :ivar type: parametrization type. One of :data:`HaloTypes`.
  
  Depending on the type of the halo, some of the following properties may be calculated during :func:`set_param`:
  
  :ivar M: mass 
  :ivar c: concentration
  :ivar Rv: virial radius
  :ivar Pots: Pots=4*pi*G*Rhos*Rs^2.
  :ivar Rhos: scale density for NFW
  :ivar Rs: scale radius 
  :ivar RScale: Rs/Rs0 for TMP profile
  :ivar PotScale: Pots/Pots0 for TMP profile
  '''
  def __init__(self, halotype=HaloTypes.NFWMC, virtype=VirTypes.C200, redshift=0., scales=None, TMPid=-1):
	'''define a halo by specifiying the parametrization, virial definition and redshift of halo
	
	halotype: halo parametrization, one of the :data:`HaloTypes` members
	
	virtype: virial definition, one of the VirTypes members
	
	redshift: redshift of halo
	
	scales: scales of halo parameters, array-like, of the same shape as parameters. default to all-ones if None. physical parameters will be params*scales
	
	TMPid: template id. only required when halotype is of template type'''
	Halo_t.__init__(self)
	#print "initing halo"
	self.set_type(halotype, virtype, redshift, scales, TMPid)

  def set_type(self, halotype=HaloTypes.NFWMC, virtype=VirTypes.C200, redshift=0., scales=None, TMPid=-1):
	'''set the parametrization, virial definition and redshift of halo
	
	halotype: halo parametrization, one of the :data:`HaloTypes` members
	
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
	try:
	  m=lib.halo_mass(r, ctypes.byref(self))
	except:
	  m=np.array([lib.halo_mass(x, ctypes.byref(self)) for x in r])
	return m
  def pot(self, r):
	'''potential'''
	try:
	  m=lib.halo_pot(r, ctypes.byref(self))
	except:
	  m=np.array([lib.halo_pot(x, ctypes.byref(self)) for x in r])
	return m
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
	    ('K', ctypes.c_double), #kinetic energy
	    ('L2', ctypes.c_double),#L^2
	    ('L', ctypes.c_double), #angular momentum
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
				  ]#: Tracer fields
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
Estimators.RBinLike.nbin=20
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

def _lib_open():
  '''initialize the library'''
  Globals.set_defaults()
  lib.alloc_integration_space()
  
def _lib_close():  
  '''finalize the library'''
  lib.free_integration_space()
  
class Tracer(Tracer_t):
  ''' Tracer: a population of tracer particles.
    
  :ivar halo: the halo (potential, :py:class:oPDF.Halo) for the tracer.
  :ivar lnL: likelihood or distance for the sample from the previous likelihood calculation, depending on estimator.
  :ivar nP: number of particles.
  :ivar mP: average particle mass.
  :ivar data: particle data, numpy record array format. It includes the following fields:
	 ('haloid', 'subid', 'flag', 'w', 'r', 'K', 'L2', 'L', 'x', 'v', 'E', 'T', 'vr', 'theta', 'rlim')
	 
  .. note::
	 - The `w` field is the particle mass in units of the average particle mass. These are all ones if no particle mass is given in the datafile.
	 - the `haloid` and `subid` fields are only filled if you have `SubID` and `HaloID` datasets in the datafile when loading. 
	 - The `E`,`theta` and `rlim` fields are the energy, phase-angle, and radial limits (peri and apo-center distances) of the orbits.These depend on the potential, and are only filled when you have done some calculation in a halo, or have filled them explicitly with :py:func:`set_phase`.
  
  .. note::
     The following members are provided for information, but do not manually assign to them. use :py:func:`radial_count` and :py:func:`radial_cut` to set them.
  
  :ivar nbin_r: number of radial bins.
  :ivar FlagRLogBin: whether radial binning is in logspace.
  :ivar RadialCount: counts in radial bins. 
  :ivar rmin: lower radial cut.
  :ivar rmax: upper radial cut.
  '''
  def __init__(self, datafile=None, rmin=None, rmax=None, shuffle=True):
	'''it loads a tracer from the datafile. 
	
	optionally, can apply radial cut given by rmin and rmax
	
	.. note:: 
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
	'''release the C-allocated memory for the tracer'''
	#never define it as __del__ and rely on the garbage collector.
	#it's dangerous. gc may never call your __del__.
	#call it yourself.
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
	
	proxy can be 'E','L','r' or 'flag'.
	
	offset, n: optional, sort n particles starting from offset.
	  n=0 means sort all particles starting from offset.'''
	if n==0:
	  n=self.nP-offset
	sort_func={'flag': lib.sort_part_flag, 'E': lib.sort_part_E, 'L': lib.sort_part_L, 'r': lib.sort_part_R}
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
	'''bin the particles radially, to be used for radial likelihood calculation.
	The histogram will be recorded in Tracer.RadialCount[].'''
	lib.count_tracer_radial(self._pointer, nbin, logscale)

  def predict_radial_count(self, nbin=100, logscale=True):
	'''predict radial counts according to oPDF.
	
	:func:`set_phase` must have been called prior to calling this.
	
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
	the parent sample.
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
	'''NFW log-likelihood. the halo should have been set to one of the NFW types before calling this.
	
	pars are the parameters to be passed to the halo.
	
	return log(likelihood)'''
	return lib.NFW_like(Param_t(*pars), self._pointer)
	
  def set_energy(self):
	'''determine particle energy in the attached halo'''
	lib.tracer_set_energy(self._pointer)

  def set_orbits(self, set_phase=True):
	'''prepare particle orbits inside the attached halo
	
	set_phase: whether to calculate the phase-angle of each particle. 
	  phase angle is needed by AD and MeanPhase estimators, but not by RBinLike
	  see estimator.need_phase for each estimator.
	'''
	lib.tracer_set_orbits(self._pointer, set_phase)

  def set_phase(self, pars, need_theta=True):
	'''prepare the phases for phase-related calculations such as like_eval or phase_density'''
	self.halo.set_param(pars)
	self.set_energy()
	self.set_orbits(need_theta)

  def like_eval(self, estimator):
	'''evaluate likelihood or fig of merit with the given estimator in the attached halo.
	one has to call :func:`set_phase` before this.
	'''
	return lib.like_eval(estimator.value, self._pointer)
  	
  def likelihood(self, pars, estimator, auto_rbin=True):
	'''calculate likelihood. automatically prepare binning, orbits and eval like.'''
	if auto_rbin and estimator==Estimators.RBinLike:
	  self.radial_count(estimator.nbin, estimator.logscale)
	self.set_phase(pars, estimator.need_phase)
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
	'''scan a likelihood surface.
	
	x,y are the vectors specifying the binning along x and y dimensions
	
	return the likelihood value z on grids, to be used for contour plots as 
	  >>> contour(x,y,z)
	'''
	if estimator==Estimators.RBinLike:
	  self.radial_count(estimator.nbin, estimator.logscale)
	z=[self.likelihood([m,c], estimator, False) for m in x for c in y]
	return np.array(z).reshape([len(y),len(x)], order="F")
  
  def scan_confidence(self, estimator, x0, ngrids=[10,10], dx=[0.5, 0.5], logscale=False, maxlike=None):
	'''scan significance levels around parameter value x0.
	
	it scans ngrids linear bins from x0-dx to x0+dx if logscale=False,
	or ngrids log bins from log10(x0)-dx to log10(x0)+dx if logscale=True.
	
	If maxlike is given, it is interpreted as the global maximum log-likelihood, and is used to determine significance for RBinLike estimator;
	otherwise the maximum likelihood is automatically scanned for RBinLike.
	'''
	if logscale:
	  xl=np.log10(x0)-dx
	  xu=np.log10(x0)+dx
	  x=np.logspace(xl[0], xu[0], ngrids[0])
	  y=np.logspace(xl[1], xu[1], ngrids[1])
	else:
	  xl=x0-dx
	  xu=x0+dx
	  x=np.linspace(xl[0], xu[0], ngrids[0])
	  y=np.linspace(xl[1], xu[1], ngrids[1])
	like=self.scan_like(estimator, x, y)
	if estimator==Estimators.RBinLike:
	  if maxlike is None:
		xfit,maxlike,status=self.dyn_fit(estimator, x0)
	  sig=Chi2Sig(2*(maxlike-like), 2)
	elif estimator==Estimators.AD:
	  sig=AD2Sig(like)
	elif estimator==Estimators.MeanPhase:
	  sig=np.sqrt(like)
	elif estimator==Estimators.MeanPhaseRaw:
	  sig=np.abs(like)
	return x,y,sig
    
  def TSprof(self, pars, proxy='L', nbin=100, estimator=Estimators.MeanPhaseRaw):
	'''calculate the likelihood inside equal-count bins of proxy.
	
	return the loglike or f.o.m. for the estimator in each bin, and the bin edges.
	
	proxy and nbin can also be of len>1; in that case, use 
	self.Views[i].Views[j].lnL and self.Views[i].Views[j].proxybin
	to get the likelihood and bins in each node'''
	self.halo.set_param(pars)
	self.set_energy()
	if proxy!='r': #r-view orbits will be updated internally.
	  self.set_orbits()
	self.create_nested_views(proxy, nbin)
	self.nested_views_like(estimator)
	ts=[self.Views[i].lnL for i in xrange(self.nView)]
	ts=np.array(ts)
	x=self.gen_bin(proxy, nbin, equalcount=True)
	return ts,x
  
  def plot_TSprof(self, pars, proxy='L', nbin=100, estimator=Estimators.MeanPhaseRaw, xtype='percent-phys', linestyle='r-'):
	'''plot the TS profile
	
	xtype: can be one of 'percent', 'physical', and 'percent-phys'. 
	  
	  when xtype='percent', plot the x-axis with percents.
	  
	  if xtype='phys', plot x-axis with physical values. 
	  
	  if xtype='percent-phys', plot xaxis in percent scale but label with physical values.
	'''
	ts,x=self.TSprof(pars, proxy, nbin, estimator)
	percents=(np.arange(nbin)+0.5)/nbin*100
	xmid=(x[:-1]+x[1:])/2
	if xtype=='percent':
	  h,=plt.plot(percents, ts, linestyle)
	  plt.xlabel('Percentile in '+proxy)
	if xtype=='percent-phys':
	  h,=plt.plot(percents, ts, linestyle)
	  plt.xticks(percents[1::3],['%.1f'%np.log10(a) for a in x][1::3])
	  plt.xlabel(r'$\log('+proxy+')$')
	if xtype=='phys':
	  h,=plt.plot(xmid, ts, linestyle)
	  plt.xlabel(r'$%s$'%proxy)
	plt.ylabel(estimator)
	if estimator==Estimators.MeanPhaseRaw:
	  plt.plot(plt.xlim(), [0, 0], 'k--')
	return h 
	
  def TSprofCum(self, pars, proxy='r', bins=100, estimator=Estimators.AD):
	'''cumulative TS profile.
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
	results will be printed on screen.
	also return minuit result and the minuit minimizer.
	
	..note:
	This is only intended for fitting the Dark Matter density profile to get the NFW parameters.
	The tracer particle mass should have been properly assigned or adjusted, 
	so that mP*number_density=physical_density.
	If you have sampled n particles from the full sample of n0 particles, 
	remember to adjust the mP of the sample to be mP0*n0/n, so that total mass is conserved.
	'''
	like=lambda m,c: -self.NFW_like([m,c])
	#profilelikelihood ratio error-def: chi-square1
	m=Minuit(like, m=x0[0], c=x0[1], print_level=0, pedantic=False, error_m=0.1, error_c=0.1, errordef=0.5,frontend=ConsoleFrontend())
	m.tol=minuittol   #default convergence edm<1e-4*tol*errordef, but we do not need that high accuracy
	result=m.migrad()
	m.print_fmin()
	m.print_matrix()
	return result,m
  
  def dyn_fit(self, estimator=Estimators.RBinLike, x0=[1,1], xtol=1e-3, ftol_abs=0.01, maxiter=500, verbose=0):
	''' 
	dynamical fit with the given estimator
	
	Parameters
		estimator(Estimator): estimator to use. select one from :data:`Estimators`.
		
		x0(array-like): initial parameter values
		
		xtol: tolerance in x to consider convergence
		
		ftol_abs: tolerance in function values to consider convergence. 
		
		convergence is reached when both dx<xtol and df<ftol_abs between subsequent steps in the search.
		
		maxiter: maximum number of iterations
		
		verbose: whether to print during each step.
	
	Returns
	  [x, fval, status_success] 
		x(array): the best fit parameter  
		
		fval(float): log-likelihood or fig of merit, depending on estimator
		
		status_success(bool): whether the search converged successfully, 1 if yes, 0 if no.
	'''
	if estimator==Estimators.RBinLike:
	  self.radial_count(estimator.nbin, estimator.logscale)
	status_success=lib.DynFit(Param_t(*x0), len(x0), xtol, ftol_abs, maxiter, verbose, estimator.value, self._pointer)	
	x=np.array(self.halo.pars[:len(x0)])
	return x,self.lnL,status_success

  def phase_density(self, proxy='E', bins=100, method='hist', logscale=False, weight=False, return_data=False):
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
	if return_data:
	  return X,Y,Z,extent,data
	else:
	  return X,Y,Z,extent
	
  def phase_image(self, pars, proxy, bins=30, logscale=True):
	'''plot an image of the particle distribution in proxy-theta space'''
	self.set_phase(pars) #determine the phase-angles with the real halo parameters x0
	X,Y,Z,extent=self.phase_density(proxy, bins=bins, logscale=logscale)
	plt.imshow(Z,extent=extent)
	plt.axis('tight')
	if logscale:
	  plt.xlabel(r'$\log('+proxy+')$')
	else:
	  plt.xlabel(r'$'+proxy+'$')
	plt.ylabel(r'$\theta$')

  #def gfmin_like(self, estimator, x0=[1,1], xtol=1e-3, ftolabs=0.01):
	#like=lambda x: lib.like_to_chi2(self.freeze_and_like(x, estimator), estimator) #the real likelihood prob
	#return fmin_gsl(like, x0, xtol=xtol, ftolabs=ftolabs, maxiter=500, full_output=True)
  
  #def gfmin_dist(self, estimator, x0=[1,1], xtol=1e-3, ftolabs=0.01):
	#like=lambda x: -self.freeze_and_like(x, estimator) #distance
	#return fmin_gsl(like, x0, xtol=xtol, ftolabs=ftolabs, maxiter=500, full_output=True)

_lib_open() #allocate integration space

if __name__=="__main__":
  datafile=rootdir+'/data/mockhalo.hdf5'
  FullSample=Tracer(datafile=datafile)
  Sample=FullSample.copy(0,1000)
  FullSample.print_data(10)
  Sample.print_data(10)
  #Sample.data[1]['r']=2
  Sample.radial_cut(1,1000)
  Sample.print_data(1)
  Sample.sort('L')
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
#_lib_close()