""" Python interface to the C code for oPDF modelling.

It wraps the C functions into python classes with ctypes.

TODO: integrate PhaseTickerFit with minuit or curvefit"""
from math import *
import numpy as np
import ctypes,os
from myutils import Chi2Sig,AD2Sig,density_of_points,get_extent,NamedValues,NamedEnum
#from myutils import fmin_gsl
#from scipy.optimize import fmin, fmin_powell
import matplotlib.pyplot as plt
from iminuit import Minuit
from iminuit.ConsoleFrontend import ConsoleFrontend
from scipy.optimize import newton,brentq,fmin,curve_fit
#import copy

oPDFdir=os.path.dirname(__file__)
if oPDFdir=='':
  oPDFdir='.'
#=======================load the library===============================
lib=ctypes.CDLL(oPDFdir+"/liboPDF.so")
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
  ''' global variables of the module. It controls numerical precision, internal units, and cosmology. 
  The numerical precision for orbit integration is controlled by Globals.tol.rel, which defaults to 1e-3 
  (Globals.tol.rel=1e-2 should already be sufficient for likelihood inference and phase calculations).'''
  _fields_=[('tol',global_tol),
			('cosmology',global_cosm),
			('units',global_units)]
  def set_defaults(self):
	'''set default global parameters, 
	including precision, cosmology and units'''
	lib.default_global_pars()
  
  def set_units(self, MassInMsunh=1.e10, LengthInKpch=1., VelInKms=1.):
	'''set system of units.
	specify Mass in Msun/h, Length in kpc/h, Velocity in km/s.
	
	:example:
	
	  If you want to use (1e10Msun, kpc, km/s) as units, and you adopt :math:`h=0.73`  in your model, then you can set the units like below
	
	  >>> h=0.73
	  >>> Globals.set_units(1e10*h,h,1)
	
	That is, to set them to (1e10h Msun/h, h kpc/h, km/s).
	
	.. note::
	   - The user should only use Globals.set_units() to change the units, which automatically updates several interal constants related to units. Never try to change the internal unit variables (e.g., Globals.units.MassInMsunh) manually.
	   -  To avoid inconsistency with units of previously loaded tracers, you must do it immediately after importing the :module:`oPDF` module if you need to call :func:`set_units`.'''
	   
	lib.set_units(MassInMsunh, LengthInKpch, VelInKms)

  def get_units(self):
	'''query the units'''
	print 'Mass  :', self.units.MassInMsunh, 'Msun/h'
	print 'Length:', self.units.LengthInKpch, 'kpc/h'
	print 'Vel   :', self.units.VelInKms, 'km/s'
	return (self.units.MassInMsunh, self.units.LengthInKpch, self.units.VelInKms)
class _NamedEnum(NamedEnum):
#just an alias of NamedEnum
#it is here just because sphinx does not auto-doc instances of imported classes,
#so we make it a native class by aliasing it.
  pass

VirTypes=_NamedEnum('TH C200 B200')
'''Collection of virial definitions. It contains 

    - ``VirTH``: the spherical collapse prediction (i.e, Bryan & Norman 98 fitting).
    - ``VirB200``: the 200 times mean density deifinition.
    - ``VirC200``: the 200 times critical density definition.
'''
Globals=globals_t.in_dll(lib, 'Globals')
'''    Collection of global variables of the module, of class :class:`globals_t`. It controls numerical precision, internal units, and cosmology. 
'''

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
			('K', ctypes.c_double), #M=KR, for isothermal profile, where K=Vc^2/G.
			('IsForbidden', ctypes.c_int), #whether the halo parameters are forbidden(e.g, negative mass)
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
lib.isNFW.restype=ctypes.c_int
lib.isNFW.argtypes=[ctypes.c_int]

HaloTypes=_NamedEnum('NFWMC NFWPotsRs NFWRhosRs TMPMC TMPPotScaleRScale CoreRhosRs CorePotsRs PointM IsothermalK')
''' Collection of halo types. It contains

    -  ``NFWMC``: NFW halo parametrized by :math:`(M,c)`
    -  ``NFWRhosRs``: NFW, :math:`(\\rho_s,r_s)`
    -  ``NFWPotsRs``: NFW, (:math:`|\\psi_s|, r_s`\ ), with :math:`|\\psi_s|=4\\pi G\\rho_s r_s^2`\ .
    -  ``CorePotsRs``: Cored Generalized NFW Potential (inner density slope=0), parametrized by (:math:`|\\psi_s|,r_s`\ )
    -  ``CoreRhosRs``: Cored GNFW, :math:`(\\rho_s,r_s)`
    -  ``TMPMC``: Template profile, :math:`(M,c)` parametrization
    -  ``TMPPotScaleRScale``: Template, :math:`\\psi_s/\\psi_{s0}, r_s/r_{s0}`
    -  ``PointM``: Point Mass at r=0
    -  ``IsothermalK``: Isothermal halo, :math:`M(r)=Kr`
'''

class Halo(Halo_t):
  '''a general halo describing the potential. It has the following properties
  
  :ivar pars: raw parameter values. do not change them manually, use :func:`set_param` to set them.
  :ivar scales: parameter scales. use :func:`set_type` to set them.
  :ivar virtype: virial definition. One of :const:`VirTypes'.
  :ivar type: parametrization type. One of :const:`HaloTypes`.
  
  Depending on the type of the halo, some of the following properties may be calculated during :func:`set_param`:
  
  :ivar M: mass 
  :ivar c: concentration
  :ivar Rv: virial radius
  :ivar Pots: :math:`\\psi_s=4\\pi G\\rho_s r_s^2`.
  :ivar Rhos: scale density for NFW
  :ivar Rs: scale radius 
  :ivar RScale: :math:`r_s/r_{s0}` for TMP profile
  :ivar PotScale: :math:`\psi_s/\psi_{s0}` for TMP profile
  
  '''
  def __init__(self, halotype=HaloTypes.NFWMC, virtype=VirTypes.C200, redshift=0., scales=None, TMPid=-1):
	'''
	:initializer:
	
	  define a halo by specifiying the parametrization, virial definition and redshift of halo
	  
	  :param halotype: halo parametrization, one of the :const:`HaloTypes` members
	  
	  :param virtype: virial definition, one of the :const:`VirTypes` members
	  
	  :param redshift: redshift of halo
	  
	  :param scales: scales of halo parameters, array-like, of the same shape as parameters. default to all-ones if None. physical parameters will be params*scales
	  
	  :param TMPid: template id. only required when halotype is of template type
	
	'''
	Halo_t.__init__(self)
	#print "initing halo"
	self.set_type(halotype, virtype, redshift, scales, TMPid)

  def set_type(self, halotype=HaloTypes.NFWMC, virtype=VirTypes.C200, redshift=0., scales=None, TMPid=-1):
	'''set the parametrization, virial definition and redshift of halo
	
	halotype: halo parametrization, one of the :data:`HaloTypes` members
	
	virtype: virial definition, one of the :const:`VirTypes` members
	
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
	'''cumulative mass profile
	:param r: array-like or float, the radius'''
	try:
	  m=lib.halo_mass(r, ctypes.byref(self))
	except:
	  m=np.array([lib.halo_mass(x, ctypes.byref(self)) for x in r])
	return m
  def pot(self, r):
	'''potential
	:param r: array-like or float, the radius'''
	try:
	  m=lib.halo_pot(r, ctypes.byref(self))
	except:
	  m=np.array([lib.halo_pot(x, ctypes.byref(self)) for x in r])
	return m
  def get_current_TMPid(self):
	'''get the id of the template currently loaded in the system.
	
	this func can be used to check whether the loaded template 
	is the template of the current halo, just in case the template does not match'''
	return lib.get_current_TMPid(ctypes.byref(self))

  def isNFW(self):
	''' return True if halo is NFW, False if not '''
	return bool(lib.isNFW(self.type))
  
  def vr_inv(self, r, E, L):
	''' inverse of radial velocity, 1/v_r, for a particle at r with binding energy E and angular momentum L'''
	vr2=2*(-E-0.5*(L/r)**2-self.pot(r))
	ivr=1./np.sqrt(vr2)
	try:
	  ivr[vr2<=0]=0.
	except:
	  if vr2<=0:
		ivr=0.
	return ivr
  
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
lib.load_tracer_particles.argtypes=[ctypes.c_char_p, Tracer_p, ctypes.c_int]
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
'''   Collection of dynamical estimators. It contains

    - ``RBinLike``: binned radial likelihood. 
    
	Use ``RBinLike.nbin`` (``integer``) and ``RBinLike.logscale`` (``True`` or ``False``) to control the number and scale of bins. 
	Since the purpose of the binning is purely to suppress shot noise, a larger number of bins is generally better, as long as it is not too noisy. On the other hand, when the likelihood contours appear too irregular, one should try reducing the number of radial bins to ensure the irregularities are not caused by shot noise. In our analysis, we have adopted 30 bins for an ideal sample of 1000 particles, and 50 bins for :math:`10^6` particles in a realistic halo, although a bin number as low as 5 could still work.
    
    - ``AD``: Anderson-Darling distance.
    - ``MeanPhaseRaw``: Normalized mean phase deviation :math:`\\bar{\\Theta}=(\\bar{\\theta}-0.5)/\\sigma_{\\theta}`\ , to be compared to a standard normal variable.
    - ``MeanPhase``: :math:`\\bar{\\Theta}^2`\ , to be compared to a chi-square variable.'''
Estimators.RBinLike.need_phase=False
Estimators.RBinLike.nbin=20 #: The number of radial bins.
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
    
  :ivar halo: the halo (potential, type :class:`Halo`) for the tracer.
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
  def __init__(self, datafile=None, rmin=None, rmax=None, shuffle=True, AddHubbleFlow=False):
	'''
	:Initializer: 
	
	loading a tracer from a datafile. 
	
	optionally, can apply radial cut given by rmin and rmax
	
	if AddHubbleFlow=True, then also add hubble flow to the loaded velocity (only support redshift 0 data now).
	
	.. note:: 
	   The datafile should contain physical positions and velocities of the tracer particles, relative to the center of the halo.
	   By default, the tracer particles will be shuffled after loading, for easy creation of subsamples by copying later.
	   To keep the original ordering of particles, set shuffle=False'''
	Tracer_t.__init__(self)
	self._pointer=ctypes.byref(self)
	self.nP=0
	self.nView=0
	self.halo=Halo()
	if datafile!=None:
	  self.load(datafile, AddHubbleFlow)
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

  def load(self, datafile, AddHubbleFlow=False):
	'''load particles from datafile'''
	if self.nP>0:
	  lib.free_tracer(self._pointer)
	lib.load_tracer_particles(datafile, self._pointer, AddHubbleFlow)
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
	
	Only particles and their radial limits are copied. The halo, RadialCounts and Views are not copied into the new sample.
	
	return the subsample'''
	newsample=Tracer()
	lib.copy_tracer_particles(int(offset), int(n), newsample._pointer, self._pointer)
	newsample.__update_array()
	return newsample
  
  def select(self, flags):
	'''select particles according to flags array, into a new sample. 
	
	.. note::
	   - Same as :func:`copy`, only particles and their radial limits are copied. The halo, RadialCounts and Views are not copied into the new sample.
	   - When doing dynamical tests, one should avoid distorting the radial distribution with any radial selection. One can still apply radial cuts, but should only do this with the :func:`radial_cut` function. So never use :func:`select` on data['r'].'''
	   
	sample=self.copy(0,0)
	sample.data['flag']=flags
	sample.squeeze() # __update_array() is called automatically
	return sample

  def squeeze(self):
	'''remove P.flag==0 (data['flag']==0) particles'''
	lib.squeeze_tracer_particles(self._pointer)
	self.__update_array()
	
  def shuffle(self, seed=100):
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

  def resample(self, seed=100):
	''' create a bootstrap sample (sampling with replacement) of the same size from tracer
	return the new sample'''
	newsample=Tracer()
	lib.resample_tracer_particles(seed, newsample._pointer, self._pointer)
	newsample.__update_array()
	return newsample

  def radial_cut(self, rmin=None, rmax=None):
	'''cut the tracer with bounds [rmin, rmax]. 
	if only rmin or rmax is given, the other bound is not changed.
	
	.. note::
	   This function not only selects particles within (rmin,rmax), but also sets the radial boundary for the dynamical model, so that only dyanmical consistency inside the selected radial range is checked. So always use this function if you want to change radial cuts. This function is automatically called when initializing a :class:`Tracer` with rmin/rmax.'''
	   
	if not ((rmin is None) and (rmax is None)):
	  if rmin is None:
		rmin=self.rmin
	  if rmax is None:
		rmax=self.rmax
	  lib.cut_tracer_particles(self._pointer, rmin, rmax)
	  self.__update_array()
	  
  def radial_count(self, nbin=10, logscale=True):
	'''bin the particles radially, to be used for radial likelihood calculation.
	The histogram will be recorded in Tracer.RadialCount[]. 
	
	.. note::
	   This function is automatically called by the relevant likelihood functions such as :func:`likelihood`,
	   :func:`dyn_fit`, :func:`scan_confidence` when the estimator is :const:`Estimators` ``.RBinLike``. In these cases,
	   `nbin` and `logscale` will be determined according to :attr:`Estimators.RBinLike.nbin` and :attr:`Estimators.RBinLike.logscale`.
	   So usually you do not need to call this function explicitly.'''
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
	'''evaluate likelihood or TS with the given estimator. 
	
	.. note::
	This is a low-leve function. One has to call :func:`set_phase` before calling this.
	Use :func:`likelihood` which combines :func:`like_eval` and :func:`set_phase` automatically, unless you do want to call them separately. 
	
	:param estimator: estimator to use for the likelihood or TS calculation.
	
	:returns: the likelihood (if estimator= :const:`Estimators` ``.RBin``) or TS value (:math:`\\bar{\\Theta}`\  for :const:`Estimators` ``.MeanPhaseRaw``, :math:`\\bar{\\Theta}^2`\  for :const:`Estimators` ``.MeanPhase``, or AD distance for :const:`Estimators` ``.AD``).
	'''
	return lib.like_eval(estimator.value, self._pointer)
  	
  def likelihood(self, pars=[1,1], estimator=Estimators.MeanPhaseRaw, auto_rbin=True):
	'''calculate likelihood or test-statistic.
	it automatically prepares orbits and then calculates the likelihood or TS
	(i.e., it combines :func:`set_phase` and :func:`like_eval`).
	
	:param pars: parameter values
	:param estimator: estimator to use for the likelihood or TS calculation.
	:param auto_rbin: whether to prepare radial bins automatically. Only relevant if you are using Estimators.RBin. default to True. set to false if you have prepared the bins manually (by calling :func:`radial_count`).
	
	:returns: the likelihood (if estimator= :const:`Estimators` ``.RBin``) or TS value (:math:`\\bar{\\Theta}`\  for :const:`Estimators` ``.MeanPhaseRaw``, :math:`\\bar{\\Theta}^2`\  for :const:`Estimators` ``.MeanPhase``, or AD distance for :const:`Estimators` ``.AD``).
	'''
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
	
	:returns: [x,y,sig,like]
	
		x,y: the scanned grid points, vectors.
		
		sig: the significance on the grids, of shape [len(x),len(y)]
		
		like: the likelihood or figure of merits on the grids. same shape as sig.
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
	return x,y,sig,like
    
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
	'''plot the TS profile in equal-count bins of proxy.
	
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
	  h,=plt.plot(np.log10(xmid), ts, linestyle)
	  plt.xlabel(r'$\log(%s)$'%proxy)
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

  def solve_meanphase(self, x0, y0=0., verbose=0):
	'''find the halo parameter x at which MeanPhaseRaw(x)=y0.
	x0 is the initial value of x to start the search.
	verbose: 0 or 1, whether to print additional convergence information.
	
	:return: (x,y,flag).
		x: solution
		y: MeanPhaseRaw(x)-y0
		success flag:
		  0: fail to converge
		  1: successfully found solution; 
		  2: failed to find solution for MeanPhaseRaw(x)-y0=0, but found minimum for (MeanPhaseRaw(x)-y0)**2.
	.. note::
		 The tracer halo type must be set before invoking this function.
	'''
	likeraw=lambda m: self.likelihood([m], Estimators.MeanPhaseRaw, False)-y0
	likeln=lambda x: likeraw(np.exp(x))
	like=lambda m: likeraw(m)**2
	#initialize the interval bracketing the root
	d=2.
	a=[x0*d,x0]
	f=[likeraw(x) for x in a]
	if abs(f[1])<abs(f[0]):
	  f=f[::-1]
	  a=a[::-1]
	  d=1./d
	niter=0
	while f[0]*f[1]>0 and niter<10:
	  a[1]=a[0]
	  a[0]*=d #push a[0] toward root
	  f[1]=f[0]
	  f[0]=likeraw(a[0])
	  niter+=1
	if niter>=10:
	  if verbose>0:
		print "Failed to bracket root after %d iterations; now search for func minimum instead."%niter
	  result=fmin(like, x0, full_output=True, disp=verbose)
	  success=2 #retreat to fmin
	  if result[4]>0:
		success=0 #failed to converge
	  return result[0][0], result[1], success
	
	#find the root
	#print a,f, np.log(a), [likeln(np.log(a[0])), likeln(np.log(a[1]))]
	x=brentq(likeln, np.log(a[0]),np.log(a[1]), rtol=1e-4)
	x0=np.sqrt(a[0]*a[1])
	#x=newton(likeln, np.log(x0),tol=1e-4)
	y=likeln(x)
	x=np.exp(x)
	#x=brentq(likeraw, a[0],a[1], rtol=1e-2)
	#x=newton(likeraw, x0,tol=1e-4)
	#y=likeraw(x)
	success=abs(y)<1e-2
	return x,y,success
  
  def tick_phase_mass(self, m0=100., verbose=0):
	''' estimate mass :math:`M(<r_c)` using the "PhaseTicker" method.
	:param m0: the initial guess of the mass, optional.
	:param verbose: whether to print diagnostic information, default no.
	
	:return: r,m,ml,mu,flag,flagl, flagu
		  r: the characteristic radius of the tracer in between rlim
		  m: the mass of the halo contained inside r
		  ml: 1-sigma lower bound on m
		  mu: 1-sigma upper bound on m
		  flag: whether m and r have converged (0:no; 1:yes; 2: no solution to the phase equation, but closest point found).
		  flagl: whether ml has converged
		  flagu: whether mu has converged
	'''
	
	rmed=np.median(self.data['r'])
	self.halo.set_type(HaloTypes.PointM)
	
	result=self.solve_meanphase(m0, verbose=verbose)
	flag1=result[-1]
	if verbose>0:
	  print result
	m=result[0]
	ml,_,flagl=self.solve_meanphase(m, -1., verbose=verbose)
	mu,_,flagu=self.solve_meanphase(m, 1., verbose=verbose)
	
	self.halo.set_type(HaloTypes.IsothermalK)
	result=self.solve_meanphase(m0/rmed, verbose=verbose)
	flag2=result[-1]
	if verbose>0:
	  print result
	k=result[0]
	#kl=self.solve_meanphase(k, -1.)[0]
	#ku=self.solve_meanphase(k, 1.)[0]
	
	r=m/k
	#rl,ru=rlim
	#rl=max(ml/ku,rlim[0])
	#ru=min(mu/kl,rlim[1])
	
	flag=flag1&flag2
	if flag>0:
	  flag=max(flag1, flag2) #in case any of them equal to 2
	  
	return r,m,ml,mu,flag,flagl,flagu 
	
  def phase_mass_select(self, xlim, proxy='r', m0=100., verbose=0):
	''' estimate mass :math:`M(<r_c)` for tracer with property "proxy" selected in between xlim, using the "PhaseTicker" method.
	:param xlim: the range of property to select the tracer.
	:param proxy: the property to select the tracer, 'r' or 'L'
	:param m0: the initial guess of the mass, optional.
	:param verbose: whether to print diagnostic information, default no.
	
	:return: r,m,ml,mu,flag,flagl, flagu
		  r: the characteristic radius of the tracer in between rlim
		  m: the mass of the halo contained inside r
		  ml: 1-sigma lower bound on m
		  mu: 1-sigma upper bound on m
		  flag: whether m and r have converged (0:no; 1:yes; 2: no solution to the phase equation, but closest point found).
		  flagl: whether ml has converged
		  flagu: whether mu has converged
	'''
	
	if proxy in 'rR':
	  S=self.copy()
	  S.radial_cut(xlim[0], xlim[1])
	else:
	  S=self.select((self.data[proxy]>xlim[0])&(self.data[proxy]<xlim[1]))
	out=S.tick_phase_mass(m0, verbose)
	S.clean()
	return out
  
  def phase_ticker_fit(self, par0=[1,1], nbin=2, proxy='r', equalcount=True):
	'''fit halo potential with phase ticker. The halo of the tracer need to be initialized to the desired type before fitting.
	
	:param par0: initial parameter values. len(par0) also gives the number of free parameters.
	:param nbin: number of bins. if nbin<len(par0), then the number of bins is set to len(par0) to avoid overfitting.
	:param proxy: the tracer property used to bin the sample into subsamples. 'r' or 'L'.
	:param equalcount: whether to use equalcount bins (each subsample has equal number of particles) or logarithmic bins in proxy.
	
	:return:
	  par: best-fit parameter
	  Cov: covariance matrix of the parameters
	  data: phase-ticker data, array of shape [nbin, 7]. each column is the fit to one bin, [r,m,ml,mu,flag,flagl, flagu], with (r,m) giving the radius and mass, (ml,mu) giving the lower and upper bound on mass, (flag, flagl, flagu) specifying whether the fit converged for mass and its lower and upper bounds (0:no; 1:yes; 2: no solution to the phase equation, but closest point found).
	  
	.. note::
	  if the code complains about curve_fit() keyword error, you need to upgrade your scipy to version 0.15.1 or newer.
	'''	
	
	nbin=max(nbin, len(par0))
	x=self.gen_bin(proxy, nbin, equalcount=equalcount)
	data=np.array([self.phase_mass_select(x[[i,i+1]], proxy) for i in xrange(len(x)-1)])
	flag=(data[:,4]==1)
	sig1=data[:,3]-data[:,1]
	sig2=data[:,1]-data[:,2]
	sig=(sig1+sig2)/2 
	def halomass(r, *pars):
	  self.halo.set_param(pars)
	  return self.halo.mass(r)
	par,Cov=curve_fit(halomass, data[flag,0], data[flag,1], sigma=sig[flag], p0=par0, absolute_sigma=1) #need scipy version >0.15.1
	return par,Cov,data
  
  def NFW_fit(self, x0=[1,1], minuittol=1):
	''' to fit an NFW density PDF with maximum likelihood.
	
	.. note::
	   You need the `iminuit <https://pypi.python.org/pypi/iminuit>`_ python package before you can use this function. If you don't have that, you need to comment out the `iminuit` related imports in the header of `oPDF.py`.
	
	:param x0: initial value of halo parameters. the interpretation of them depends on the halotype and scales of the tracer's halo. see Tracer.halo of :class:`Tracer` and halo.type, halo.scales of :class:`halo`. 
	:param minuittol: tolerance of minuit to consider convergence. Convergence is defined when the estimated distance to minimum edm<1e-4*minuittol*0.5

	:return: results will be printed on screen.
			 also return minuit result and the minuit minimizer.
			 Please consult the `iminuit <https://pypi.python.org/pypi/iminuit>`_ documentation for the `iminuit` outputs.
			  
	.. note::
	   This is only intended for fitting the Dark Matter density profile to get the NFW parameters.
	   The tracer particle mass should have been properly assigned or adjusted, 
	   so that mP*number_density=physical_density.
	   If you have sampled n particles from the full sample of n0 particles, 
	   remember to adjust the mP of the sample to be mP0*n0/n, so that total mass is conserved.
	'''
	if not self.halo.isNFW():
	  print "Error: not an NFW halo. use Tracer.halo.set_type() to set halo type to NFW before NFW_fit()"
	  
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
	'''plot an image of the particle distribution in proxy-theta space
	:param pars: parameters specifying the potential
	:param proxy: proxy to use for the proxy-theta plot, 'E' or 'L'.
	:param bins: binning in proxy. if an integer, create the input number of bins. If an array, use the array as the bins.
	:param logscale: True or False, whether to bin in logscale or not when bins is an integer.'''
	
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
  datafile=oPDFdir+'/data/mockhalo.hdf5'
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