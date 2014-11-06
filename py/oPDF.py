from math import *
import numpy as np
import ctypes,os,ConfigParser,h5py
from myutils import fmin_gsl,density_of_points,get_extent
from scipy.optimize import fmin, fmin_powell
from scipy.stats import norm,chi2
from iminuit import Minuit
from iminuit.ConsoleFrontend import ConsoleFrontend
import copy

lib=ctypes.CDLL("../liboPDF.so")

if os.uname()[1]=='Medivh':
    rootdir='/work/Projects/DynDistr/'
else:
    rootdir='/gpfs/data/jvbq85/DynDistr/'

MaxNPar=10
Param_t=ctypes.c_double*MaxNPar

class Enum(object):
  def __init__(self, names):
    for number, name in enumerate(names.split()):
      setattr(self, name, number)
HaloTypes=Enum('NFWMC NFWPotsRs NFWRhosRs TMPMC TMPPotScaleRScale CoreRhosRs CorePotsRs')
Estimators=Enum('RBinLike AD MeanPhase MeanPhaseRaw')
VirTypes=Enum('VirTH VirC200 VirB200')

#=============C complex datatypes=====================
class Particle_t(ctypes.Structure):
  _fields_=[('haloid', ctypes.c_int),
		('subid', ctypes.c_int),
		#('strmid', ctypes.c_int),
		('flag', ctypes.c_int),
		('w', ctypes.c_double),
	    ('r', ctypes.c_double),
	    ('K', ctypes.c_double),
	    ('L2', ctypes.c_double),
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
				   ('FlagUseWeight', ctypes.c_int),
				   ('P', Particle_p),
				   ('nbin_r', ctypes.c_int),
				   ('FlagRLogBin', ctypes.c_int),
				   ('RadialCount', ctypes.POINTER(ctypes.c_double)),
				   ('rmin', ctypes.c_double),
				   ('rmax', ctypes.c_double),
				   ('proxybin', ctypes.c_double*2),
				   ('halo', Halo_p),
				   ('nView', ctypes.c_int),
				   ('ViewType', ctypes.c_char),
				   ('Views', Tracer_p)
				  ]
#=======================prototype the library==========================
#globals.h
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
			('units',global_units),
			('virtype',ctypes.c_int)]
  def set_defaults(self):
	'''set default global parameters, 
	including precision, cosmology and units'''
	lib.default_global_pars()
  
  def set_units(self, MassInMsunh=0.73e10, LengthInKpch=0.73, VelInKms=1.):
	'''set system of units,
	specify Mass in Msun/h, Length in kpc/h, Velocity in km/s'''
	lib.set_units(MassInMsunh, LengthInKpch, VelInKms)
Globals=globals_t.in_dll(lib, 'Globals')

#halo.h
class Halo_t(ctypes.Structure):
  #all quantities are physical
  _fields_=[('z', ctypes.c_double),
			('M', ctypes.c_double),
			('c', ctypes.c_double),
			('Rv',ctypes.c_double),
			('Pots',ctypes.c_double),#-4*pi*G*rhos*rs^2, the potential at r=0
			('Rs',ctypes.c_double),
			('Rhos',ctypes.c_double),
			('Ms',ctypes.c_double),#4*pi*rs^3*rhos
			('RScale',ctypes.c_double),#for TMP profile, Rs/Rs0
			('PotScale',ctypes.c_double),#for TMP profile, Pots/Pots0
			('TMPid',ctypes.c_int),#for TMP profile
			('virtype',ctypes.c_int),
			('type',ctypes.c_int)
			]
  def set_type(self, halotype=HaloTypes.NFWMC, virtype=VirTypes.VirC200, redshift=0.):
	'''set the parametrization, virial definition and redshift of halo
	halotype: halo parametrization, one of the HaloTypes members
	virtype: virial definition, one of the VirTypes members
	redshift: redshift of halo'''
	lib.halo_set_type(halotype, virtype, redshift, ctypes.byref(self))
  def set_param(self, pars=[1.,1.]):
	'''set the parameters of the halo
	pars: parameters describing the halo'''
	lib.halo_set_param(Param_t(pars), ctypes.byref(self))
Halo_p=ctypes.POINTER(Halo_t) 
lib.halo_set_type.restype=None
lib.halo_set_type.argtypes=[ctypes.c_int, ctypes.c_int, ctypes.c_double, Halo_p]
lib.halo_set_param.restype=None
lib.halo_set_param.argtypes=[Param_t, Halo_p]
lib.halo_mass.restype=ctypes.c_double
lib.halo_mass.argtypes=[ctypes.c_double, Halo_p]
lib.halo_pot.restype=ctypes.c_double
lib.halo_pot.argtypes=[ctypes.c_double, Halo_p]

#models.h
lib.alloc_integration_space.restype=None
lib.alloc_integration_space.argtypes=[]
lib.free_integration_space.restype=None
lib.free_integration_space.argtypes=[]

openlib=lib.alloc_integration_space
closelib=lib.free_integration_space

lib.like_eval.restype=ctypes.c_double
lib.like_eval.argtypes=[ctypes.c_int, Tracer_p]
lib.nested_views_like.restype=ctypes.c_double
lib.nested_views_like.argtypes=[ctypes.c_int, Tracer_p]
lib.jointLE_like.restype=ctypes.c_double
lib.jointLE_like.argtypes=[ctypes.c_int, ctypes.c_int, ctypes.c_int, Tracer_p]
lib.jointE_like.restype=ctypes.c_double
lib.jointE_like.argtypes=[ctypes.c_int, ctypes.c_int, Tracer_p]
lib.predict_radial_count.restype=None
lib.predict_radial_count.argtypes=[ctypes.POINTER(ctypes.c_double), ctypes.c_int, ctypes.c_int, Tracer_p]



#
lib.create_nested_views.restype=None
lib.create_nested_views.argtypes=[lib.Param_t, ctypes.POINTER(ctypes.c_int), ctypes.c_char_p, Tracer_p];
lib.nested_views_FChi2.restype=ctypes.c_double
lib.nested_views_FChi2.argtypes=[lib.Param_t, ctypes.c_int, Tracer_p]

#halos
lib.define_halo.restype=None
lib.define_halo.argtypes=[lib.Param_t]
lib.halo_pot.restype=ctypes.c_double
lib.halo_pot.argtypes=[ctypes.c_double]
lib.halo_mass.restype=ctypes.c_double
lib.halo_mass.argtypes=[ctypes.c_double]
lib.comoving_virial_radius.argtypes=[ctypes.c_double, ctypes.c_double, ctypes.c_int]
lib.comoving_virial_radius.restype=ctypes.c_double
lib.decode_NFWprof.argtypes=[ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_int, ctypes.POINTER(NFWHalo_t)]
lib.decode_NFWprof.restype=None
lib.decode_TemplateProf.argtypes=[ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_int, ctypes.POINTER(NFWHalo_t)]
lib.decode_TemplateProf.restype=None
lib.NFW_mass.restype=ctypes.c_double
lib.NFW_mass.argtypes=[ctypes.c_double]
lib.NFW_like.restype=ctypes.c_double
lib.NFW_like.argtypes=[lib.Param_t, Tracer_p]

#tracers
lib.load_tracer_particles.restype=None
lib.load_tracer_particles.argtypes=[ctypes.c_char_p, Tracer_p]
lib.cut_tracer_particles.restype=None
lib.cut_tracer_particles.argtypes=[Tracer_p, ctypes.c_double, ctypes.c_double]
#extern void cut_tracer(Tracer_t *Sample, double rmin, double rmax);
lib.shuffle_tracer_particles.restype=None
lib.shuffle_tracer_particles.argtypes=[ctypes.c_ulong, Tracer_p]
lib.squeeze_tracer_particles.restype=None
lib.squeeze_tracer_particles.argtypes=[Tracer_p]
#extern void squeeze_tracer(Tracer_t *Sample);
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
#extern void count_tracer_radial(Tracer_t *Sample, int nbin);
lib.free_tracer_rcounts.restype=None
lib.free_tracer_rcounts.argtypes=[Tracer_p]
lib.sort_part_flag.restype=None
lib.sort_part_flag.argtypes=[Particle_p, ctypes.c_int]
#extern void sort_part_flag(Particle_t *P, int nP);
lib.sort_part_E.restype=None
lib.sort_part_E.argtypes=[Particle_p, ctypes.c_int]
#extern void sort_part_L(Particle_t *P, int nP);
lib.sort_part_L.restype=None
lib.sort_part_L.argtypes=[Particle_p, ctypes.c_int]
#extern void sort_part_E(Particle_t *P, int nP);
lib.sort_part_R.restype=None
lib.sort_part_R.argtypes=[Particle_p, ctypes.c_int]
lib.create_tracer_views.restype=None
lib.create_tracer_views.argtypes=[Tracer_p, ctypes.c_int, ctypes.c_char]
lib.free_tracer_views.restype=None
lib.free_tracer_views.argtypes=[Tracer_p]
lib.init_tracer.restype=None
lib.init_tracer.argtypes=[Tracer_p]
#extern void init_tracer(Tracer_t *Sample);
lib.make_sample.restype=None
lib.make_sample.argtypes=[ctypes.c_int, ctypes.c_int, Tracer_p, Tracer_p]
#extern void make_sample(int offset, int size, Tracer_t *Sample, Tracer_t *FullSample);
lib.free_tracer.restype=None
lib.free_tracer.argtypes=[Tracer_p]
#extern void free_tracer(Tracer_t *Sample);
lib.init_potential_spline.restype=None
lib.init_potential_spline.argtypes=[]
lib.free_potential_spline.restype=None
lib.free_potential_spline.argtypes=[]
lib.eval_potential_spline.restype=ctypes.c_double
lib.eval_potential_spline.argtypes=[ctypes.c_double]

lib.mock_stars.argtypes=[ctypes.c_char, ctypes.c_int, Tracer_p]
lib.mock_stars.restype=None
lib.save_mockstars.argtypes=[ctypes.c_char, Tracer_p]
lib.save_mockstars.restype=None

#wenting's lib
prototype=ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double)
lib.dataprob=prototype(("dataprob",lib))
PARTYPEXV=ctypes.c_double*6
prototype=ctypes.CFUNCTYPE(ctypes.c_double, PARTYPEXV)
lib.dataprob6d=prototype(("dataprob6d",lib))

#elike=lambda m,c,estimator: likefunc(PARTYPE(m,c),estimator)
#freeze_energy=lambda m,c: freezefunc(PARTYPE(m,c))
#elike2=lambda m,c,estimator: likefunc2(PARTYPE(m,c),estimator) #freeze_and_like()
#====================python wrapper classes=========================================

class NFWHalo(object):
  '''halo related functions and variables''' 
  def __init__(self):
	self.halo=NFWHalo_t.in_dll(lib, 'Halo')
	self.M0=ctypes.c_double.in_dll(lib, 'HaloM0')
	self.C0=ctypes.c_double.in_dll(lib, 'HaloC0')
	self.Rhos0=ctypes.c_double.in_dll(lib, 'HaloRhos0')
	self.Rs0=ctypes.c_double.in_dll(lib, 'HaloRs0')
	self.ProfID=ctypes.c_int.in_dll(lib, 'HaloProfID')
	self.FitType=ctypes.c_int.in_dll(lib, 'HaloFitType')
	
  def define_halo(self, pars):
	lib.define_halo(lib.Param_t(*pars))
	
  def pot(self,r):
	return lib.halo_pot(r)
  
  def comoving_virial_radius(m,z,virtype='c200'):
	vt={'sc':0,'c200':1,'b200':2}
	return lib.comoving_virial_radius(m,z,vt[virtype])
  
  def mass(self, r):
	return lib.halo_mass(r)
 	
class Tracer(Tracer_t):
  def __init__(self, halo=None, **newoptions):
	'''if called with no par, create an empty tracer; otherwise load according to config
	configuration according to get_config(halo). newoptions will overide get_config.'''
	Tracer_t.__init__(self)
	#self.load_tracer_particles()
	#self.rmin=min(self.data['r'])
	#self.rmax=xxx
	#if rmin, then self.rmin=xx; self.cut_radial()
	#self.shuffle()
	self.nP=0  #init state
	self.nView=0
	self.FlagUseWeight=0
	#self._as_parameter_=ctypes.byref(self) #this makes it possible to directly pass this as pointer arg
	self._pointer=ctypes.byref(self) #or ctypes.cast(ctypes.pointer(self), Tracer_p)? or ctypes.pointer(self)
	if halo!=None or newoptions!={}:
	  options=lib.get_config(halo)
	  options.update(newoptions)
	  print options
	  for k,v in options.iteritems():
		os.environ[k]=str(v)
	  self.auto_load()
	  
  def clean(self):
	''' never defined it as __del__ and rely on the garbage collector.
	it's dangerous. gc may never call your __del__.
	call it yourself.'''
	#print "cleaning Tracer ", id(self)
	lib.free_tracer(self._pointer)
	#print self.nView

  def __enter__(self):
	return self
  
  def __exit__(self, type, value, traceback):
	self.clean()
	
  def __update_array(self):
	'''this should be called whenever the pointer self.P has changed
	to point the numpy array to the new memory'''
	Parr=(Particle_t*self.nP).from_address(ctypes.addressof(self.P.contents)) #struct array
	self.data=np.frombuffer(Parr, np.dtype(Parr))[0] #the numpy array
  
  def print_data(self, i=0):
	print self.nP, '%0x'%ctypes.addressof(self.P.contents)
	print self.P[i].x[0], self.P[i].r
	print self.data[i]['x'][0], self.data[i]['r']
	print '-----------'
	lib.print_tracer_particle(self._pointer, i)
	print '============='

  def wenting_like_conditional(self, pars=[1,1]):
	"wenting's likelihood for the mocks, with other parameters conditioned at true values"
	lnL=lib.wenting_like(lib.Param_t(pars[0],pars[1],1,1,1,1), self._pointer)
	#print pars, lnL
	return lnL
  
  def wenting_like_marginal(self, pars=[1,1]):
	'''wenting's likelihood for the mocks, with other parameters marginalized.
	this is too slow to use...'''
	like=lambda c,d,e,f: -lib.wenting_like(lib.Param_t(pars[0],pars[1],c,d,e,f), self._pointer)
	m=Minuit(like, c=1,d=1,e=1,f=1, fix_d=False, fix_e=0, fix_f=0, limit_c=[-0.7,1.4], limit_d=[0.1,2], limit_e=[0.1,10], limit_f=[0.1,10], print_level=3, pedantic=False, errordef=1, frontend=ConsoleFrontend())
	#m.set_strategy(0)
	m.tol=10   #default convergence edm<1e-4*tol*errordef, but we do not need that high accuracy
	m.migrad()
	print pars, -m.fval
	return -m.fval
  
  def freeze_energy(self, pars=[1,1]):
	lib.freeze_energy(lib.Param_t(*pars), self._pointer)

  def likelihood(self, pars=[1,1], estimator=10):
	return lib.likelihood(lib.Param_t(*pars), estimator, self._pointer)
  
  def like_init(self, pars=[1,1], estimator=10):
	lib.like_init(lib.Param_t(*pars), estimator, self._pointer)

  def like_eval(self, pars=[1,1], estimator=10):
	return lib.like_eval(lib.Param_t(*pars), estimator, self._pointer)
	
  def freeze_and_like(self, pars=[1,1], estimator=10):
	y=lib.freeze_and_like(lib.Param_t(*pars), estimator, self._pointer)
	#print pars, y
	return y
  
  def jointE_FChi2(self, pars=[1,1], estimator=10, nbinE=10):
	return lib.jointE_FChi2(lib.Param_t(*pars), estimator, nbinE, self._pointer)
	
  def jointLE_FChi2(self, pars=[1,1], estimator=10, nbinL=10, nbinE=10):
	return lib.jointLE_FChi2(lib.Param_t(*pars), estimator, nbinL, nbinE, self._pointer)

  def create_nested_views(self, pars=[1,1], viewtypes='EL', nbins=[10,10]):
	try:
	  nbins=list(nbins)
	except:
	  nbins=[nbins]
	lib.create_nested_views(lib.Param_t(*pars), (ctypes.c_int*(len(nbins)+1))(*nbins), ctypes.c_char_p(viewtypes), self._pointer)
  
  def nested_views_Chi2(self, pars=[1,1], estimator=10):
	return lib.nested_views_Chi2(lib.Param_t(*pars), estimator, self._pointer)
  
  def nested_views_FChi2(self, pars=[1,1], estimator=10):
	return lib.nested_views_FChi2(lib.Param_t(*pars), estimator, self._pointer)
  
  def joint_FChi2(self, pars, estimator, viewtypes, nbins):
	'''nestviews and like'''
	self.create_nested_views(pars, viewtypes, nbins)
	return self.nested_views_FChi2(pars, estimator)
  
  def TSprof(self, pars=[1,1], estimator=10, viewtypes='L', nbins=100):
	if viewtypes=='L2':
	  viewtypes='L'
	with self.copy(0,0) as S: #make a copy to avoid messing up the weights
	  S.joint_FChi2(pars, estimator, viewtypes, nbins)
	  ts=[S.Views[i].lnL for i in xrange(S.nView)] #these are like_eval, not chi2
	ts.append(np.nan)
	ts=np.array(ts)
	x,p=self.gen_bin(viewtypes, nbins, equalcount=True)
	return ts,x
  
  def TSprofCum(self, pars=[1,1], estimator=8, proxy='r', bins=100):
	'''cumulative TS profile
	reuturn bin edges, ts, counts'''
	self.freeze_energy(pars)
	self.sort(proxy)
	if proxy=='E':
	  self.data[...]=self.data[-1::-1] #reverse the order
	#if not isinstance(bins, int): 
	  #if bins[0]>self.data[proxy][0]:#always start from the min value of the proxy
		#tmp,bins=bins,[self.data[proxy][0]-1]
		#bins.extend(list(tmp))
	n,bins=np.histogram(self.data[proxy],bins)
	if proxy=='E':
	  n=n[-1::-1].cumsum()
	else:
	  n=n.cumsum()
	ts=[]
	if proxy!='r':
	  self.like_init(pars,estimator)
	for i,x in enumerate(bins[1:]):
	  if n[i]==0:
		y=np.nan
	  else:
		with self.copy(0,n[i]) as s:
		  if proxy!='r':
			y=s.like_eval(pars, estimator)
		  else:
			s.rmax=x
			y=s.likelihood(pars,estimator)
	  ts.append(y)
	if proxy=='E':
	  ts=ts[-1::-1]
	return np.array(bins),np.array(ts),np.array(n)

  def scan_like(self, estimator, x, y):
	'''scan a likelihood surface to be used for contour plots as contour(x,y,z)'''
	like=lambda x: self.freeze_and_like(x, estimator)
	z=[like([m,c]) for m in x for c in y]
	return x,y,np.array(z).reshape([len(y),len(x)], order="F")
  
  def scan_Chi2(self, estimator, x, y):
	'''scan a likelihood surface to be used for contour plots as contour(x,y,z)'''
	like=lambda x: lib.like_to_chi2(self.freeze_and_like(x, estimator), estimator) #the real likelihood prob translated to chi-square values
	z=[like([m,c]) for m in x for c in y]
	return x,y,np.array(z).reshape([len(y),len(x)], order="F")
  
  def phase_density(self, proxy='E', bins=100, method='hist', logscale=False):
	'''estimate density in proxy-theta space'''
	if logscale:
	  f=self.data[proxy]>0
	  data=np.array((np.log10(self.data[proxy][f]), self.data['theta'][f]))
	  w=self.data['w'][f]
	else:
	  data=np.array((self.data[proxy], self.data['theta']))
	  w=self.data['w']
	if self.FlagUseWeight:
	  X,Y,Z=density_of_points(data, bins=bins, method=method, weights=w)
	else:
	  X,Y,Z=density_of_points(data, bins=bins, method=method)
	extent=get_extent(X,Y)
	return X,Y,Z,extent,data
  
  def fmin_FixBinIter(self, estimator,proxy,nbins,x0=[1,1],itertol=0.01, maxiter=50, xtol=1e-3, ftolabs=0.1, ftolrel=1e-2):
	''' itertol: tolerance for iteration
	xtol,ftolabs,ftolrel: tolerance for each fmin() step
	'''
	x=x0
	x0=[x[0]+1,x[1]]
	iter=0
	while abs(x0[0]-x[0])>itertol or abs(x0[1]-x[1])>itertol:
	  x0=x
	  self.create_nested_views(x0, proxy, nbins)
	  ftol=min(ftolabs/abs(self.nested_views_FChi2(x0, estimator)), ftolrel)
	  result=fmin_powell(self.nested_views_FChi2, x0, args=(estimator,), xtol=xtol, ftol=ftol, maxiter=1000, maxfun=5000, full_output=True)
	  x=result[0]
	  print result[0], result[1], result[-1]
	  iter+=1
	  if iter>maxiter:
		print "Warning: maxiter=%d reached in fmin_FixBinIter with"%maxiter,proxy,nbins
		result=list(result)
		result[-1]=3 #top level maxiter has reached.
		return tuple(result)
	return result #converged
  
  def fmin_jointLE(self, estimator, nbinL, nbinE, x0=[1,1], xtol=1e-3, ftolabs=0.01, ftolrel=1e-2):
	ftol=min(ftolabs/abs(self.jointLE_FChi2(x0,estimator,nbinL,nbinE)), ftolrel)
	x=fmin_powell(self.jointLE_FChi2, x0, args=(estimator, nbinL, nbinE), xtol=xtol, ftol=ftol, maxiter=500, maxfun=1000, full_output=True)
	print '\t ', x[0]
	return x
  
  def fmin_like(self, estimator, x0=[1,1], xtol=1e-3, ftolabs=0.01, ftolrel=1e-2):
	like=lambda x: lib.like_to_chi2(self.freeze_and_like(x, estimator), estimator) #the real likelihood prob
	ftol=min(ftolabs/abs(like(x0)), ftolrel)
	x=fmin_powell(like, x0, xtol=xtol, ftol=ftol, maxiter=500, maxfun=1000, full_output=True)
	print '\t ', x[0]
	return x
  
  def fmin_dist(self, estimator, x0=[1,1],xtol=1e-3, ftolabs=0.01, ftolrel=1e-2):
	like=lambda x: -self.freeze_and_like(x, estimator) #distance
	ftol=min(ftolabs/abs(like(x0)), ftolrel)
	x=fmin_powell(like, x0, xtol=xtol, ftol=ftol, maxiter=500, maxfun=1000, full_output=True)
	print '\t ', x[0]
	return x
  
  def gfmin_FixBinIter(self, estimator,proxy,nbins,x0=[1,1],itertol=0.01, maxiter=50, xtol=1e-3, ftolabs=0.01):
	''' itertol: tolerance for iteration
	xtol,ftolabs,ftolrel: tolerance for each fmin_gsl() step
	'''
	x=x0
	x0=[x[0]+1,x[1]]
	iter=0
	while abs(x0[0]-x[0])>itertol or abs(x0[1]-x[1])>itertol:
	  x0=x
	  self.create_nested_views(x0, proxy, nbins)
	  result=fmin_gsl(self.nested_views_FChi2, x0, args=(estimator,), xtol=xtol, ftolabs=ftolabs, maxiter=1000,  full_output=True)
	  x=result[0]
	  print result[0], result[1], result[-1]
	  iter+=1
	  if iter>maxiter:
		print "Warning: maxiter=%d reached in gfmin_FixBinIter with"%maxiter,proxy,nbins
		result=list(result)
		result[-1]=3 #top level maxiter has reached.
		return tuple(result)
	return result #converged
  
  def gfmin_jointLE(self, estimator, nbinL, nbinE, x0=[1,1], xtol=1e-3, ftolabs=0.01):
	return fmin_gsl(self.jointLE_FChi2, x0, args=(estimator, nbinL, nbinE), xtol=xtol, ftolabs=ftolabs, maxiter=500, full_output=True)
  
  def gfmin_like(self, estimator, x0=[1,1], xtol=1e-3, ftolabs=0.01):
	like=lambda x: lib.like_to_chi2(self.freeze_and_like(x, estimator), estimator) #the real likelihood prob
	return fmin_gsl(like, x0, xtol=xtol, ftolabs=ftolabs, maxiter=500, full_output=True)
  
  def gfmin_dist(self, estimator, x0=[1,1], xtol=1e-3, ftolabs=0.01):
	like=lambda x: -self.freeze_and_like(x, estimator) #distance
	return fmin_gsl(like, x0, xtol=xtol, ftolabs=ftolabs, maxiter=500, full_output=True)
  
  def minuit_like(self, estimator, x0=[1,1], minuittol=1):
	"""too difficult for minuit to work. just use this to estimate the error matrix? still just fantasy. discard it.
	set a huge tol, say, 1e10, to avoid walking away"""
	like=lambda m,c: lib.like_to_chi2(self.freeze_and_like([m,c], estimator),estimator)
	#profilelikelihood ratio error-def: chi-square1
	m=Minuit(like, m=x0[0], c=x0[1], print_level=0, pedantic=False, error_m=0.1, error_c=0.1, errordef=1, limit_m=[0.1,10],limit_c=[0.1,10],frontend=ConsoleFrontend())
	m.tol=minuittol   #default convergence edm<1e-4*tol*errordef, but we do not need that high accuracy
	result=m.migrad()
	m.print_fmin()
	m.print_matrix()
	return result,m
  
  def minuit_jointLE(self, estimator, nbinL, nbinE, x0=[1,1], minuittol=1):
	"""too difficult for minuit to work. just use this to estimate the error matrix.
	set a huge tol, say, 1e10, to avoid walking away"""
	like=lambda m,c: self.jointLE_FChi2([m,c], estimator, nbinL, nbinE)
	#profilelikelihood ratio error-def: chi-square1
	m=Minuit(like, m=x0[0],c=x0[1], print_level=0, pedantic=False, error_m=0.1, error_c=0.1, errordef=1, limit_m=[0.1,10],limit_c=[0.1,10],frontend=ConsoleFrontend())
	m.tol=minuittol   #default convergence edm<1e-4*tol*errordef, but we do not need that high accuracy
	result=m.migrad()
	m.print_fmin()
	m.print_matrix()
	return result,m
  
  def minuit_FixBinIter(self, estimator,proxy,nbins,x0=[1,1],tol=0.01, minuittol=1):
	'''set tol and minuittol to 1e10 to only get error and avoid walking away'''
	like= lambda m,c: self.nested_views_FChi2([m,c], estimator) #note: nested_views_Chi2() does not work here.
	m=Minuit(like, m=x[0], c=x[1], print_level=0, pedantic=False, error_m=0.1, error_c=0.1, errordef=1, limit_m=[0.1,10],limit_c=[0.1,10],frontend=ConsoleFrontend())
	m.tol=minuittol
	x=x0
	x0=[x[0]+1,x[1]]
	while abs(x0[0]-x[0])>tol or abs(x0[1]-x[1])>tol:
	  x0=x
	  self.create_nested_views(x0, proxy, nbins)
	  result=m.migrad()
	  x=[m.values['m'],m.values['c']]
	  m.print_fmin()
	m.print_matrix()
	return result,m
  
  def NFW_like(self, pars=[1,1]):
	return lib.NFW_like(lib.Param_t(*pars), self._pointer)
	
  def minuit_NFWlike(self, x0=[1,1], minuittol=1):
	''' to fit a NFW density PDF with ML '''
	like=lambda m,c: -self.NFW_like([m,c])
	#profilelikelihood ratio error-def: chi-square1
	m=Minuit(like, m=x0[0], c=x0[1], print_level=0, pedantic=False, error_m=0.1, error_c=0.1, errordef=0.5,frontend=ConsoleFrontend())
	m.tol=minuittol   #default convergence edm<1e-4*tol*errordef, but we do not need that high accuracy
	result=m.migrad()
	m.print_fmin()
	m.print_matrix()
	return result,m
  
  def predict_radial_count(self, nbin=100, logscale=True):
	n=np.empty(nbin,dtype='f8')
	lib.predict_radial_count(n.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), nbin, logscale, self._pointer)
	#lib.predict_radial_count(n.ctypes.data, nbin) //also work?
	return n
  
  def load(self, infile):
	lib.load_tracer_particles(ctypes.c_char_p(infile), self._pointer)
	self.__update_array()
	
  def copy(self, offset=0, size=0):
	newsample=Tracer()
	lib.copy_tracer_particles(int(offset), int(size), newsample._pointer, self._pointer)
	newsample.__update_array()
	return newsample
  
  def select(self, flags):
	'''select particles according to flags array, to create a subsample'''
	sample=self.copy(0,0)
	sample.data['flag']=flags
	sample.squeeze() # __update_array() is called automatically
	return sample
	
  def create_views(self, n=10, proxy='L'):
	lib.create_tracer_views(self._pointer, n, proxy)
	
  def destroy_views(self):
	lib.free_tracer_views(self._pointer)

  def radial_cut(self, rmin, rmax):
	lib.cut_tracer_particles(self._pointer, rmin, rmax)
	self.__update_array()
	
  def shuffle(self, seed):
	lib.shuffle_tracer_particles(ctypes.c_ulong(seed), self._pointer)
	
  def resample(self, seed):
	newsample=Tracer()
	lib.resample_tracer_particles(seed, newsample._pointer, self._pointer)
	newsample.__update_array()
	return newsample

  def squeeze(self):
	lib.squeeze_tracer_particles(self._pointer)
	self.__update_array()
	
  def radial_count(self, nbin=None, logscale=True):
	if nbin is None:
	  nbin=lib.NumRadialCountBin.value
	lib.count_tracer_radial(self._pointer, nbin, logscale)
	  
  def sort(self, proxy, offset=0,n=0):
	if n==0:
	  n=self.nP
	sort_func={'flag': lib.sort_part_flag, 'E': lib.sort_part_E, 'L2': lib.sort_part_L, 'r': lib.sort_part_R}
	sort_func[proxy](ctypes.byref(Particle_t.from_buffer(self.P[offset])), n)

  def auto_load(self):
	'''high level init; it does the following:
	1) config internal parameters according to env; 2)load the data; 3)radial_cut 4)shuffle 5)count_radial
	also tries to free() if already loaded'''
	if self.nP>0:
	  lib.free_tracer(self._pointer)
	lib.init_tracer(self._pointer)
	self.__update_array()
	
  def sample(self, sampleid=0):
	'''high level sampling; 
	copy() and radial_count()'''
	newsample=Tracer()
	size=lib.SubSampleSize.value
	lib.make_sample(sampleid*size, size, newsample._pointer, self._pointer)
	newsample.__update_array()
	return newsample
	
  def gen_bin(self,bintype,nbin=30,logscale=True, equalcount=False):
	if bintype=='L':
	  bintype='L2'
	proxy=np.copy(self.data[bintype])
	n=nbin+1
	if equalcount:
		x=np.array(np.percentile(proxy,list(np.linspace(0,100,n))))
	else:  
		if logscale:  
			x=np.logspace(np.log10(proxy[proxy>0].min()),np.log10(proxy.max()),n)
		else:
			x=np.linspace(proxy.min(),proxy.max(),n)
	return x,proxy
  
  def assign_TSmap(self, TS, x):
	 #TODO
	 pass

  def filter_TSmap(self, TS, x, lim):
	'''select particles with TS in the range [lim[0],lim[1]], 
	where TS is a map with pixel boundaries in x, 
	and x is an ordereddict
	return: update P['flag'] in place.'''
	p=((TS>lim[0])*(TS<lim[1])).nonzero()
	b=x.keys()
	xx=x.values()
	f=np.zeros(self.nP)
	for pix in zip(p[0],p[1]):
		f=f+(self.P[b[0]]>xx[0][pix[0]])*(self.P[b[0]]<xx[0][pix[0]+1])*(self.P[b[1]]>xx[1][pix[1]])*(self.P[b[1]]<xx[1][pix[1]+1])
	f=f>0
	self.P['flag'][:]=f
	print "Found %d pixels, %d particles"%(p[0].shape[0],np.sum(self.P['flag']))

lib.open()

if __name__=="__main__":
  #lib.open() # common initialization
  FullSample=Tracer('AqA4',DynSIZE=100)
  Sample=FullSample.sample()
  FullSample.print_data(10)
  Sample.print_data(10)
  Sample.data[1]['r']=2
  Sample.rmin=10
  Sample.print_data(1)
  Sample.sort('L2')
  Sample.print_data(1)
  print Sample.jointLE_like([1,1])
  FullSample.clean()
  Sample.clean()
  
  #good practice: use with!
  with Tracer('Mock') as FullSample:
	with FullSample.sample() as Sample:
	  FullSample.print_data()
	  Sample.print_data()
	  Sample.data[1]['r']=2
	  Sample.rmin=10
	  Sample.print_data()

  #lib.close()