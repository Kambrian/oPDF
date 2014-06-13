from math import *
import numpy as np
import ctypes,os,ConfigParser,h5py

#=============C complex datatypes=====================
class Particle_t(ctypes.Structure):
  _fields_=[('flag', ctypes.c_int),
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
Tracer_t._fields_=[('nP', ctypes.c_int),
				   ('P', Particle_p),
				   ('nbin_r', ctypes.c_int),
				   ('RadialCount', ctypes.POINTER(ctypes.c_int)),
				   ('rmin', ctypes.c_double),
				   ('rmax', ctypes.c_double),
				   ('nView', ctypes.c_int),
				   ('ViewType', ctypes.c_char),
				   ('Views', Tracer_p)
				  ]

class NFWHalo_t(ctypes.Structure):
  """all properties are physical"""
  _fields_=[('z', ctypes.c_double),
	    ('M', ctypes.c_double),
	    ('c', ctypes.c_double),
	    ('Rv', ctypes.c_double), 
	    ('Rs', ctypes.c_double),
	    ('Rhos', ctypes.c_double),
	    ('Pots', ctypes.c_double),#-4*pi*G*rhos*rs^2, the potential at r=0
	    ('virtype', ctypes.c_int)
	    ]  

#=======================load the library==========================
lib=ctypes.CDLL("./libdyn.so")
#general
lib.MaxNPar=10
lib.ParType=ctypes.c_double*lib.MaxNPar
lib.MODEL_TOL_BIN=ctypes.c_double.in_dll(lib,'MODEL_TOL_BIN')
lib.MODEL_TOL_BIN_ABS=ctypes.c_double.in_dll(lib,'MODEL_TOL_BIN_ABS')
lib.MODEL_TOL_REL=ctypes.c_double.in_dll(lib,'MODEL_TOL_REL')
lib.SubSampleSize=ctypes.c_int.in_dll(lib,'SubSampleSize')
lib.NumRadialCountBin=ctypes.c_int.in_dll(lib,'NumRadialCountBin')
lib.alloc_integration_space.restype=None
lib.alloc_integration_space.argtypes=[]
lib.free_integration_space.restype=None
lib.free_integration_space.argtypes=[]
lib.like_to_chi2.restype=ctypes.c_double
lib.like_to_chi2.argtypes=[ctypes.c_double, ctypes.c_int]
lib.NameList={0:'f(E,L)',4:'RBin',8:'AD',9:'Resultant',10:'Mean',11:'KS',12:'Kuiper',13:'CosMean', 14:'RawMean'}
if os.uname()[1]=='Medivh':
    lib.rootdir='/work/Projects/DynDistr/'
else:
    lib.rootdir='/gpfs/data/jvbq85/DynDistr/'
    
lib.open=lib.alloc_integration_space
lib.close=lib.free_integration_space

def get_config(halo):
  if halo==None:
	return {}
  c=ConfigParser.ConfigParser()
  c.optionxform=str
  c.read('DataFiles.cfg')
  try:
	options=dict(c.items(halo))
  except:
	options=dict(c.defaults())
	print "Warning: no config for %s; using defaults."%halo
  return options    

lib.get_config=get_config

#models
lib.like_init.restype=None
lib.like_init.argtypes=[lib.ParType, ctypes.c_int, Tracer_p]
lib.like_eval.restype=ctypes.c_double
lib.like_eval.argtypes=[lib.ParType, ctypes.c_int, Tracer_p]
lib.likelihood.restype=ctypes.c_double
lib.likelihood.argtypes=[lib.ParType, ctypes.c_int, Tracer_p]
lib.freeze_energy.restype=None
lib.freeze_energy.argtypes=[lib.ParType, Tracer_p]
lib.freeze_and_like.restype=ctypes.c_double
lib.freeze_and_like.argtypes=[lib.ParType, ctypes.c_int, Tracer_p]
lib.jointLE_like.restype=ctypes.c_double
lib.jointLE_like.argtypes=[lib.ParType, ctypes.c_int, ctypes.c_int, ctypes.c_int, Tracer_p]
#extern double jointE_like(double pars[], int estimator, int nbin, Tracer_t *Sample);
lib.jointE_like.restype=ctypes.c_double
lib.jointE_like.argtypes=[lib.ParType, ctypes.c_int, ctypes.c_int, Tracer_p]

lib.predict_radial_count.restype=None
lib.predict_radial_count.argtypes=[ctypes.POINTER(ctypes.c_double), ctypes.c_int, Tracer_p]
#halos
lib.define_halo.restype=None
lib.define_halo.argtypes=[lib.ParType]
lib.halo_pot.restype=None
lib.halo_pot.argtypes=[ctypes.c_double]
lib.comoving_virial_radius.argtypes=[ctypes.c_double, ctypes.c_double, ctypes.c_int]
lib.comoving_virial_radius.restype=ctypes.c_double

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
lib.count_tracer_radial.argtypes=[Tracer_p, ctypes.c_int]
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
	halo=NFWHalo_t.in_dll(lib, 'Halo')
	M0=ctypes.c_double.in_dll(lib, 'M0')
	C0=ctypes.c_double.in_dll(lib, 'C0')
	Rhos0=ctypes.c_double.in_dll(lib, 'Rhos0')
	Rs0=ctypes.c_double.in_dll(lib, 'Rs0')
	
  def define_halo(self, pars):
	lib.define_halo(lib.ParType(*pars))
	
  def pot(self,r):
	return lib.halo_pot(r)
  
  def comoving_virial_radius(m,z,virtype='c200'):
	vt={'sc':0,'c200':1,'b200':2}
	return lib.comoving_virial_radius(m,z,vt[virtype])
  
class Tracer(Tracer_t):
  def __init__(self, halo=None, **newoptions):
	'''if called with no par, create an empty tracer; otherwise load according to config
	configuration according to get_config(halo). newoptions will overide get_config.'''
	Tracer_t.__init__(self)
	self.nP=0  #init state
	self.nView=0
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
	
  def freeze_energy(self, pars=[1,1]):
	lib.freeze_energy(lib.ParType(*pars), self._pointer)

  def likelihood(self, pars=[1,1], estimator=10):
	return lib.likelihood(lib.ParType(*pars), estimator, self._pointer)
  
  def freeze_and_like(self, pars=[1,1], estimator=10):
	return lib.freeze_and_like(lib.ParType(*pars), estimator, self._pointer)
  
  def jointE_like(self, pars=[1,1], estimator=10, nbinE=10):
	return lib.jointE_like(lib.ParType(*pars), estimator, nbinE, self._pointer)
	
  def jointLE_like(self, pars=[1,1], estimator=10, nbinL=10, nbinE=10):
	return lib.jointLE_like(lib.ParType(*pars), estimator, nbinL, nbinE, self._pointer)
		
  def predict_radial_count(self, nbin=100):
	n=np.empty(nbin,dtype='f8')
	lib.predict_radial_count(n.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), nbin)
	#lib.predict_radial_count(n.ctypes.data, nbin) //also work?
	return n
  
  def load(self, infile):
	lib.load_tracer_particles(ctypes.c_char_p(infile), self._pointer)
	self.__update_array()
	
  def copy(self, offset=0, size=0):
	newsample=Tracer()
	lib.copy_tracer_particles(offset, size, newsample._pointer, self._pointer)
	newsample.__update_array()
	return newsample
  
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
	
  def radial_count(self, nbin=None):
	if nbin is None:
	  nbin=lib.NumRadialCountBin
	lib.count_tracer_radial(self._pointer, nbin)
	  
  def sort(self, proxy, offset=0,n=0):
	if n==0:
	  n=self.nP
	sort_func={'flag': lib.sort_part_flag, 'E': lib.sort_part_E, 'L2': lib.sort_part_L}
	sort_func[proxy](ctypes.byref(Particle_t.from_buffer(Sample.P[offset])), n)

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

class SubData(object):
  def __init__(self, infile):
    f=h5py.File(infile,'r')
    self.x=f['/x'][...]
    self.v=f['/v'][...]
    self.m=f['/m'][...]
    f.close()
    self.r=np.sqrt(np.sum(self.x**2,1))
    self.vr=np.sum(self.v*self.x,1)/self.r
    self.vr[np.isnan(self.vr)]=0
    self.K=np.sum(self.v**2,1)/2.
    self.L2=np.sum(np.cross(self.x, self.v)**2,1)
    self.Nsub=self.x.shape[0]
    
  def eval_energy(self):
	"""calc binding energy for each sub. the halo potential must be initialized before calling this func,
	via freeze_energy() or define_halo()"""
	self.E=np.array([-k-lib.halo_pot(r) for k,r in zip(self.K,self.r)])
  
if __name__=="__main__":
  lib.open() # common initialization
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

  lib.close()