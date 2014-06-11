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
  _fields_=[('nP', ctypes.c_int),
			('nbin_r', ctypes.c_int),
			('RadialCount', ctypes.POINTER(ctypes.c_int)),
			('rmin', ctypes.c_double),
			('rmax', ctypes.c_double),
			('P', ctypes.POINTER(Particle_t))
			]
Tracer_p=ctypes.POINTER(Tracer_t)

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
lib.SUBSAMPLE_SIZE=ctypes.c_int.in_dll(lib,'SUBSAMPLE_SIZE')
lib.NumRadialCountBin=ctypes.c_int.in_dll(lib,'NumRadialCountBin')
lib.alloc_integration_space.restype=None
lib.alloc_integration_space.argtypes=[]
lib.free_integration_space.restype=None
lib.free_integration_space.argtypes=[]
lib.NameList={0:'f(E,L)',4:'RBin',8:'AD',9:'Resultant',10:'Mean',11:'KS',12:'Kuiper',13:'CosMean', 14:'RawMean'}

#models
lib.like_init.restype=None
lib.like_init.argtypes=[lib.ParType, ctypes.c_int, Tracer_p]
lib.like_eval.restype=ctypes.c_double
lib.like_eval.argtypes=[lib.ParType, ctypes.c_int, Tracer_p]
lib.likelihood.restype=ctypes.c_double
lib.likelihood.argtypes=[lib.ParType, ctypes.c_int, Tracer_p]
lib.freeze_energy.restype=None
lib.freeze_energy.argtypes=[lib.ParType]
lib.freeze_and_like.restype=ctypes.c_double
lib.freeze_and_like.argtypes=[lib.ParType, ctypes.c_int, Tracer_p]

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
lib.load_tracer.restype=None
lib.load_tracer.argtypes=[ctypes.c_char_p, Tracer_p]
lib.shuffle_tracer.restype=None
lib.shuffle_tracer.argtypes=[ctypes.c_ulong, Tracer_p]
lib.resample_tracer.restype=None
lib.resample_tracer.argtypes=[ctypes.c_ulong, Tracer_p, Tracer_p]
lib.copy_tracer.restype=None
lib.copy_tracer.argtypes=[ctypes.c_int, ctypes.c_int, Tracer_p, Tracer_p]
#extern void copy_tracer(int offset, int sample_size, Tracer_t *Sample, Tracer_t *FullSample);
lib.squeeze_tracer.restype=None
lib.squeeze_tracer.argtypes=[Tracer_p]
#extern void squeeze_tracer(Tracer_t *Sample);
lib.free_tracer.restype=None
lib.free_tracer.argtypes=[Tracer_p]
#extern void free_tracer(Tracer_t *Sample);
lib.print_tracer.restype=None
lib.print_tracer.argtypes=[Tracer_p]
#extern void print_tracer(Tracer_t *Sample);
lib.sort_tracer_flag.restype=None
lib.sort_tracer_flag.argtypes=[Tracer_p]
#extern void sort_tracer_flag(Tracer_t *Sample);
lib.cut_tracer.restype=None
lib.cut_tracer.argtypes=[Tracer_p, ctypes.c_double, ctypes.c_double]
#extern void cut_tracer(Tracer_t *Sample, double rmin, double rmax);
lib.count_tracer_radial.restype=None
lib.count_tracer_radial.argtypes=[Tracer_p, ctypes.c_int]
#extern void count_tracer_radial(Tracer_t *Sample, int nbin);
lib.init_tracer.restype=None
lib.init_tracer.argtypes=[Tracer_p]
#extern void init_tracer(Tracer_t *Sample);
lib.make_sample.restype=None
lib.make_sample.argtypes=[ctypes.c_int, Tracer_p, Tracer_p]
#extern void make_sample(int sample_id, Tracer_t *Sample, Tracer_t *FullSample);

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
	self._as_parameter_=ctypes.byref(self) #this makes it possible to directly pass this as pointer arg
	self._pointer=ctypes.byref(self) #or ctypes.cast(ctypes.pointer(self), Tracer_p)? or ctypes.pointer(self)
	if halo!=None or newoptions!={}:
	  options=get_config(halo)
	  options.update(newoptions)
	  print options
	  for k,v in options.iteritems():
		os.environ[k]=str(v)
	  self.auto_load()
	  
  def __del__(self):
	print "deleting Tracer ", id(self)
	lib.free_tracer(self._pointer)

  def __update_array(self):
	'''this should be called whenever the pointer self.P has changed
	to point the numpy array to the new memory'''
	Parr=(Particle_t*self.nP).from_address(ctypes.addressof(self.P.contents)) #struct array
	self.data=np.frombuffer(Parr, np.dtype(Parr))[0] #the numpy array
  
  def print_data(self):
	print self.nP, '%0x'%ctypes.addressof(self.P.contents)
	print self.P[0].x[0], self.P[0].r
	print self.P[1].x[0], self.P[1].r
	print self.data[1]['x'][0], self.data[1]['r']
	print '-----------'
	lib.print_tracer(self._pointer)
	print '============='
	
  def freeze_energy(self, pars):
	lib.freeze_energy(lib.ParType(pars))

  def likelihood(self, pars, estimator):
	return lib.likelihood(lib.ParType(*pars), estimator, self._pointer)
  
  def freeze_and_like(self, pars, estimator):
	return lib.freeze_and_like(lib.ParType(*pars), estimator, self._pointer)

  def predict_radial_count(self, nbin=100):
	n=np.empty(nbin,dtype='f8')
	lib.predict_radial_count(n.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), nbin)
	#lib.predict_radial_count(n.ctypes.data, nbin) //also work?
	return n
  
  def load(self, infile):
	lib.load_tracer(ctypes.c_char_p(infile), self._pointer)
	self.__update_array()
	
  def copy(self, offset=0, size=0):
	newsample=Tracer()
	lib.copy_tracer(offset, size, newsample._pointer, self._pointer)
	newsample.__update_array()
	return newsample
  
  def radial_cut(self, rmin, rmax):
	lib.cut_tracer(self._pointer, rmin, rmax)
	self.__update_array()
	
  def shuffle(self, seed):
	lib.shuffle_tracer(ctypes.c_ulong(seed), self._pointer)
	
  def resample(self, seed):
	newsample=Tracer()
	lib.resample_tracer(seed, newsample._pointer, self._pointer)
	newsample.__update_array()
	return newsample

  def squeeze(self):
	lib.squeeze_tracer(self._pointer)
	self.__update_array()
	
  def radial_count(self, nbin=None):
	if nbin is None:
	  nbin=lib.NumRadialCountBin
	lib.count_tracer_radial(self._pointer, nbin)
	  
  def sort_flag(self):
	lib.sort_tracer_flag(self._pointer)

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
	lib.make_sample(sampleid, newsample._pointer, self._pointer)
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
	 return 

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
	
if os.uname()[1]=='Medivh':
    rootdir='/work/Projects/DynDistr/'
else:
    rootdir='/gpfs/data/jvbq85/DynDistr/'    

if __name__=="__main__":
  FullSample=Tracer('AqA4',DynSIZE=100)
  Sample=FullSample.sample()
  FullSample.print_data()
  Sample.print_data()
  Sample.data[1]['r']=2
  Sample.rmin=10
  Sample.print_data()
  del FullSample,Sample
  
  FullSample2=Tracer('Mock')
  Sample2=FullSample2.sample()
  FullSample.print_data()
  Sample.print_data()
  Sample.data[1]['r']=2
  Sample.rmin=10
  Sample.print_data()
  del FullSample,Sample