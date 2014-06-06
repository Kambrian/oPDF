from ctypes import *
from math import *
import numpy as np
import os,ConfigParser,h5py

NameList={0:'f(E,L)',4:'RBin',8:'AD',9:'Resultant',10:'Mean',11:'KS',12:'Kuiper',13:'CosMean'}

#load the library
lib=CDLL("./libdyn.so")
MaxNPar=10
PARTYPE=c_double*MaxNPar
prototype=CFUNCTYPE(c_double, PARTYPE, c_int)
likefunc=prototype(("likelihood",lib))#,paramflags)
likefunc2=prototype(("freeze_and_like",lib))
prototype=CFUNCTYPE(None, PARTYPE)
freezefunc=prototype(("freeze_energy",lib))
init=lib.init
#init.argtypes=[c_int]
#init.restype=c_int
select_particles=lib.select_particles
select_particles.argtypes=[c_int]
lib.squeeze_data.restype=c_int #avoid using this directly. use P.squeeze_data() instead.
#squeeze_data=lib.squeeze_data
#squeeze_data.restype=c_int
free_data=lib.free_data
lib.halo_pot.argtypes=[c_double]
lib.halo_pot.restype=c_double
lib.like_init.argtypes=[PARTYPE,c_int]
like_init=lambda m,c,estimator: lib.like_init(PARTYPE(m,c), estimator)
lib.predict_radial_count.argtypes=[POINTER(c_double), c_int]
def predict_radial_count(nbin=100):
	n=np.empty(nbin,dtype='f8')
	lib.predict_radial_count(n.ctypes.data_as(POINTER(c_double)), nbin)
	return n

#wenting's lib
prototype=CFUNCTYPE(c_double, c_double, c_double, c_double)
dataprob=prototype(("dataprob",lib))
PARTYPEXV=c_double*6
prototype=CFUNCTYPE(c_double, PARTYPEXV)
dataprob6d=prototype(("dataprob6d",lib))


elike=lambda m,c,estimator: likefunc(PARTYPE(m,c),estimator)
freeze_energy=lambda m,c: freezefunc(PARTYPE(m,c))
elike2=lambda m,c,estimator: likefunc2(PARTYPE(m,c),estimator) #freeze_and_like()

def gen_par(p):
  "convert parameter array p into PARTYPE() array"
  par=PARTYPE()
  for i,x in enumerate(p):
    par[i]=x
  return par

def pack_to_dict(p):
  "pack parameter values p[n] into a dict(a=,b=,c=,...) to be used as default pars for minuit"
  d=dict()
  for i,pp in enumerate(p):
    d[parname_repo[i]]=pp
  return d

class NFWHalo(Structure):
  """all properties are physical"""
  _fields_=[('z', c_double),
	    ('M', c_double),
	    ('c', c_double),
	    ('Rv', c_double), 
	    ('Rs', c_double),
	    ('Rhos', c_double),
	    ('Pots', c_double),#-4*pi*G*rhos*rs^2, the potential at r=0
	    ('virtype', c_int)
	    ]
Halo=NFWHalo.in_dll(lib, 'Halo')
comoving_virial_radius=lib.comoving_virial_radius
comoving_virial_radius.argtypes=[c_double, c_double, c_int]
comoving_virial_radius.restype=c_double

class Particle(Structure):
  _fields_=[('flag', c_int),
	    ('r', c_double),
	    ('K', c_double),
	    ('L2', c_double),
	    ('x', c_double*3),
	    ('v', c_double*3),
	    ('E', c_double),
	    ('T', c_double),
	    ('vr', c_double),
	    ('theta', c_double),
	    ('rlim', c_double*2)
	    ]  

class ParticleData:
  def __init__(self):
    self.SIZE=c_int.in_dll(lib,'SUBSAMPLE_SIZE')
    self.R_MIN=c_double.in_dll(lib,'R_MIN')
    self.R_MAX=c_double.in_dll(lib,'R_MAX')
    self.nP=c_int.in_dll(lib,'nP')
    #self.StructP=(Particle*self.nP.value).in_dll(lib,'P')
    self.P2P=POINTER(Particle).in_dll(lib,'P')
    self.StructP=(Particle*self.nP.value).from_address(addressof(self.P2P.contents))
    self.P=np.frombuffer(self.StructP,np.dtype(self.StructP))[0]
    #P=from_buffer(cast(P2P, POINTER(Particle*nP.value)).contents)

  def squeeze_data(self):
	  lib.squeeze_data()
	  self.__init__()
	  
  def print_data(self):
    print self.nP.value, '%0x'%addressof(self.P2P.contents)
    print self.P2P[0].x[0], self.P2P[0].r
    print self.P2P[1].x[0], self.P2P[1].r
    print '-----------'
    lib.print_data()
    print '============='

  def gen_bin(self,bintype,nbin=30,logscale=True, equalcount=False):
	proxy=np.copy(self.P[bintype])
	n=nbin+1
	if bintype=='E':
		proxy=-proxy
	if equalcount:
		x=np.array(np.percentile(proxy,list(np.linspace(0,100,n))))
	else:  
		if logscale:  
			x=np.logspace(np.log10(proxy[proxy>0].min()),np.log10(proxy.max()),n)
		else:
			x=np.linspace(proxy.min(),proxy.max(),n)
	return x,proxy
  
  def assign_TSmap(self, TS, x):
    '''select particles with TS in the range [lim[0],lim[1]], 
    where TS is a map with pixel boundaries in x, 
    and x is an ordereddict
    return: update P['flag'] in place.'''
    p=((TS>lim[0])*(TS<lim[1])).nonzero()
    b=x.keys()
    xx=x.values()
    f=np.zeros(self.nP.value)
    for pix in zip(p[0],p[1]):
		f=f+(self.P[b[0]]>xx[0][pix[0]])*(self.P[b[0]]<xx[0][pix[0]+1])*(self.P[b[1]]>xx[1][pix[1]])*(self.P[b[1]]<xx[1][pix[1]+1])
    f=f>0
    print "Found %d particles"%np.sum(f)  
    self.P['flag'][:]=f
    
  def filter_TSmap(self, TS, x, lim):
		'''select particles with TS in the range [lim[0],lim[1]], 
		where TS is a map with pixel boundaries in x, 
		and x is an ordereddict
		return: update P['flag'] in place.'''
		self.P['flag']=0
		p=((TS>lim[0])*(TS<lim[1])).nonzero()
		b=x.keys()
		xx=x.values()
		for pix in zip(p[0],p[1]):
			ff=np.ones(self.nP.value)
			for i in [0,1]:
				if b[i]=='E':
					ff=ff*(self.P[b[i]]<-xx[i][pix[i]])*(self.P[b[i]]>-xx[i][pix[i]+1])
				else:
					ff=ff*(self.P[b[i]]>xx[i][pix[i]])*(self.P[b[i]]<xx[i][pix[i]+1])
			self.P['flag'][ff>0]=1
		print "Found %d pixels, %d particles"%(p[0].shape[0],np.sum(self.P['flag']))

class SubData:
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
	self.E=np.array([k+lib.halo_pot(r) for k,r in zip(self.K,self.r)])
def get_config(halo):
  c=ConfigParser.ConfigParser()
  c.optionxform=str
  c.read('DataFiles.cfg')
  for x in c.options(halo):
    os.environ[x]=c.get(halo,x)

if os.uname()[1]=='Medivh':
    rootdir='/work/Projects/DynDistr/'
else:
    rootdir='/gpfs/data/jvbq85/DynDistr/'    
    
if __name__=="__main__":
  get_config('AqA4')
  init()
  select_particles(0)
  P=ParticleData()
  P.print_data()
  P.P[1]['r']=2
  P.R_MIN.value=10
  P.print_data()
  free_data()
  
  get_config('Mock')
  init()
  select_particles(0)
  P=ParticleData()
  P.print_data()
  free_data()