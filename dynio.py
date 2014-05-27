from ctypes import *
from math import *
import numpy as np
import os,ConfigParser

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
squeeze_data=lib.squeeze_data
squeeze_data.restype=c_int
free_data=lib.free_data

#wenting's lib
prototype=CFUNCTYPE(c_double, c_double, c_double, c_double)
dataprob=prototype(("dataprob",lib))
PARTYPEXV=c_double*6
prototype=CFUNCTYPE(c_double, PARTYPEXV)
dataprob6d=prototype(("dataprob6d",lib))


#neglike=lambda m,c,estimator: -likefunc(PARTYPE(m,c),estimator)
#like=lambda m,c,estimator: likefunc(PARTYPE(m,c),estimator)
freeze_energy=lambda m,c: freezefunc(PARTYPE(m,c))
#like2=lambda m,c,estimator: likefunc2(PARTYPE(m,c),estimator) #freeze_and_like()
#neglike2=lambda m,c,estimator: -likefunc2(PARTYPE(m,c),estimator)

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
  import os
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