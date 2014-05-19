from ctypes import *
from math import *
import numpy as np

#load the library
lib=CDLL("./libdyn.so")
MaxNPar=10
PARTYPE=c_double*MaxNPar
prototype=CFUNCTYPE(c_double, PARTYPE)
likefunc=prototype(("likelihood",lib))#,paramflags)
likefunc2=prototype(("freeze_and_like",lib))
prototype=CFUNCTYPE(None, PARTYPE)
freezefunc=prototype(("freeze_energy",lib))
prototype=CFUNCTYPE(c_double, c_double, c_double, c_double)
dataprob=prototype(("dataprob",lib))
PARTYPEXV=c_double*6
prototype=CFUNCTYPE(c_double, PARTYPEXV)
dataprob6d=prototype(("dataprob6d",lib))
init=lib.init
init.argtypes=[c_int]
init.restype=c_int
select_particles=lib.select_particles
select_particles.argtypes=[c_int]

neglike=lambda m,c: -likefunc(PARTYPE(m,c))
like=lambda m,c: likefunc(PARTYPE(m,c))
freeze_energy=lambda m,c: freezefunc(PARTYPE(m,c))
like2=lambda m,c: likefunc2(PARTYPE(m,c)) #freeze_and_like()
neglike2=lambda m,c: -likefunc2(PARTYPE(m,c))

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

if __name__=="__main__":
  init(-1)
  select_particles(0)
  P=ParticleData()
  P.print_data()
  P.P[1]['r']=2
  P.print_data()