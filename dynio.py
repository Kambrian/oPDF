from ctypes import *
from math import *

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