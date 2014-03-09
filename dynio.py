from ctypes import *
from math import *

#load the library
lib=CDLL("./libdyn.so")
MaxNPar=10
PARTYPE=c_double*MaxNPar
prototype=CFUNCTYPE(c_double, PARTYPE)
likefunc=prototype(("likelihood",lib))#,paramflags)
init=lib.init
#GAMA_init.argtypes=[c_char_p]

neglike=lambda m,c: -likefunc(PARTYPE(m,c))
like=lambda m,c: likefunc(PARTYPE(m,c))

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