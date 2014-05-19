#to do statistical test rather than max-like
#differ in errordef and dof, ts, sigma etc.
import sys,os
from dynio import *
from numpy import random
os.environ['OMP_NUM_THREADS']='24'

fullpars='mc'
if len(sys.argv)>1:
  freepars=sys.argv[1]
else:
  freepars=fullpars

mid=init(-1)
outfile='/gpfs/data/jvbq85/DynDistr/data/mockfitHighPrec_mc/fit%d_'%mid+freepars+'.dat'
import shutil
try:
  shutil.move(outfile,outfile+'.bk')
except IOError as e:  
  print "IO Warning:{0} {1}".format(e.strerror, outfile)

from iminuit import Minuit
from iminuit.ConsoleFrontend import ConsoleFrontend

fixnames=fullpars.translate(None,freepars)
fixdict={} #gen a dict for fixed pars
for x in fixnames: 
  fixdict['fix_'+x]=True
x=10.**(random.rand(2)-0.5)  
defaults={'m':x[0],'c':x[1]}
defaults.update(fixdict)
for sampleid in range(750):
  print sampleid
  select_particles(sampleid)
  m=Minuit(neglike2,print_level=0,errordef=1.0, limit_m=[0.1,10],limit_c=[0.1,10], pedantic=False, frontend=ConsoleFrontend(),**defaults)
  m.set_strategy(2);
  m.tol=0.1 #default convergence edm<1e-4*tol*errordef, but we do not need that high accuracy
  result=m.migrad()
  L=-m.fval
  L1=-neglike2(1,1)
  with open(outfile,'a') as vlog:
    for par in ['m','c']:
      if par in m.values:
	vlog.write("%f\t"%(m.values[par]))
    vlog.write("%f\t%f\t"%(L,L1))
    try:
      vlog.write("%f\t"%m.matrix(True)[0][1]) #newly added, correlation coef
    except:
      vlog.write("nan\t")
    vlog.write("%d\n"%result[0]['is_valid'])