#to do statistical test rather than max-like
#differ in errordef and dof, ts, sigma etc.
import sys,os
from dynio import *
from numpy import random
os.environ['OMP_NUM_THREADS']='24'

mid=init(-1)
freepars='mc'
outfile='/gpfs/data/jvbq85/DynDistr/data/mockfitFmin_mc/fit%d_'%mid+freepars+'.dat'
import shutil
try:
  shutil.move(outfile,outfile+'.bk')
except IOError as e:  
  print "IO Warning:{0} {1}".format(e.strerror, outfile)

from scipy.optimize import *
like= lambda x:-likefunc2(PARTYPE(x[0],x[1]))

for sampleid in range(750):
  print sampleid
  select_particles(sampleid)
  result=fmin(like, 10**(random.rand(2)-0.5), xtol=0.001, ftol=1e-4, maxiter=1000, maxfun=5000, full_output=True)
  L=-result[1]
  L1=-neglike2(1,1)
  with open(outfile,'a') as vlog:
    vlog.write("%f\t%f\t"%(result[0][0],result[0][1]))
    vlog.write("%f\t%f\t"%(L,L1))
    vlog.write("%d\n"%(result[-1]==0))