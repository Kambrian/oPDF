#to do statistical test rather than max-like
#differ in errordef and dof, ts, sigma etc.
import sys,os
from dynio import *
os.environ['OMP_NUM_THREADS']='24'
import numpy as np
from scipy.optimize import *

xlm=[-0,5,0.5]
nbin=100
sampleid=0

mid=init(-1)
outdir='/gpfs/data/jvbq85/DynDistr/data/scan%d'%sampleid
if not os.path.exists(outdir):
  os.makedirs(outdir)
outfile=outdir+'/model%d'%mid+'.dat'

select_particles(sampleid)

for x in np.logspace(xlm[0],xlm[1],nbin):
  for y in np.logspace(xlm[0],xlm[1],nbin):
    with open(outfile,'a') as vlog:
      vlog.write("%f\t%f\t%f\n"%(x,y,like2(x,y)))