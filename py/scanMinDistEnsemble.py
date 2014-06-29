#scan mock ensembles with min-dist
import matplotlib
matplotlib.use('Agg')
from dynio import *
from myutils import *
import os,sys
#from scipy.stats import chi2

#os.environ['OMP_NUM_THREADS']='24'

estimator=int(sys.argv[1])   

elist=[8,10]
if not estimator in elist:
  print "Err: estimator must be one of",elist,", not %d"%estimator
  raise estimator

nbinL=100
nbinE=1
npart=1000
init_par=[2,2]

outdir=lib.rootdir+'/data/MinDistEnsemble'
if not os.path.exists(outdir):
  os.makedirs(outdir)
outfile=outdir+'/fit%d'%estimator+'.dat'
import shutil
try:
  shutil.move(outfile,outfile+'.bk')
except IOError as e:  
  print "IO Warning:{0} {1}".format(e.strerror, outfile)

lib.open()
with Tracer('Mock') as FullSample:
  for sampleid in range(750):
	with FullSample.copy(sampleid*npart,npart) as Sample:
	  result=Sample.fmin_dist(estimator, init_par)
	  L1=-Sample.freeze_and_like([1,1], estimator)
	  L=result[1]
	  with open(outfile,'a') as vlog:
		vlog.write("%f\t%f\t"%(result[0][0],result[0][1]))
		vlog.write("%f\t%f\t"%(L,L1))
		vlog.write("%d\n"%(result[-1])) #this differs from previous version

lib.close()