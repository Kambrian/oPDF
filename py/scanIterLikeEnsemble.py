#scan mock ensembles with likelihoods
import matplotlib
matplotlib.use('Agg')
from dynio import *
from myutils import *
import os,sys
#from scipy.stats import chi2

#os.environ['OMP_NUM_THREADS']='24'

try: 
  estimator=int(sys.argv[1]) #use 10 or 16
  if estimator not in [8,10,16]:
	raise
  proxy=sys.argv[2]
  if proxy not in ['EL','LE','E']:
	raise
except:
  print "Incorrect usage.\n Example: %s 10 EL (or LE or E)\n Now exit."%sys.argv[0]
  raise

npart=1000
#init_par=[2,2]
nbin=10
if proxy=='E':
  nbinE=nbin*nbin
  nbinL=1  
  nbins=[nbinE]
elif proxy=='LE':
  nbinE=nbin
  nbinL=nbin
  nbins=[nbinL,nbinE]
elif proxy=='EL':
  nbinL=nbin
  nbinE=nbin
  nbins=[nbinE,nbinL]

outdir=lib.rootdir+'/data/IterLikeEnsembleIniRand'
if not os.path.exists(outdir):
  os.makedirs(outdir)
outfile=outdir+'/fit%d'%estimator+proxy+'.dat'
import shutil
try:
  shutil.move(outfile,outfile+'.bk')
except IOError as e:  
  print "IO Warning:{0} {1}".format(e.strerror, outfile)

lib.open()
with Tracer('Mock') as FullSample:
  for sampleid in xrange(750):
	init_par=10**(np.random.rand(2)-0.5)
	with FullSample.copy(sampleid*npart,npart) as Sample:
	  result=Sample.gfmin_FixBinIter(estimator, proxy, nbins, init_par, maxiter=10)
	  L1=Sample.nested_views_FChi2([1,1], estimator)
	  L=result[1]
	  with open(outfile,'a') as vlog:
		vlog.write("%f\t%f\t"%(result[0][0],result[0][1]))
		vlog.write("%f\t%f\t"%(L,L1))
		vlog.write("%d\n"%(result[-1])) #this differs from previous version

lib.close()