#scan mock ensembles with likelihoods
import matplotlib
matplotlib.use('Agg')
from dynio import *
from myutils import *
import os,sys
#from scipy.stats import chi2

#os.environ['OMP_NUM_THREADS']='24'

estimator=int(sys.argv[1])   #4 or 8(L conditioned) or 10 or 16 or 0 (wenting_like_conditional)

epool=[0,-4,4,8,10,16]
if not estimator in epool:
  print "Err: estimator must be one of",epool,", not %d"%estimator
  raise estimator

if estimator==-4:
  nbinL=10
else:
  nbinL=100
nbinE=1
npart=1000
#init_par=[2,2]
#for estimator==4:
nbin_r=30
FlagRBinLog=1

outdir=lib.rootdir+'/data/MaxLikeEnsembleIniRand'
if not os.path.exists(outdir):
  os.makedirs(outdir)
if estimator==-4:
  outfile=outdir+'/fit4L.dat'
else:
  outfile=outdir+'/fit%d'%estimator+'.dat'
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
	  if estimator==0:
		like=lambda x: -2*Sample.wenting_like_conditional(x)
		result=fmin_gsl(like, init_par, xtol=0.001, full_output=True)
		L1=like([1,1])
	  elif estimator==4:
		Sample.radial_count(nbin_r,FlagRBinLog)
		result=Sample.gfmin_like(estimator, init_par)
		L1=-2*Sample.freeze_and_like([1,1], estimator)
	  elif estimator==-4:
		lib.NumRadialCountBin=nbin_r
		result=Sample.gfmin_jointLE(-estimator, nbinL, nbinE, init_par)
		L1=Sample.jointLE_FChi2([1,1], -estimator, nbinL, nbinE)
	  else:
		result=Sample.gfmin_jointLE(estimator, nbinL, nbinE, init_par)
		L1=Sample.jointLE_FChi2([1,1], estimator, nbinL, nbinE)
	  L=result[1]
	  with open(outfile,'a') as vlog:
		vlog.write("%f\t%f\t"%(result[0][0],result[0][1]))
		vlog.write("%f\t%f\t"%(L,L1))
		vlog.write("%d\n"%(result[-1])) #this differs from previous version

lib.close()