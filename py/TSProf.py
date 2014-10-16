import matplotlib
matplotlib.use('Agg')
from dynio import *
from myutils import *
import h5py,os,sys
from scipy.stats import chi2
#plt.ion()

#useweight=int(sys.argv[1])
useweight=0
estimator=int(sys.argv[1])
rmin=10
rmax=300
pars=[1,1]
npart=0 
nbin=50

outfile=lib.rootdir+'/plots/paper/star/TSprof'+lib.NameList[estimator]+'.hdf5'

lib.open()
f=h5py.File(outfile,'w')
for halo in 'ABCDE':
  with Tracer(halo+'2',DynRMIN=rmin,DynRMAX=rmax,DynDataFile=halo+'2star.hdf5') as F:
	with F.select(F.data['subid']<=0) as C:
	  Sample=C.copy(0,int(npart))
  Sample.FlagUseWeight=useweight
  lib.init_potential_spline()
  Sample.freeze_energy(pars)
  for proxy in ['r','E','L2']:
	print proxy
	bins=Sample.gen_bin(proxy, nbin, equalcount=True)[0]
	x,tscum,count=Sample.TSprofCum([1,1], estimator, proxy, bins)
	tsdiff,x=Sample.TSprof(pars, estimator, viewtypes=proxy, nbins=nbin)
	tsdiff=tsdiff[:-1]
	if estimator==8:
	  tsdiff=AD2Sig(-tsdiff)
	  tscum=AD2Sig(-tscum)
	f.create_dataset('/'+halo+'/'+proxy+'/x', data=x)
	f.create_dataset('/'+halo+'/'+proxy+'/tscum', data=tscum)
	f.create_dataset('/'+halo+'/'+proxy+'/tsdiff', data=tsdiff)
	plt.plot(tsdiff)
  Sample.clean()
  lib.free_potential_spline()
  
f.close()  
lib.close()

f=h5py.File(outfile,'r')

fig=plt.figure(figsize=(8,12))
for i,proxy in enumerate(['r','E','L2']):
  plt.subplot(3,1,i+1)
  for halo in 'ABCDE':  
	ts=f['/'+halo+'/'+proxy+'/tscum'][...]
	n=len(ts)
	x=(np.arange(n)+1.)/n
	if proxy=='E':
	  ts=ts[-1::-1]
	plt.plot(x*100, ts/np.sqrt(x))
  plt.plot(plt.xlim(),[0,0],'k--')
  if proxy=='E':
	plt.xlabel('-E Percentile')
  else:
	plt.xlabel(proxy[0]+' Percentile')
plt.subplot(312)
plt.ylabel(r'Cumulative $D[\sigma]$')
fig.subplots_adjust(hspace=0.25, bottom=0.1)
plt.subplot(311)
plt.legend(list('ABCDE'), ncol=3, loc='upper left')
plt.savefig(lib.rootdir+'/plots/paper/star/TSprofCum'+lib.NameList[estimator]+'.scaled.eps') #rasterize=True, dpi=300

fig=plt.figure(figsize=(8,12))
for i,proxy in enumerate(['r','E','L2']):
  plt.subplot(3,1,i+1)
  for halo in 'ABCDE':  
	ts=f['/'+halo+'/'+proxy+'/tsdiff'][...]
	n=len(ts)
	plt.plot((np.arange(n)+1.)/n*100, ts)
  plt.plot(plt.xlim(),[0,0],'k--')
  plt.xlabel(proxy[0]+' Percentile')
plt.subplot(312)
plt.ylabel(r'$D[\sigma]$')  
fig.subplots_adjust(hspace=0.25, bottom=0.1)
plt.subplot(311)
plt.legend(list('ABCDE'), ncol=3, loc='upper left')
plt.savefig(lib.rootdir+'/plots/paper/star/TSprof'+lib.NameList[estimator]+'.eps') #rasterize=True, dpi=300

f.close()