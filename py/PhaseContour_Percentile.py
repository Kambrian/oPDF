# todo: stream-cleaned version
# star 500kpc
from dynio import *
from myutils import *
import h5py,os,sys
from scipy.stats import chi2
plt.ion()

def plot_subhalos(sample, sub, proxy, nPmin=1000):
  '''plot subhaloes as circles in phase space'''
  with sub.select(sub.mP*sub.data['w']>nPmin) as S:
	p=np.sort(sample.data[proxy])
	x=np.searchsorted(p, S.data[proxy])/float(sample.nP)
	plt.scatter(x, S.data['theta'], marker='o', edgecolor='r', facecolor='none', s=(S.mP*S.data['w'])/300)
	
def phase_contour(sample, proxy, bins, method, logscale, levels=5, percents=None):
  '''percents will overide levels if present'''
  X,Y,Z,_,_=sample.phase_density(proxy, bins, method, logscale)
  #Z[Z==0]=np.nan;
  if percents is None:
	cs=plt.contour(X,Y,Z,levels)
	levels=cs.levels
  else:
	h,h0,levels=percentile_contour(X,Y,Z, percents=percents)
  #plt.axis('tight')
  return levels
  
def phase_imshow(sample, proxy, bins, method='hist'):
  X,Y,Z,extent,data=sample.phase_density(proxy, bins, method, False)
  Z[Z==0]=np.nan;
  plt.imshow(Z,extent=[0,1,0,1])
  plt.axis('tight')
  sk=skeleton(data[0],data[1],nbin=bins[0])
  nbin=len(bins[0])-1
  plt.plot((np.arange(nbin)+0.5)/nbin, sk['y']['mean'], 'w-')
  plt.plot([0,1], [0.5, 0.5], 'k-')

halo=sys.argv[1] #'A4'
npart=int(1e6)
rmin=1
rmax=300

#lib.open()
with Tracer(halo,DynRMIN=rmin,DynRMAX=rmax,DynDataFile=halo+'DM.hdf5') as FullSample:
  Sample=FullSample.copy(0,npart)
Sub=Tracer(halo, DynRMIN=rmin, DynRMAX=rmax, DynDataFile=halo+'sub.hdf5')  
lib.init_potential_spline()

Sample.freeze_energy();
Sample.like_init()
Sub.freeze_energy();
Sub.like_init()

Sample.AD=-Sample.like_eval(estimator=8)
Sample.Mean=Sample.like_eval(estimator=10)

SampleClean=Sample.select(Sample.data['subid']<=0)
SampleClean.AD=-SampleClean.like_eval(estimator=8)
SampleClean.Mean=SampleClean.like_eval(estimator=10)

SampleSub=Sample.select(Sample.data['subid']>0)
SampleSub.AD=-SampleSub.like_eval(estimator=8)
SampleSub.Mean=SampleSub.like_eval(estimator=10)
 
nbin=[50, 50]
method='hist'
proxy='E'
logscale=False
bins=[Sample.gen_bin(proxy, nbin[0], equalcount=True)[0],np.linspace(0, 1, nbin[1]+1)]

plt.figure(figsize=(18,7))
for i,sample in enumerate([Sample, SampleClean, SampleSub]):
  plt.subplot(1,3,i+1)
  phase_imshow(sample, proxy, bins, method)
  #plt.title(r'$\bar{\Theta}_{'+halo+'}=$%.1f'%sample.Mean)
  plt.title(r'$D_{'+halo+'}=$%.1f'%sample.AD)
  #lvls=phase_contour(sample, proxy, bins, method, logscale, lvls, percents)
  #percents=None
plot_subhalos(Sample, Sub, proxy)
plt.subplot(132)
plt.xlabel(r'$E$ percentile')
plt.subplot(131)
plt.ylabel(r'$\theta$')
#plt.subplot(133)
#plt.colorbar()
plt.savefig(lib.rootdir+'/plots/paper/extra/PhaseDistr'+proxy+'Percent.'+halo+'.R%d.pdf'%rmax)#, rasterize=True, dpi=300)

proxy='L2'
bins=[Sample.gen_bin(proxy, nbin[0], equalcount=True)[0],np.linspace(0, 1, nbin[1]+1)]
plt.figure(figsize=(18,7))
for i,sample in enumerate([Sample, SampleClean, SampleSub]):
  plt.subplot(1,3,i+1)
  phase_imshow(sample, proxy, bins, method)
  #plt.title(r'$\bar{\Theta}_{'+halo+'}=$%.1f'%sample.Mean)
  #plt.title(r'$D_{'+halo+'}=$%.1f'%sample.AD)
plot_subhalos(Sample, Sub, proxy)
plt.subplot(132)
plt.xlabel(r'$L$ percentile')
plt.subplot(131)
plt.ylabel(r'$\theta$')
plt.savefig(lib.rootdir+'/plots/paper/extra/PhaseDistr'+proxy+'Percent.'+halo+'.R%d.pdf'%rmax)#, rasterize=True, dpi=300)

#lib.free_potential_spline()
#lib.close()