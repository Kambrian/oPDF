# todo: stream-cleaned version
# star 500kpc
from dynio import *
from myutils import *
import h5py,os,sys
from scipy.stats import chi2
plt.ion()

def plot_subhalos(sample, proxy, logscale, nPmin=1000):
  '''plot subhaloes as circles in phase space'''
  for i in range(1,Sub.nP):
	mP=Sub.mP*Sub.P[i].w
	if mP>nPmin:
	  if logscale and Sub.data[proxy][i]>0:
		if proxy=='L2':
		  plt.plot(np.log10(Sub.data[proxy][i])/2, Sub.data['theta'][i], 'o', mec='r', mfc='none', markersize=(mP**(1./3))/3)
		else:  
		  plt.plot(np.log10(Sub.data[proxy][i]), Sub.data['theta'][i], 'o', mec='r', mfc='none', markersize=(mP**(1./3))/3)
	  else:
		plt.plot((Sub.data[proxy][i]), Sub.data['theta'][i], 'o', mec='r', mfc='none', markersize=(mP**(1./3))/3)
	
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
  
def phase_imshow(sample, proxy, bins, method='hist', logscale=True):
  X,Y,Z,extent,data=sample.phase_density(proxy, bins, method, logscale)
  #Z[Z==0]=np.nan;
  extent=list(extent)
  if proxy=='L2':
	extent[0]/=2
	extent[1]/=2
  plt.imshow(Z,extent=extent)
  plt.axis('tight')
  if sample.FlagUseWeight:
	if logscale:
	  w=sample.data['w'][sample.data[proxy]>0]
	else:
	  w=sample.data['w']
  else:
	w=None
  sk=skeleton(data[0],data[1],nbin=50,weights=w)
  if proxy=='L2':
	plt.plot(sk['x']['median']/2, sk['y']['mean'], 'w-')
  else:
	plt.plot(sk['x']['median'], sk['y']['mean'], 'w-')
  #yerr=1/np.sqrt(12.*sk['x']['hist']) #scaled by bin count
  #plt.plot(sk['x']['median'], sk['y']['mean']+yerr, 'w--')
  #plt.plot(sk['x']['median'], sk['y']['mean']-yerr, 'w--')
  plt.plot(extent[:2], [0.5, 0.5], 'k-')

halo=sys.argv[1] #'A4'
useweight=0 #int(sys.argv[2])
cleanmethod='sub' #sys.argv[3]
if cleanmethod not in ['sub','branch']:
  print 'cleanmethod must be sub or branch. now exit.'
  exit()
npart=1e6
rmin=10
rmax=300

lib.open()

with Tracer(halo,DynRMIN=rmin,DynRMAX=rmax,DynDataFile='old/'+halo+'star.hdf5') as FullSample:
  #if not useweight:
	#tmp=FullSample.select((FullSample.data['haloid']>0))
	#FullSample.clean()
	#FullSample=tmp
  if npart>FullSample.nP:
	npart=0  
  Sample=FullSample.copy(0,npart)
Sample.FlagUseWeight=useweight
Sub=Tracer(halo, DynRMIN=rmin, DynRMAX=rmax, DynDataFile=halo+'sub.hdf5')  
lib.init_potential_spline()

Sample.freeze_energy();
Sample.like_init()
Sub.freeze_energy();
Sub.like_init()

Sample.AD=-Sample.like_eval(estimator=8)
Sample.Mean=Sample.like_eval(estimator=10)

if cleanmethod=='sub':
  SampleClean=Sample.select(Sample.data['subid']<=0)
  SampleSub=Sample.select(Sample.data['subid']>0)
else:
  SampleClean=Sample.select(Sample.data['haloid']==0) #branch clean
  SampleSub=Sample.select(Sample.data['haloid']>0)  
SampleClean.AD=-SampleClean.like_eval(estimator=8)
SampleClean.Mean=SampleClean.like_eval(estimator=10)
SampleSub.AD=-SampleSub.like_eval(estimator=8)
SampleSub.Mean=SampleSub.like_eval(estimator=10)
 
nbin=[80, 50]
method='hist'
proxy='E'
#lvls=10
#percents=None
#percents=np.linspace(0.05,1,10) #initial levels determined by percents
logscale=True
#extent=[4,5.5]
extent={'A2':[4,5.3],'A4':[4.2,5.2],'B2':[4.2,5.1],'B4':[4,5],'C2':[4,5.5],'D2':[4,5.5],'E2':[4,5.5]}[halo]
extent.extend([0,1])
bins=[np.linspace(extent[0], extent[1], nbin[0]+1), np.linspace(extent[2], extent[3], nbin[1]+1)]
plt.figure(figsize=(18,7))
for i,sample in enumerate([Sample, SampleClean, SampleSub]):
  plt.subplot(1,3,i+1)
  phase_imshow(sample, proxy, bins, method, logscale)
  #plt.title(r'$\bar{\Theta}_{'+halo+'}=$%.1f'%sample.Mean)
  plt.title(r'$D_{'+halo+'}=$%.1f'%sample.AD)
  #lvls=phase_contour(sample, proxy, bins, method, logscale, lvls, percents)
  #percents=None
plot_subhalos(Sub, proxy, logscale)
plt.subplot(132)
plt.xlabel(r'$\log(E[{\rm km/s}]^2)$')
plt.subplot(131)
plt.ylabel(r'$\theta$')
#plt.subplot(133)
#plt.colorbar()
if cleanmethod=='sub':
  plt.savefig(lib.rootdir+'/plots/paper/extra/PhaseDistr'+proxy+'.'+halo+'.R%d.Weight%d.star.Dup.pdf'%(rmax,useweight))#, rasterize=True, dpi=300)
else:
  plt.savefig(lib.rootdir+'/plots/paper/extra/PhaseDistr'+proxy+'.'+halo+'.R%d.Weight%d.BranchClean.star.Dup.pdf'%(rmax,useweight))#, rasterize=True, dpi=300)

proxy='L2'
#extent=[3,10]
extent={'A2':[3.6,10],'A4':[6.5,10],'B2':[3,9.5],'B4':[6,9.5],'C2':[6,10],'D2':[6,10],'E2':[6,10]}[halo]
extent.extend([0,1])
bins=[np.linspace(extent[0], extent[1], nbin[0]+1), np.linspace(extent[2], extent[3], nbin[1]+1)]
plt.figure(figsize=(18,7))
for i,sample in enumerate([Sample, SampleClean, SampleSub]):
  plt.subplot(1,3,i+1)
  phase_imshow(sample, proxy, bins, method, logscale)
  #plt.title(r'$\bar{\Theta}_{'+halo+'}=$%.1f'%sample.Mean)
  #plt.title(r'$D_{'+halo+'}=$%.1f'%sample.AD)
plot_subhalos(Sub, proxy, logscale)
plt.subplot(132)
plt.xlabel(r'$\log(L[{\rm kpc\cdot km/s}])$')
plt.subplot(131)
plt.ylabel(r'$\theta$')
if cleanmethod=='sub':
  plt.savefig(lib.rootdir+'/plots/paper/extra/PhaseDistr'+proxy+'.'+halo+'.R%d.Weight%d.star.Dup.pdf'%(rmax,useweight))#, rasterize=True, dpi=300)
else:
  plt.savefig(lib.rootdir+'/plots/paper/extra/PhaseDistr'+proxy+'.'+halo+'.R%d.Weight%d.BranchClean.star.Dup.pdf'%(rmax,useweight))#, rasterize=True, dpi=300)


lib.free_potential_spline()
FullSample.clean()
Sample.clean()
SampleClean.clean()
SampleSub.clean()
Sub.clean()
lib.close()