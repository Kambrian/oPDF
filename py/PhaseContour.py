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
		plt.plot(np.log10(Sub.data[proxy][i]), Sub.data['theta'][i], 'o', mec='r', mfc='none', markersize=(mP**(1./3))/5)
	  else:
		plt.plot((Sub.data[proxy][i]), Sub.data['theta'][i], 'o', mec='r', mfc='none', markersize=(mP**(1./3))/5)
	
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
  plt.imshow(Z,extent=extent)
  plt.axis('tight')
  sk=skeleton(data[0],data[1],nbin=50)
  plt.plot(sk['x']['median'], sk['y']['mean'], 'w-')
  yerr=1/np.sqrt(12.*sk['x']['hist']) #scaled by bin count
  plt.plot(sk['x']['median'], sk['y']['mean']+yerr, 'w--')
  plt.plot(sk['x']['median'], sk['y']['mean']-yerr, 'w--')
  plt.plot(extent[:2], [0.5, 0.5], 'k-')

halo=sys.argv[1] #'A4'
npart=int(1e6)
rmin=1.1
rmax=100

lib.open()
#with Tracer(halo,DynRMIN=rmin,DynRMAX=rmax,DynDataFile=halo+'DM.hbt.hdf5') as FullSample:
  #SampleHBT=FullSample.copy(0,npart)
#SubHBT=Tracer(halo, DynRMIN=rmin, DynRMAX=rmax, DynDataFile=halo+'sub.hbt.hdf5')

#lib.init_potential_spline()
#SubHBT.freeze_energy();
#SubHBT.like_init()
#SampleHBT.freeze_energy()
#SampleHBT.like_init()

#SampleHBT.AD=-SampleHBT.like_eval(estimator=8)
#SampleHBT.Mean=SampleHBT.like_eval(estimator=10)

#SampleHBTClean=SampleHBT.select(SampleHBT.data['subid']<=0)
#SampleHBTClean.AD=-SampleHBTClean.like_eval(estimator=8)
#SampleHBTClean.Mean=SampleHBTClean.like_eval(estimator=10)

#SampleHBTSub=SampleHBT.select(SampleHBT.data['subid']>0)
#SampleHBTSub.AD=-SampleHBTSub.like_eval(estimator=8)
#SampleHBTSub.Mean=SampleHBTSub.like_eval(estimator=10)

#SampleHBTStrm=SampleHBT.select(SampleHBT.data['strmid']>0)
#SampleHBTStrmClean=SampleHBT.select(SampleHBT.data['strmid']<=0)
#SampleHBTStrmClean.AD=-SampleHBTStrmClean.like_eval(estimator=8)
#SampleHBTStrm.AD=-SampleHBTStrm.like_eval(estimator=8)
#SampleHBTStrmClean.Mean=SampleHBTStrmClean.like_eval(estimator=10)
#SampleHBTStrm.Mean=SampleHBTStrm.like_eval(estimator=10)
#lib.free_potential_spline()

with Tracer(halo+'N',DynRMIN=rmin,DynRMAX=rmax,DynDataFile=halo+'DM.hdf5') as FullSample:
  Sample=FullSample.copy(0,npart)
Sub=Tracer(halo+'N', DynRMIN=rmin, DynRMAX=rmax, DynDataFile=halo+'sub.hdf5')  
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
 
nbin=[80, 50]
method='hist'
proxy='E'
#lvls=10
#percents=None
#percents=np.linspace(0.05,1,10) #initial levels determined by percents
logscale=True
extent={'A2':[4.2,5.2],'A4':[4.2,5.2],'B2':[4,5],'B4':[4,5],'C2':[4,5.5],'D2':[4,5.5],'E2':[4,5.5]}[halo]
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
plt.xlabel(r'$\log10(E)$')
plt.subplot(131)
plt.ylabel(r'$\theta$')
plt.savefig(lib.rootdir+'/plots/paper/extra/PhaseDistr'+proxy+'.'+halo+'.R%d.pdf'%rmax)#, rasterize=True, dpi=300)

proxy='L2'
extent={'A2':[6.5,10],'A4':[6.5,10],'B2':[6,9.5],'B4':[6,9.5],'C2':[6,10],'D2':[6,10],'E2':[6,10]}[halo]
extent.extend([0,1])
bins=[np.linspace(extent[0], extent[1], nbin[0]+1), np.linspace(extent[2], extent[3], nbin[1]+1)]
plt.figure(figsize=(18,7))
for i,sample in enumerate([Sample, SampleClean, SampleSub]):
  plt.subplot(1,3,i+1)
  phase_imshow(sample, proxy, bins, method, logscale)
  #plt.title(r'$\bar{\Theta}_{'+halo+'}=$%.1f'%sample.Mean)
  plt.title(r'$D_{'+halo+'}=$%.1f'%sample.AD)
plot_subhalos(Sub, proxy, logscale)
plt.subplot(132)
plt.xlabel(r'$\log10(L^2)$')
plt.subplot(131)
plt.ylabel(r'$\theta$')
plt.savefig(lib.rootdir+'/plots/paper/extra/PhaseDistr'+proxy+'.'+halo+'.R%d.pdf'%rmax)#, rasterize=True, dpi=300)

#plt.subplot(211)
#f=Sample.data['E']>0
#X,Y,Z=density_of_points(np.array([np.log10(Sample.data['E'][f]), Sample.data['theta'][f]]), bins=nbin, method=method)
#extent=get_extent(X,Y)
#bins=[np.linspace(extent[0], extent[1], nbin[0]+1), np.linspace(extent[2], extent[3], nbin[1]+1)]
#Z[Z==0]=np.nan; plt.imshow(Z, extent=extent)
#h,h0,lvls=percentile_contour(X,Y,Z, percents=percents, linestyles='dashed', alpha=0.5)
#cs=plt.contour(X,Y,Z, linestyles='dashed', alpha=0.5);lvls=cs.levels
#f=(Sample.data['E']>0)*(Sample.data['subid']>0)
#plt.scatter(np.log10(Sample.data['E'][f]), Sample.data['theta'][f], '.', c=Sample.data['subid'][f])
#plt.contour(*density_of_points(np.array([np.log10(Sample.data['E'][f]), Sample.data['theta'][f]]), bins=bins, method=method))#,
			#levels=lvls, linestyles='solid')
#for i in range(300):
	#plt.plot(np.log10(Sublist.data['E'][i]), Sublist.data['theta'][i], 'o', markerfacecolor='w', markersize=40./(i+1.)**0.6)
#plt.axis(extent)




#lib.free_potential_spline()
#FullSample.clean()
#lib.close()