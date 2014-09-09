# todo: stream-cleaned version
# star 500kpc
from dynio import *
from myutils import *
import h5py,os,sys
from scipy.stats import chi2
plt.ion()

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
  plt.plot(extent[:2], [0., 0.], 'k-')

halo='AqA4'
npart=int(1e6)

lib.open()
FullSample=Tracer(halo+'allN',DynRMIN=1.1,DynRMAX=300)
lib.init_potential_spline()
FullSample.freeze_energy();
Sample=FullSample.copy(0,npart)
Sample.like_init()

SampleClean=Sample.select(Sample.data['subid']<=0)
SampleSub=Sample.select(Sample.data['subid']>0)
SampleStrm=Sample.select(Sample.data['strmid']>0)
SampleStrmClean=Sample.select(Sample.data['strmid']<=0)

Sample.AD=-Sample.like_eval(estimator=8)
SampleClean.AD=-SampleClean.like_eval(estimator=8)
SampleSub.AD=-SampleSub.like_eval(estimator=8)
SampleStrmClean.AD=-SampleStrmClean.like_eval(estimator=8)
SampleStrm.AD=-SampleStrm.like_eval(estimator=8)

Sample.Mean=Sample.like_eval(estimator=10)
SampleClean.Mean=SampleClean.like_eval(estimator=10)
SampleSub.Mean=SampleSub.like_eval(estimator=10)
SampleStrmClean.Mean=SampleStrmClean.like_eval(estimator=10)
SampleStrm.Mean=SampleStrm.like_eval(estimator=10)
 
nbin=[100, 50]
method='hist'
proxy='E'
#lvls=10
#percents=None
#percents=np.linspace(0.05,1,10) #initial levels determined by percents
for logscale,extent in [(False,[0, 2e5, -0.5, 0.5]), (True, [4,5.5,-0.5,0.5])]:
  bins=[np.linspace(extent[0], extent[1], nbin[0]+1), np.linspace(extent[2], extent[3], nbin[1]+1)]
  plt.figure()
  for i,sample in enumerate([Sample, SampleClean, SampleSub, Sample, SampleStrmClean, SampleStrm]):
	plt.subplot(2,3,i+1)
	phase_imshow(sample, proxy, bins, method, logscale)
	#lvls=phase_contour(sample, proxy, bins, method, logscale, lvls, percents)
	#percents=None
  
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