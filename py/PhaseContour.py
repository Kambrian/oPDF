#using nested_views rather than jointEL
from dynio import *
from myutils import *
import h5py,os,sys
from scipy.stats import chi2
plt.ion()

halo='AqA4'
npart=int(1e5)

lib.open()
FullSampleAll=Tracer(halo+'allN',DynRMIN=1.1,DynRMAX=499)
All=FullSampleAll.copy(0,npart)
FullSampleFoF=Tracer(halo+'N',DynRMIN=1.1,DynRMAX=499)
FoF=FullSampleFoF.copy(0,npart)
#with Tracer(halo+'subN',DynRMIN=1.1,DynRMAX=499) as FullSample:
  #Sub=FullSample.copy(0,npart)
#with Tracer(halo+'strmN',DynRMIN=1.1,DynRMAX=499) as FullSample:
  #Strm=FullSample.copy(0,npart)  
with Tracer(halo+'sublistN',DynRMIN=1.1,DynRMAX=499) as FullSample:
  Sublist=FullSample.copy(0,0)
  
percents=np.linspace(0.05,1,10)  
lib.init_potential_spline()
FullSampleAll.freeze_energy(); FullSampleFoF.freeze_energy()
#FullSampleAll.like_init(); FullSampleFoF.like_init()
for Sample,tag in [(All,'all'), (FoF,'fof'), (Sublist,'sat')]: #(Sub,'main'), (Strm, 'strm')
  Sample.freeze_energy()
  Sample.like_init()
  plt.figure()
  plt.subplot(211)
  plt.title(tag)
  percentile_contour(*density_of_points(np.array([Sample.data['E'], Sample.data['theta']])), percents=percents) #, method='hist', nbin=20)
  #plt.subplot(212)
  #percentile_contour(np.array([np.log10(Sample.data['L2']), Sample.data['theta']]), percents=percents) #, method='hist', nbin=20)

##############################
n=int(sum(FullSampleFoF.data['E']>0)*1.0/sum(FullSampleAll.data['E']>0)*sum(All.data['E']>0))
nbin=50
method='hist'
Sample=All
plt.figure()
#plt.subplot(211)
f=Sample.data['E']>0
X,Y,Z=density_of_points(np.array([np.log10(Sample.data['E'][f]), Sample.data['theta'][f]]), bins=nbin, method=method)
extent=get_extent(X,Y)
bins=[np.linspace(extent[0], extent[1], nbin+1), np.linspace(extent[2], extent[3], nbin+1)]
#plt.imshow(Z, extent=extent)
#h,h0,lvls=percentile_contour(X,Y,Z, percents=percents, linestyles='dashed', alpha=0.5)
cs=plt.contour(X,Y,Z, linestyles='dashed', alpha=0.5);lvls=cs.levels
f=FoF.data['E']>0
plt.contour(*density_of_points(np.array([np.log10(FoF.data['E'][f]), FoF.data['theta'][f]])[:,0:n], bins=bins, method=method),
			levels=lvls, linestyles='solid')
for i in range(300):
	plt.plot(np.log10(Sublist.data['E'][i]), Sublist.data['theta'][i], 'o', markerfacecolor='w', markersize=40./(i+1.)**0.6)
plt.axis(extent)


lib.free_potential_spline()
#FullSample.clean()
#lib.close()