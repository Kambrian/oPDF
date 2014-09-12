#using nested_views rather than jointEL
from dynio import *
from myutils import *
import h5py,os,sys
from scipy.stats import chi2
plt.ion()

halo=sys.argv[1] #'A4'
npart=int(1e5)

lib.open()
FullSampleAll=Tracer(halo+'N',DynRMIN=1.1,DynRMAX=499)
All=FullSampleAll.copy(0,npart)
FullSampleAll.clean()

FoF=All.select(All.data['haloid']==0)
Bk=All.select(All.data['haloid']!=0) #ungrouped particle dominate over other halo particles.
#Other=All.select(All.data['haloid']>0)

r=np.logspace(0,np.log10(500),50)
bid=np.digitize(FoF.data['r'], r)
rFoF=np.array([(FoF.data['x']/FoF.data['r'].reshape([-1,1]))[bid==i].mean(axis=0) for i in xrange(len(r)-1)])
bid=np.digitize(All.data['r'], r)
rAll=np.array([(All.data['x']/All.data['r'].reshape([-1,1]))[bid==i].mean(axis=0) for i in xrange(len(r)-1)])
plt.plot(r[:-1], np.sqrt(np.sum(rFoF**2,axis=1)), 'r-')#mean value of directional vectors in each bin
plt.plot(r[:-1], np.sqrt(np.sum(rAll**2,axis=1)), 'g--')
plt.legend(('FoF','All'),loc='upper center')
plt.xscale('log')
plt.xlabel('r')
plt.ylabel(r'$|<\vec{e_r}>|$')
plt.ylim([0,0.5])
plt.savefig(lib.rootdir+'/plots/paper/extra/FoFDirectionProf.'+halo+'.eps')#, rasterize=True, dpi=300)

plt.figure()
bins=np.linspace(-500,500, 100)
x,y,z=density_of_points(All.data['x'][:,:2].T, bins, 'hist')
percentile_contour(x,y,z, percents=[0.1,0.3,0.5,0.7,0.9])
plt.title(halo)
plt.savefig(lib.rootdir+'/plots/paper/extra/DensityMap.'+halo+'.eps')#, rasterize=True, dpi=300)
#plt.contour(*density_of_points(Other.data['x'][:,:2].T, bins, 'kde'), levels=cs.levels, cmap=plt.cm.summer)
#plt.plot(Bk.data['x'][:,0], Bk.data['x'][:,1], 'r.', alpha=0.5)
#plt.gray()
#percentile_contour(*density_of_points(Bk.data['x'][:,:2].T, bins, 'kde'), percents=[0.1,0.3,0.5,0.7])
#h,h0,_=percentile_contour(*density_of_points(FoF.data['x'][:,:2].T, bins, 'hist'), percents=[0.1,0.3,0.5,0.7,0.9, 1], colors='r', alpha=0.5)
#lib.free_potential_spline()


lib.close()