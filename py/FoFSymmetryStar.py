#using nested_views rather than jointEL
from dynio import *
from myutils import *
import h5py,os,sys
from scipy.stats import chi2
plt.ion()

halo=sys.argv[1] #'A4'
npart=int(1e5)
rmax=50

lib.open()
with Tracer(halo,DynDataFile=halo+'star.hdf5') as F:
  F.select(F.data['subid']<=0)
  Star=F.copy(0,npart)
with Tracer(halo) as F:
  F.select(F.data['subid']<=0)
  DM=F.copy(0,npart)

r=np.logspace(0,np.log10(Star.rmax),50)
bid=np.digitize(Star.data['r'], r)
rStar=np.array([(Star.data['x']/Star.data['r'].reshape([-1,1]))[bid==i].mean(axis=0) for i in xrange(len(r)-1)])
bid=np.digitize(DM.data['r'], r)
rDM=np.array([(DM.data['x']/DM.data['r'].reshape([-1,1]))[bid==i].mean(axis=0) for i in xrange(len(r)-1)])
plt.plot(r[:-1], np.sqrt(np.sum(rStar**2,axis=1)), 'r-')#mean value of directional vectors in each bin
plt.plot(r[:-1], np.sqrt(np.sum(rDM**2,axis=1)), 'g--')
plt.legend(('Tag','DM'),loc='upper center')
plt.xscale('log')
plt.xlabel('r')
plt.ylabel(r'$|<\vec{e_r}>|$')
plt.ylim([0,0.5])
plt.savefig(lib.rootdir+'/plots/paper/extra/Haloes/DirectionProf.'+halo+'.tag.eps')#, rasterize=True, dpi=300)

plt.figure()
bins=np.linspace(-rmax,rmax, 100)
x,y,z=density_of_points(Star.data['x'][:,1:].T, bins, 'kde')
percentile_contour(x,y,z, percents=[0.1,0.3,0.5,0.7,0.9])
DM.radial_cut(1,rmax)
x,y,z=density_of_points(DM.data['x'][:,1:].T, bins, 'kde')
percentile_contour(x,y,z, percents=[0.1,0.3,0.5,0.7], colors=([1,0,0,0.3],), linestyles='solid')
plt.title(halo)
plt.savefig(lib.rootdir+'/plots/paper/extra/Haloes/DensityMap.'+halo+'.tag.yz.pdf')#, rasterize=True, dpi=300)
#plt.contour(*density_of_points(Other.data['x'][:,:2].T, bins, 'kde'), levels=cs.levels, cmap=plt.cm.summer)
#plt.plot(Bk.data['x'][:,0], Bk.data['x'][:,1], 'r.', alpha=0.5)
#plt.gray()
#percentile_contour(*density_of_points(Bk.data['x'][:,:2].T, bins, 'kde'), percents=[0.1,0.3,0.5,0.7])
#h,h0,_=percentile_contour(*density_of_points(Star.data['x'][:,:2].T, bins, 'hist'), percents=[0.1,0.3,0.5,0.7,0.9, 1], colors='r', alpha=0.5)
#lib.free_potential_spline()

Star.clean()
DM.clean()

lib.close()