import matplotlib
#matplotlib.use('Agg')
from dynio import *
from myutils import *
import h5py,os,sys
from scipy.stats import chi2
plt.ion()

def plot_star(halo,dims=[0,1]):
  with Tracer(halo,DynDataFile=halo+'star.hdf5') as F:
	f=F.data['subid']<=0
	bins=np.linspace(-300,300,100)
	percents=[0.3,0.5,0.9] #np.arange(0.5,0.92,0.2)
	plt.figure()
	x,y,z=density_of_points(F.data['x'][f].take(dims, axis=1), bins, method='hist')
	lvls=percent2level(percents,z)
	cs=plt.contour(x,y,z, lvls)
	#plt.clabel(cs, inline=1, fmt=dict(zip(lvls, map(str, percents))))
	plt.plot(F.data['x'][~f,dims[0]], F.data['x'][~f,dims[1]], 'r.')
	#x,y,z=density_of_points(F.data['x'][~f,:2], bins, method='hist')
	#percentile_contour(x,y,z, percents, colors='g')
	plt.title(halo)
	plt.xlabel('x/kpc')
	plt.ylabel('y/kpc')
	plt.savefig(lib.rootdir+'/plots/paper/extra/Haloes/'+halo+'starMap.eps', rasterize=True)

for h in 'ABCD':
  plot_star(h+'2')

plot_star('E2',[1,2])  