#using nested_views rather than jointEL
from dynio import *
from myutils import *
import h5py,os,sys,itertools
from scipy.stats import chi2
plt.ion()

nbin=100

colors=itertools.cycle(plt.cm.jet(np.linspace(0.,1.,6)))

lib.open()
H=NFWHalo()
h=[]
for halo in 'ABCDE':
  with Tracer(halo+'2', DynDataFile=halo+'2star.hdf5') as F:
	with F.select(F.data['subid']<=0) as S:
	  rmax=300
	  xbin=np.logspace(np.log10(0.5), np.log10(rmax), nbin)
	  countM,tmp=np.histogram(S.data['r'], np.hstack([0., xbin]))#dM
	  Mstar=countM.cumsum()
	  
	  lib.init_potential_spline()
	  Mdm=np.array(map(H.mass, xbin))
	  lib.free_potential_spline()
	  
	  color=colors.next()
	  htmp,=plt.plot(xbin, Mstar/float(Mstar[-1]), '-', color=color)
	  plt.plot(xbin, Mdm/Mdm[-1],'--', color=color)
	  h.append(htmp)

lib.close()
plt.xscale('log')
plt.xlabel(r'r/kpc')
plt.ylabel(r'$f(<r)$')
plt.legend(h,list('ABCDE'),loc=0)
plt.savefig(lib.rootdir+'/plots/paper/StarMassProfile.eps')