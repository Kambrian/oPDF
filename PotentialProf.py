#using nested_views rather than jointEL
from dynio import *
from myutils import *
import h5py,os,sys
from scipy.stats import chi2
plt.ion()

halo=sys.argv[1] #'B4'
nbin=100
npart=0 #int(1e6)

G=43007.1
lib.open()
FullSample=Tracer(halo+'N',DynRMIN=0,DynRMAX=500)
Sample=FullSample.copy(0,npart)
Sample.mP=FullSample.mP*FullSample.nP/Sample.nP
FullSample.clean()

#lib.init_potential_spline()
xbin=np.logspace(np.log10(0.5), np.log10(500), nbin)
xcen=xbin/np.sqrt(xbin[1]/xbin[0])
vol=np.diff(np.hstack([0., xbin])**3)*np.pi*4/3

countM,tmp=np.histogram(Sample.data['r'], np.hstack([0., xbin]))#dM
countR,tmp=np.histogram(Sample.data['r'], np.hstack([0., xbin]), weights=1./Sample.data['r'])#dM/R
density=countM*Sample.mP/vol
pot=countM.cumsum()/xbin+countR.sum()-countR.cumsum()
pot*=G*Sample.mP
density_cum=countM.cumsum()/xbin**3/(4*np.pi/3)*Sample.mP
#pad with bin 0
xbin=np.hstack([0., xbin])
pot=np.hstack([countR.sum()*G*Sample.mP, pot])
density_cum=np.hstack([density_cum[0], density_cum])


c0=get_config(halo+'N')
Halo=NFWHalo()
Halo.M0.value=float(c0['DynM0'])
Halo.C0.value=float(c0['DynC0'])
Halo.define_halo([1,1])
potNFW=-np.array(map(Halo.pot,xbin))

#iref=-1
iref=np.abs(xbin-Halo.halo.Rs).argmin()
plt.plot(xbin, pot-pot[iref]+potNFW[iref], 'gx')

plt.plot(xbin, potNFW, 'k')
plt.loglog()

plt.xlabel('R')
plt.ylabel(r'$\psi$')
plt.legend(('Data','NFWfit'))
#plt.savefig(lib.rootdir+'/plots/paper/extra/DensityProf'+halo+'ALLvsFoF.eps') #rasterize=True, dpi=300
print 'R'
print ','.join(['{:f}'.format(i) for i in xbin])
print 'Pot'
print ','.join(['{:f}'.format(i) for i in pot])
print 'AvDensity'
print ','.join(['{:g}'.format(i) for i in density_cum])

lib.init_potential_spline()
xnew=np.logspace(-1,3,50)
Halo.define_halo([1,2])
potNew=-np.array(map(Halo.pot, xnew))
Halo.define_halo([2,1])
potNew2=-np.array(map(Halo.pot, xnew))
plt.figure()
plt.plot(xnew, potNew, 'ro')
plt.plot(xnew, potNew2, 'gs')
plt.plot(xbin, pot, 'k-')
plt.loglog()
#Sample.clean()
#lib.free_potential_spline()
#lib.close()