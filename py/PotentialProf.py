#using nested_views rather than jointEL
from dynio import *
from myutils import *
import h5py,os,sys
from scipy.stats import chi2
plt.ion()

halo='AqA2'
nbin=100

G=43007.1
lib.open()
FullSample=Tracer(halo,DynRMIN=0,DynRMAX=500)
#lib.init_potential_spline()
xbin=np.logspace(0, np.log10(500), nbin)
xcen=xbin/np.sqrt(xbin[1]/xbin[0])
vol=np.diff(np.hstack([0., xbin])**3)*np.pi*4/3

c0=get_config(halo+'N')
Halo=NFWHalo()
Halo.M0.value=float(c0['DynM0'])
Halo.C0.value=float(c0['DynC0'])
Halo.define_halo([1,1])
potNFW=-np.array(map(Halo.pot,xbin))

countM,tmp=np.histogram(FullSample.data['r'], np.hstack([0., xbin]))#dM
countR,tmp=np.histogram(FullSample.data['r'], np.hstack([0., xbin]), weights=1./FullSample.data['r'])#dM/R
density=countM*FullSample.mP/vol
pot=countM.cumsum()/xbin+countR.sum()-countR.cumsum()
pot*=G*FullSample.mP
#iref=-1
iref=np.abs(xbin-Halo.halo.Rs).argmin()
plt.plot(xbin, pot-pot[iref]+potNFW[iref], 'gx')

plt.plot(xbin, potNFW, 'k')
plt.loglog()

plt.xlabel('R')
plt.ylabel(r'$\psi$')
plt.legend(('Data','NFWfit'))
#plt.savefig(lib.rootdir+'/plots/paper/extra/DensityProf'+halo+'ALLvsFoF.eps') #rasterize=True, dpi=300
print ','.join(['{:f}'.format(i) for i in xbin])
print
print ','.join(['{:f}'.format(i) for i in pot])

#FullSample.clean()
#lib.free_potential_spline()
#lib.close()