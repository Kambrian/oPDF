import matplotlib
#matplotlib.use('Agg')
from dynio import *
from myutils import *
import h5py,os,sys
from scipy.stats import chi2
plt.ion()



halo='A2'
npart=int(1e6)

lib.open()

nbin=100
xbin=np.logspace(np.log10(1), np.log10(300), nbin)
xcen=(xbin*np.sqrt(xbin[1]/xbin[0]))[1:]
vol=np.diff(xbin**3)*np.pi*4/3

sample=Tracer(halo,DynDataFile='Null.hdf5')
sample.clean()
Halo=NFWHalo()
Halo.define_halo([1,1])
Rs=Halo.halo.Rs
Rv=Halo.halo.Rv
Rhos=Halo.halo.Rhos
potNFW=-np.array(map(Halo.pot,xbin))
massNFW=np.array(map(Halo.mass, xbin))
densNFW=np.diff(massNFW)/vol
Halo.define_halo([0.919,0.688])
Rs2=Halo.halo.Rs
Rv2=Halo.halo.Rv
Rhos2=Halo.halo.Rhos
potNFW2=-np.array(map(Halo.pot,xbin))
massNFW2=np.array(map(Halo.mass, xbin))
densNFW2=np.diff(massNFW2)/vol
lib.init_potential_spline()  
Halo.define_halo([1,1])
pot=-np.array(map(Halo.pot,xbin))
mass=np.array(map(Halo.mass,xbin))
dens=np.diff(mass)/vol
lib.free_potential_spline()
#iref=-1
plt.figure()
iref=np.abs(xbin-Rv).argmin()
plt.plot(xbin, potNFW/(pot-pot[iref]+potNFW[iref]), 'r-')
plt.plot(xbin, potNFW2/(pot-pot[iref]+potNFW2[iref]), 'g--')
#plt.plot(xbin, massNFW/mass, color=color)
#plt.plot(xbin, mass/xbin, 'x', color=color);plt.plot(xbin,massNFW/xbin,'-',color=color)
#plt.plot(xbin, potNFW, 'k',xbin, (pot-pot[iref]+potNFW[iref]), 'gx')
#plt.semilogx()
plt.plot(plt.xlim(), [1,1], 'k:')
plt.xlabel(r'$R$[kpc]')
plt.ylabel(r'$\psi/\psi_0$')
plt.legend(('NFW:True params','NFW:Best fit'),loc=4)
#plt.ylim([0.8,1.3])
plt.xlim([1,300])
plt.savefig(lib.rootdir+'/plots/paper/extra/A2Pot.eps')


with Tracer(halo) as F:
  with F.select(F.data['subid']<=0) as S:
	#S.mP*=float(F.nP)/S.nP
	countM,tmp=np.histogram(S.data['r'], xbin)#dM
	densClean=countM*S.mP/vol
	
plt.figure()
#plt.plot(xcen, densNFW/dens, color=color)
plt.plot(xcen, dens*xcen**2, 'x', color='k');
#plt.plot(xcen, densClean*xcen**2, 'o', color='c');
plt.plot(xcen,densNFW*xcen**2,'-',color='r')
plt.plot(xcen,densNFW2*xcen**2,'--',color='g')
#plt.plot(xcen, densNFW2, 'r--')
plt.loglog()
#plt.ylim([0.5,2])
plt.xlabel(r'$R$[kpc]')
plt.ylabel(r'$\rho R^2 [10^{10}M_\odot/\rm{kpc}]$')
plt.legend(('Data', 'NFW:True params','NFW:Best fit'),loc=2)
plt.savefig(lib.rootdir+'/plots/paper/extra/A2Dens.eps')




	