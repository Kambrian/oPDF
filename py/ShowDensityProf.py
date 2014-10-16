from dynio import *
from myutils import *
import h5py,os,sys
from scipy.stats import chi2
plt.ion()

lib.open()

nbin=100
xbin=np.logspace(np.log10(1), np.log10(500), nbin)
xcen=(xbin*np.sqrt(xbin[1]/xbin[0]))[1:]
vol=np.diff(xbin**3)*np.pi*4/3

plt.figure()
for i,(halo,color) in enumerate([('A2','r'),('B2','g',),('C2','b'),('D2','c'),('E2','m')]):
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
  densNFW2=Halo.halo.Rhos/(xcen/Rs)/(1+xcen/Rs)**2
  lib.init_potential_spline()  
  Halo.define_halo([1,1])
  pot=-np.array(map(Halo.pot,xbin))
  mass=np.array(map(Halo.mass,xbin))
  dens=np.diff(mass)/vol
  lib.free_potential_spline()
  #iref=-1
  plt.subplot(121)
  iref=np.abs(xbin-Rv).argmin()
  plt.plot(xbin, potNFW/(pot-pot[iref]+potNFW[iref]), color=color)
  #plt.plot(xbin, massNFW/mass, color=color)
  #plt.plot(xbin, mass/xbin, 'x', color=color);plt.plot(xbin,massNFW/xbin,'-',color=color)
  #plt.plot(xbin, potNFW, 'k',xbin, (pot-pot[iref]+potNFW[iref]), 'gx')
  plt.semilogx()
  plt.xlabel('R')
  plt.ylabel(r'$\psi$')
  #plt.legend(('Data','NFWfit'))
  plt.ylim([0.8,1.3])

  plt.subplot(122)#consider subhalo removed smooth potential???????????????????????????????????????????????????????
  #plt.plot(xcen, densNFW/dens, color=color)
  iref=nbin/2
  plt.plot(xcen, 10**(-0.4*i)/dens[iref]*dens*xcen**2, 'x', color=color);plt.plot(xcen,10**(-0.4*i)/dens[iref]*densNFW*xcen**2,'-',color=color)
  #plt.plot(xcen, densNFW2, 'r--')
  plt.loglog()
  #plt.ylim([0.5,2])
  plt.xlabel('R')
  plt.ylabel(r'$\rho$')
  #plt.legend(('Data','NFWfit'))
plt.subplot(121)
plt.legend(('A','B','C','D','E'))  