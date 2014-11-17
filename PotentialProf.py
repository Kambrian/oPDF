""" This file creates the potential and cumulative density profile templates for a given DM file"""
from oPDF import *
from myutils import *
import h5py,os,sys
from scipy.stats import chi2
plt.ion()


DMfile=rootdir+'/../../data/A4DM.hdf5'
#real parameters, for comparison with analytical profile:
M0=183.8
C0=15.07

nbin=100
npart=0 #int(1e6)

FullSample=Tracer(DMfile,rmin=0,rmax=500)
Sample=FullSample.copy(0,npart)
Sample.mP=FullSample.mP*FullSample.nP/Sample.nP #rescale particle mass to account for selection
FullSample.clean()

xbin=np.logspace(np.log10(0.5), np.log10(500), nbin)
xcen=xbin/np.sqrt(xbin[1]/xbin[0])
vol=np.diff(np.hstack([0., xbin])**3)*np.pi*4/3

countM,tmp=np.histogram(Sample.data['r'], np.hstack([0., xbin]))#dM
countR,tmp=np.histogram(Sample.data['r'], np.hstack([0., xbin]), weights=1./Sample.data['r'])#dM/R
density=countM*Sample.mP/vol
pot=countM.cumsum()/xbin+countR.sum()-countR.cumsum() #M(<r)/r+\int_r^rmax dM/r
pot*=Globals.units.Const.G*Sample.mP
density_cum=countM.cumsum()/xbin**3/(4*np.pi/3)*Sample.mP
#pad with bin 0
xbin=np.hstack([0., xbin])
pot=np.hstack([countR.sum()*Globals.units.Const.G*Sample.mP, pot])
density_cum=np.hstack([density_cum[0], density_cum])

halo=Halo()
halo.set_param([M0,C0])
potNFW=-halo.pot(xbin)

#iref=-1
iref=np.abs(xbin-halo.Rs).argmin()
plt.plot(xbin, pot-pot[iref]+potNFW[iref], 'gx')
plt.plot(xbin, potNFW, 'k')
plt.loglog()

plt.xlabel('R')
plt.ylabel(r'$\psi$')
plt.legend(('Data','NFWfit'))
#plt.savefig(lib.rootdir+'/plots/paper/extra/DensityProf'+halo+'ALLvsFoF.eps') #rasterize=True, dpi=300
print 'Profile template to be added to C/TemplateData.h:'
print 'R'
print ','.join(['{:f}'.format(i) for i in xbin])
print 'Pot'
print ','.join(['{:f}'.format(i) for i in pot])
print 'AvDensity'
print ','.join(['{:g}'.format(i) for i in density_cum])

## Now recompile and try the newly added template
TMPid=0 #id of the newly added template 

xnew=np.logspace(-1,3,50)
tmphalo=Halo(halotype=HaloTypes.TMPMC, TMPid=0)
tmphalo.set_param([M0,C0])
potNew=-tmphalo.pot(xnew)
tmphalo2=Halo(halotype=HaloTypes.TMPMC, TMPid=0)
tmphalo2.set_param([2*M0,C0])
potNew2=-tmphalo2.pot(xnew)
plt.figure()
plt.plot(xnew, potNew, 'ro')
plt.plot(xnew, potNew2, 'gs')
plt.plot(xbin, pot, 'k-')
plt.loglog()
#Sample.clean()
#lib.free_potential_spline()
#lib.close()