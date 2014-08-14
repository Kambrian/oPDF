#using nested_views rather than jointEL
from dynio import *
from myutils import *
import h5py,os,sys
from scipy.stats import chi2
plt.ion()

#estimator=8
#proxy='LE'
#nbins=[10,10]
npart=1000

lib.open()
FullSample=Tracer('Mock')
Sample=FullSample.copy(5000,npart)


def plot_mass(pars, linestyle='-'):
  print pars
  r=np.logspace(0,2,20)
  Halo=NFWHalo()
  Halo.define_halo(pars)
  y=np.array(map(Halo.mass, r))
  plt.loglog(r, y, linestyle)

plot_mass([0.6013295 ,  2.22171445], 'r-')
plot_mass([1,1], 'k-')
plot_mass([2.16523294471, 0.476724073476], 'g-')
plt.legend((r'$0.6M_{\rm true},\,2.2c_{\rm true}$',r'$M_{\rm true},\,c_{\rm true}$',r'$2.2M_{\rm true},\,0.5c_{\rm true}$'), loc='lower right')
plt.xlabel(r'$R$[kpc]')
plt.ylabel(r'$M(<R)[10^{10} M_\odot]$')

#plot_mass(10**np.array([-0.0102, -0.03025]), 'k:')
#plot_mass(10**np.array([0.0178,0.0102]), 'k--')
#plot_mass(10**np.array([-0.2393, 0.3367]), 'r:')
#plot_mass(10**np.array([0.3163, -0.3385]), 'g:')


plot_mass(10**np.array([-0.0306, -0.0714]), 'k:')
plot_mass(10**np.array([-0.2755, 0.3142]), 'r:')
plot_mass(10**np.array([0.2959, -0.3859]), 'g:')
plt.xlim([10,100])
plt.ylim([4,100])
r=NFWHalo().Rs0.value
plt.plot([r, r], plt.ylim(), 'c--')
plt.text(r*1.05, 80, r'$r_s$')
r=np.median(Sample.data['r'])
plt.plot([r,r], plt.ylim(), 'c-')
plt.text(r*1.05, 80, r'$r_{\rm half}$')
plt.savefig(lib.rootdir+'/plots/paper/MassDegeneracy.eps')




  
Sample.clean()
FullSample.clean()
lib.close()