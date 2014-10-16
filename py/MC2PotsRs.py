#using nested_views rather than jointEL
from dynio import *
from myutils import *
import h5py,os,sys
from scipy.stats import chi2
plt.ion()

infile=lib.rootdir+'/plots/paper/extra/'+sys.argv[1]
#infile=lib.rootdir+'/plots/paper/extra/RBin50/fit.dat'
usetemplate=int(sys.argv[2])
outfile=infile+'.PsiR'
data=np.loadtxt(infile,delimiter=',')
data=data.reshape(5,-1,2) #halo, model, m-c

lib.open()
Halo=NFWHalo()
for haloid,halo in enumerate('ABCDE'):
  lib.load_config(halo+'2')
  lib.init_potential_spline()
  for i in xrange(data.shape[1]):
	data[haloid, i, :]=Halo.halo.mc2PotR(data[haloid, i, :], usetemplate)
  lib.free_potential_spline()
lib.close()
np.savetxt(outfile, data.reshape(-1,2), delimiter=',')