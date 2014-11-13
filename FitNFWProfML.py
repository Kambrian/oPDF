import matplotlib
#matplotlib.use('Agg')
from dynio import *
from myutils import *
import h5py,os,sys
from scipy.stats import chi2
plt.ion()


npart=int(1e6)

for halo in 'ABCDE':
  with Tracer(halo+'2') as F:
	with F.copy(0,npart) as S:
	  S.mP*=float(F.nP)/S.nP
	  result,m=S.minuit_NFWlike()
	  print 'Halo'+halo, m.values['m'], m.errors['m'], m.values['c'], m.errors['c']
