import matplotlib
#matplotlib.use('Agg')
from dynio import *
from myutils import *
import h5py,os,sys
from scipy.stats import chi2
plt.ion()

halo=sys.argv[1]
estimator=8
npart=1e5

#lib.open()
with Tracer(halo) as Raw:
  with Raw.copy(0, npart) as F:
	F.freeze_energy()
	F.like_init()
	with F.select(F.data['subid']>0) as Sub:
	  with F.select(F.data['subid']<=0) as Clean:
		TSf=F.like_eval(estimator=estimator)
		TSs=Sub.like_eval(estimator=estimator)
		TSc=Clean.like_eval(estimator=estimator)
		TS=(TSs*np.sqrt(Sub.nP)+TSc*np.sqrt(Clean.nP))/np.sqrt(F.nP)

print halo, TSf, TSc, TSs, TS
