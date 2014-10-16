# todo: stream-cleaned version
# star 500kpc
from dynio import *
from myutils import *
import h5py,os,sys
from scipy.stats import chi2
plt.ion()

def plot_subhalos(sample, sub, proxy, nPmin=1000):
  '''plot subhaloes as circles in phase space'''
  with sub.select(sub.mP*sub.data['w']>nPmin) as S:
	p=np.sort(sample.data[proxy])
	x=np.searchsorted(p, S.data[proxy])/float(sample.nP)
	plt.scatter(x, S.data['theta'], marker='o', edgecolor='r', facecolor='none', s=(S.mP*S.data['w'])/300)
	
def phase_contour(sample, proxy, bins, method, logscale, levels=5, percents=None):
  '''percents will overide levels if present'''
  X,Y,Z,_,_=sample.phase_density(proxy, bins, method, logscale)
  #Z[Z==0]=np.nan;
  if percents is None:
	cs=plt.contour(X,Y,Z,levels)
	levels=cs.levels
  else:
	h,h0,levels=percentile_contour(X,Y,Z, percents=percents)
  #plt.axis('tight')
  return levels
  
def phase_imshow(sample, proxy, bins, method='hist'):
  X,Y,Z,extent,data=sample.phase_density(proxy, bins, method, False)
  Z[Z==0]=np.nan;
  plt.imshow(Z,extent=[0,1,0,1])
  plt.axis('tight')
  sk=skeleton(data[0],data[1],nbin=bins[0])
  nbin=len(bins[0])-1
  plt.plot((np.arange(nbin)+0.5)/nbin, sk['y']['mean'], 'w-')
  plt.plot([0,1], [0.5, 0.5], 'k-')

halo=sys.argv[1] #'A4'
npart=0 #int(1e6)
rmin=1
rmax=300

#lib.open()
with Tracer(halo,DynRMIN=rmin,DynRMAX=rmax,DynDataFile=halo+'DM.hdf5') as FullSample:
  Sample=FullSample.copy(0,npart)
Star=Tracer(halo, DynRMIN=rmin, DynRMAX=rmax, DynDataFile=halo+'star.hdf5')  
lib.init_potential_spline()

Sample.freeze_energy();
#Sample.like_init()
Star.freeze_energy();
#Sub.like_init()

#e=np.sort(Sample.data['E'])
#l=np.sort(Sample.data['L2'])
#ie=np.argsort(Sample.data['E'])
#il=np.argsort(Sample.data['L2'])
#x=np.searchsorted(e, Star.data['E'])
#y=np.searchsorted(l, Star.data['L2'])

#plt.plot(x,y,'.')
#bins=np.linspace(0,Sample.nP,100)
#X,Y,Z=density_of_points(np.array((Star.data['E'],Star.data['L2'])), bins=100, method='hist')
#plt.contour(Z,extent=[0,1,0,1])

TMP=Sample.select(Sample.data['E']>0)
Sample.clean()
Sample=TMP
TMP=Star.select(Star.data['E']>0)
Star.clean()
Star=TMP
nbin=20
ebins=np.logspace(np.log10(Sample.data['E'].min()),np.log10(Sample.data['E'].max()),nbin)
lbins=np.logspace(np.log10(Sample.data['L2'].min()),np.log10(Sample.data['L2'].max()),nbin)
Counts,x,y=np.histogram2d(Star.data['E'],Star.data['L2'], bins=[ebins,lbins])
ids=[]
np.random.seed(10)
for i in xrange(nbin-1):
  print i
  for j in xrange(nbin-1):
	cand=np.where((Sample.data['E']>=x[i])&(Sample.data['E']<x[i+1])&(Sample.data['L2']>=y[j])&(Sample.data['L2']<y[j+1]))
	if (Counts[i,j]>0)&(Counts[i,j]<=len(cand)):
	  ids.extend(cand[np.random.randint(0, len(cand), Counts[i,j])])
	elif Counts[i,j]>len(cand):
	  print 'Fail', Counts[i,j], len(cand)
f=np.zeros(Sample.nP, dtype=bool)
f[ids]=1
NewStar=Sample.select(f)
	
