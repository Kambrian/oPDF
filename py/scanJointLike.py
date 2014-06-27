#using nested_views to scan fixed-bin joint likelihood
import matplotlib
matplotlib.use('Agg')
from dynio import *
from myutils import *
import h5py,os,sys
from scipy.stats import chi2

try: 
  estimator=int(sys.argv[1]) #use 10 or 16
  proxy=sys.argv[2]
  if proxy not in ['EL','LE','E']:
	raise
except:
  print "Incorrect usage.\n Example: %s EL (or LE or E)\n Now exit."%sys.argv[0]
  raise

halo='Mock'
nbin=10
IterFit=True  #whether to use iterative fit
binpar=[2.,2.] #values to freeze energy if not IterFit; if IterFit then this is the initial value for the fit
npart=1000 #number of particles
nx=30 #scan grid
x=np.logspace(-0.3,0.3,nx)
if proxy=='E':
  nbinE=nbin*nbin
  nbinL=1  
  nbins=[nbinE]
elif proxy=='LE':
  nbinE=nbin
  nbinL=nbin
  nbins=[nbinL,nbinE]
elif proxy=='EL':
  nbinL=nbin
  nbinE=nbin
  nbins=[nbinE,nbinL]

outdir=lib.rootdir+'/plots/scan'+halo+'%dZoom/'%npart
if not os.path.exists(outdir):
  os.makedirs(outdir)
  
lib.open()
FullSample=Tracer(halo)
Sample=FullSample.copy(5000,npart)
FullSample.clean()

if IterFit:
  x0=Sample.fmin_FixBinIter(estimator, proxy, nbins, binpar, maxiter=10)
  binpar=x0[0]
Sample.create_nested_views(binpar, proxy, nbins)

like=lambda x: Sample.nested_views_FChi2(x, estimator)
y=[like([m,c]) for m in x for c in x]
xmin=fmin(like, [1,1], xtol=0.001, ftol=1e-4, maxiter=1000, maxfun=5000, full_output=True)
print min(y),xmin[1]
fval=min(min(y),xmin[1])	
y=np.array(y).reshape([nx,nx], order="F")-fval

name=lib.NameList[estimator]
if proxy!='':
  if IterFit:
	name+='Iter'+'|'+proxy
  else:
	name+='|'+proxy+'fix_%.1f_%.1f'%tuple(binpar)
	
outfile=outdir+name.replace('|','_')+'.hdf5'
f=h5py.File(outfile,'w')
mm=[m for m in x for c in x] #the first for is top layer, the second nested, so c varies first
mm=np.log10(np.array(mm).reshape([nx,nx], order="F"))
cc=[c for m in x for c in x]
cc=np.log10(np.array(cc).reshape([nx,nx], order="F"))
f.create_dataset('/logm',data=mm)
f.create_dataset('/logc',data=cc)
dset=f.create_dataset('/ts',data=y)
dset.attrs['proxies']=proxy
dset.attrs['nbinE']=nbinE
dset.attrs['nbinL']=nbinL
sig=P2Sig(chi2.sf(y,2))
dset=f.create_dataset('/sig_like',data=sig)
dset.attrs['xmin']=np.log10(xmin[0])
f.close()

def im_plot(infile, flagsave=True):
  f=h5py.File(infile,'r')
  mm=f['/logm'][...]
  cc=f['/logc'][...]
  dm=mm[0,1]-mm[0,0]
  dc=cc[1,0]-cc[0,0]
  extent=[mm.ravel().min()-dm/2,mm.ravel().max()+dm/2,cc.ravel().min()-dc/2,cc.ravel().max()+dc/2]
  plt.figure() 
  sig=f['/sig_like'][...]
  plt.imshow(sig, extent=extent, cmap=plt.cm.summer) #the true extent should b
  lvls=[0.1,1,2,3,4]
  cs=plt.contour(mm,cc,sig, levels=lvls)
  plt.clabel(cs, inline=1, fmt=r'%.0f$\sigma$')
  #plt.contour(mm,cc,y, levels=[sig.min()+0.1], colors='y')
  xmin=f['/sig_like'].attrs['xmin']
  plt.plot(xmin[0],xmin[1],'y+')
  print sig.argmin()
  print mm.ravel()[sig.ravel().argmin()],cc.ravel()[sig.ravel().argmin()]
  plt.plot(mm.ravel()[sig.ravel().argmin()], cc.ravel()[sig.ravel().argmin()],'ys') #this is slightly offset
  plt.title(name)
  plt.xlabel(r'$\log(M/M_0)$')
  plt.ylabel(r'$\log(c/c_0)$')
  plt.plot(plt.xlim(),[0,0],'k:',[0,0],plt.ylim(),'k:')
  if flagsave:
	plt.savefig(outdir+name.replace('|','_')+'.eps')
	
im_plot(outfile)

Sample.clean()
lib.close()