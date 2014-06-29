#scan the likelihoods
import matplotlib
matplotlib.use('Agg')
from dynio import *
from myutils import *
import h5py,os,sys
from scipy.stats import chi2

estimator=int(sys.argv[1])
proxy='' #one of '','E','L','LE'
if len(sys.argv)==3:
  proxy=sys.argv[2]

halo='Mock'
npart=1000 #number of particles
nbin=10 #E and L bin, or nbin**2 for single proxies
nx=30 #scan grid
x=np.logspace(-0.3,0.3,nx)
xx=x
#xx=np.logspace(-0.3,0.3,nx)
#for RBin estimator
FlagRBinLog=1
nbin_r=50
lib.NumRadialCountBin.value=nbin_r

outdir=lib.rootdir+'/plots/scan'+halo+'%dZoom/'%npart
if not os.path.exists(outdir):
  os.makedirs(outdir)
  
name=lib.NameList[estimator]
if estimator==4:
  namelog={0:'Lin', 1:'Log'}
  name+=namelog[FlagRBinLog]+'%d'%nbin_r
elif estimator==0:
  name+='_condition'
if proxy!='':
  name+='|'+proxy  
outfile=outdir+name.replace('|','_')+'.hdf5'
f=h5py.File(outfile,'w')
mm=[m for m in x for c in xx] #the first for is top layer, the second nested, so c varies first
mm=np.log10(np.array(mm).reshape([nx,nx], order="F"))
cc=[c for m in x for c in xx]
cc=np.log10(np.array(cc).reshape([nx,nx], order="F"))
#if 'logm' not in f.keys():
f.create_dataset('/logm',data=mm)
f.create_dataset('/logc',data=cc)

lib.open()
FullSample=Tracer(halo)
Sample=FullSample.copy(5000,npart)
if estimator==4:
  Sample.radial_count(nbin_r, FlagRBinLog)

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
			  
if proxy=='':
  nbinE=1
  nbinL=1
  if estimator==0:
	like=lambda x: -2*Sample.wenting_like_conditional(x)
  elif estimator==4:
	like=lambda x: -2*Sample.freeze_and_like(x, estimator)
  else:
	print "error: use scanMinDist instead."
	raise estimator
elif proxy=='E':
  nbinE=nbin**2
  nbinL=1  
  like=lambda x: Sample.jointE_FChi2(x, nbinE=nbinE, estimator=estimator)
elif proxy=='L':
  nbinE=1
  nbinL=nbin #*nbin  
  like=lambda x: Sample.jointLE_FChi2(x, nbinL=nbinL, nbinE=nbinE, estimator=estimator)
elif proxy=='LE' or proxy=='EL':
  nbinE=nbin
  nbinL=nbin  
  like=lambda x: Sample.jointLE_FChi2(x, nbinL=nbinL, nbinE=nbinE, estimator=estimator)
else:
  print "error: unknown proxy", proxy
  raise proxy

y=[like([m,c]) for m in x for c in xx]
xmin=fmin(like, [1,1], xtol=0.001, ftol=1e-4, maxiter=1000, maxfun=5000, full_output=True)
print min(y),xmin[1]
fval=min(min(y),xmin[1])	
y=np.array(y).reshape([nx,nx], order="F")-fval
dset=f.create_dataset('/ts',data=y)
dset.attrs['nbinE']=nbinE
dset.attrs['nbinL']=nbinL
sig=P2Sig(chi2.sf(y,2))
dset=f.create_dataset('/sig_like',data=sig)
dset.attrs['xmin']=np.log10(xmin[0])
f.close()
im_plot(outfile)

FullSample.clean()
Sample.clean()
lib.close()