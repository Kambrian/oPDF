#scan with min dist estimator
import matplotlib
matplotlib.use('Agg')
from dynio import *
from myutils import *
import h5py,os,sys
from scipy.stats import chi2

estimator=int(sys.argv[1]) #8 or 10

halo='Mock'
npart=1000 #number of particles
nbin=10 #E and L bin, or nbin**2 for single proxies
nx=50 #scan grid
x=np.logspace(-0.5,0.5,nx)

outdir=lib.rootdir+'/plots/scan'+halo+'%dZoom/'%npart
if not os.path.exists(outdir):
  os.makedirs(outdir)
  
name=lib.NameList[estimator]+'DistPotsRs'
outfile=outdir+name.replace('|','_')+'.hdf5'
f=h5py.File(outfile,'w')
mm=[m for m in x for c in x] #the first for is top layer, the second nested, so c varies first
mm=np.log10(np.array(mm).reshape([nx,nx], order="F"))
cc=[c for m in x for c in x]
cc=np.log10(np.array(cc).reshape([nx,nx], order="F"))
#if 'logm' not in f.keys():
f.create_dataset('/logm',data=mm)
f.create_dataset('/logc',data=cc)

lib.open()
FullSample=Tracer(halo)
Sample=FullSample.copy(5000,npart)
#Sample.radial_count(30)

def im_plot(infile, flagsave=True):
  ''' nbin=nbinE*nbinL: total number of bins used '''
  f=h5py.File(infile,'r')
  mm=f['/logm'][...]
  cc=f['/logc'][...]
  dm=mm[0,1]-mm[0,0]
  dc=cc[1,0]-cc[0,0]
  extent=[mm.ravel().min()-dm/2,mm.ravel().max()+dm/2,cc.ravel().min()-dc/2,cc.ravel().max()+dc/2]
  plt.figure()
  sig=f['/sig_dist'][...]
  plt.imshow(sig, extent=extent, cmap=plt.cm.summer) #the true extent should b
  lvls=[0.1,1,2,3,4]
  cs=plt.contour(mm,cc,sig, levels=lvls)
  plt.clabel(cs, inline=1, fmt=r'%.0f$\sigma$')
  #plt.contour(mm,cc,y, levels=[sig.min()+0.1], colors='y')
  xmin=f['/sig_dist'].attrs['xmin']
  plt.plot(xmin[0],xmin[1],'rx')
  print sig.argmin()
  print mm.ravel()[sig.ravel().argmin()],cc.ravel()[sig.ravel().argmin()]
  plt.plot(mm.ravel()[sig.ravel().argmin()], cc.ravel()[sig.ravel().argmin()],'ro') #this is slightly offset
  
  sig=f['/sig_like'][...]
  plt.imshow(sig, extent=extent, cmap=plt.cm.summer) #the true extent should b
  lvls=[0.1,1,2,3,4]
  cs=plt.contour(mm,cc,sig, levels=lvls, linestyles='dashed')
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
			  
y=[Sample.freeze_and_like([m,c], estimator=estimator) for m in x for c in x]
y=np.array(y).reshape([nx,nx], order="F")
if estimator==8:
  sig=AD2Sig(-y)
elif estimator==10:
  sig=P2Sig(chi2.sf(-y, 1))
else:
  print "wrong estimator=%d"%estimator
  raise estimator

dset=f.create_dataset('/ts',data=y)
dset.attrs['nbinE']=1
dset.attrs['nbinL']=1
f.create_dataset('/sig_dist',data=sig)
xmin=Sample.gfmin_dist(estimator,[1,1])
print xmin
print np.log10(xmin[0])
f['/sig_dist'].attrs['xmin']=np.log10(xmin[0])
siglike=np.array([lib.like_to_chi2(x, estimator) for x in y.flat]).reshape(y.shape)
siglike=siglike-siglike.ravel().min()
siglike=P2Sig(chi2.sf(siglike, 2)) #likeratio
f.create_dataset('/sig_like', data=siglike)
xmin=Sample.gfmin_like(estimator,[1,1])
print xmin
print np.log10(xmin[0])
f['/sig_like'].attrs['xmin']=np.log10(xmin[0])
f.close()

im_plot(outfile, name)

FullSample.clean()
Sample.clean()
lib.close()