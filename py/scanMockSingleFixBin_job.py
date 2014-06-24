#using nested_views rather than jointEL
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
binpar=[1.,1.] #m,c values to freeze the bins
npart=1000 #number of particles
nx=20 #scan grid
x=np.logspace(-0.3,0.3,nx)
if proxy=='E':
  nbinE=nbin*nbin
  nbinL=1  
  nbins=[nbinE]
else:
  nbinE=nbin
  nbinL=nbin
  nbins=[nbin,nbin]

outdir=lib.rootdir+'/plots/scan'+halo+'%dZoom/'%npart
if not os.path.exists(outdir):
  os.makedirs(outdir)
  
name=lib.NameList[estimator]
if proxy!='':
  name+='|'+proxy+'fix_%.1f_%.1f'%tuple(binpar)
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
Sample=FullSample.copy(0,npart)
FullSample.clean()
Sample.create_nested_views(binpar, proxy, nbins)

def im_plot(mm,cc,y,nbin, name, flagsave=True):
  ''' nbin=nbinE*nbinL: total number of bins used '''
  dm=mm[0,1]-mm[0,0]
  dc=cc[1,0]-cc[0,0]
  extent=[mm.ravel().min()-dm/2,mm.ravel().max()+dm/2,cc.ravel().min()-dc/2,cc.ravel().max()+dc/2]
  plt.figure()
  plt.imshow(y, extent=extent, cmap=plt.cm.summer) #the true extent should b
  cs=plt.contour(mm,cc,y)
  plt.clabel(cs, inline=1)
  percents=[0.683, 0.954, 0.997]
  levels=chi2.ppf(percents, 2) #always use like-rat error def
  y=y-y.min()
  cs=plt.contour(mm,cc,y, levels=levels, colors='r')
  plt.clabel(cs, inline=1, fmt={levels[0]:r'1$\sigma$', levels[1]:r'2$\sigma$', levels[2]:r'3$\sigma$'})
  print y.argmin()
  print mm.ravel()[y.ravel().argmin()],cc.ravel()[y.ravel().argmin()]
  plt.plot(mm.ravel()[y.ravel().argmin()], cc.ravel()[y.ravel().argmin()],'ro') #this is slightly offset
  plt.title(name)
  plt.xlabel(r'$\log(M/M_0)$')
  plt.ylabel(r'$\log(c/c_0)$')
  if flagsave:
	plt.savefig(outdir+name.replace('|','_')+'.eps')

y=[Sample.nested_views_Flike([m,c], estimator) for m in x for c in x]
y=np.array(y).reshape([nx,nx], order="F")
dset=f.create_dataset('/ts',data=y)
dset.attrs['proxies']=proxy
dset.attrs['nbinE']=nbinE
dset.attrs['nbinL']=nbinL
f.close()
im_plot(mm,cc,y, nbinE*nbinL, name)

Sample.clean()
lib.close()