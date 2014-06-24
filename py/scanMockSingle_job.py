#scan with binned min dist estimator
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
nx=20 #scan grid
x=np.logspace(-0.2,0.2,nx)

outdir=lib.rootdir+'/plots/scan'+halo+'%dZoom/'%npart
if not os.path.exists(outdir):
  os.makedirs(outdir)
  
name=lib.NameList[estimator]
if proxy!='':
  name+='|'+proxy
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
#Sample.radial_count(30)

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
  levels=chi2.ppf(percents, nbin)
  try:
	if name[0:4] in ['ADGE','ADNM','ADBN']:
	  levels=chi2.ppf(percents, 2)
	  y=y-y.min()
	  print name
  except:
	pass
  cs=plt.contour(mm,cc,y, levels=levels, colors='r')
  plt.clabel(cs, inline=1, fmt={levels[0]:r'1$\sigma$', levels[1]:r'2$\sigma$', levels[2]:r'3$\sigma$'})
  plt.contour(mm,cc,y, levels=[y.min()+0.1], colors='y')
  print y.argmin()
  print mm.ravel()[y.ravel().argmin()],cc.ravel()[y.ravel().argmin()]
  plt.plot(mm.ravel()[y.ravel().argmin()], cc.ravel()[y.ravel().argmin()],'ro') #this is slightly offset
  plt.title(name)
  plt.xlabel(r'$\log(M/M_0)$')
  plt.ylabel(r'$\log(c/c_0)$')
  if flagsave:
	plt.savefig(outdir+name.replace('|','_')+'.eps')
			  
if proxy=='':
  nbinE=1
  nbinL=1
  y=[lib.like_to_chi2(Sample.freeze_and_like([m,c], estimator=estimator),estimator) for m in x for c in x]
elif proxy=='E':
  nbinE=nbin**2
  nbinL=1  
  y=[Sample.jointE_Flike([m,c], nbinE=nbinE, estimator=estimator) for m in x for c in x]  
elif proxy=='L':
  nbinE=1
  nbinL=nbin**2  
  y=[Sample.jointLE_Flike([m,c], nbinL=nbinL, nbinE=nbinE, estimator=estimator) for m in x for c in x]  
elif proxy=='LE' or proxy=='EL':
  nbinE=nbin
  nbinL=nbin  
  y=[Sample.jointLE_Flike([m,c], nbinL=nbinL, nbinE=nbinE, estimator=estimator) for m in x for c in x]
else:
  print "error: unknown proxy", proxy
  raise proxy  
y=np.array(y).reshape([nx,nx], order="F")
dset=f.create_dataset('/ts',data=y)
dset.attrs['nbinE']=nbinE
dset.attrs['nbinL']=nbinL
f.close()
im_plot(mm,cc,y, nbinE*nbinL, name)

FullSample.clean()
Sample.clean()
lib.close()