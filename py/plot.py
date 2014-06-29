import sys,os
from dynio import *
os.environ['OMP_NUM_THREADS']='32'

#from matplotlib.pyplot import *
#import matplotlib.pyplot as plt
from numpy import *
from scipy.stats import chi2,norm
from myutils import *
import triangle
import glob,itertools
import h5py

plt.ion()

def fig_DegeneracyDistributionRad(nbin=50, logscale=True, cumulative=False, flagsave=False):
  npart=1000
  lib.open()
  FullSample=Tracer('Mock')
  Sample=FullSample.copy(5000,npart)
  FullSample.clean()
  if logscale:
	xbin=np.logspace(np.log10(Sample.rmin), np.log10(Sample.rmax), nbin+1)
  else:
	xbin=np.linspace(Sample.rmin, Sample.rmax, nbin+1)
  count,tmp=np.histogram(Sample.data['r'], xbin)
  
  plt.figure()
  pars=[1,1]
  Sample.freeze_energy(pars)
  Sample.like_init(pars)
  pred=Sample.predict_radial_count(nbin, logscale)
  if cumulative:
	plt.plot(xbin, np.hstack((0.,pred)).cumsum()/Sample.nP, 'r-')
  else:
	plt.step(xbin,hstack((pred,nan)), where='post', color='r')
	
  #pars=[ 0.82423562, 1.28484559]
  pars=[0.6013295 ,  2.22171445]
  Sample.freeze_energy(pars)
  Sample.like_init(pars)
  pred=Sample.predict_radial_count(nbin, logscale)
  if cumulative:
	plt.plot(xbin, np.hstack((0.,pred)).cumsum()/Sample.nP, 'g-')
  else:
	plt.step(xbin,hstack((pred,nan)), where='post', color='g')
  
  if cumulative:
	plt.plot(xbin, np.hstack((0.,count)).cumsum()/Sample.nP, 'k--')
	plt.ylabel(r'$P(<r)$')
	plt.legend((r'$M_{\rm true},\, c_{\rm true}$', r'$0.6 M_{\rm true},\, 2.2 c_{\rm true}$', 'Data'),loc='lower right')
	plt.ylim([0,1])
  else:
	plt.step(xbin,hstack((count,nan)), where='post', color='k', linestyle='--', linewidth=2)
	plt.ylabel('Counts')
	plt.legend((r'$M_{\rm true},\, c_{\rm true}$', r'$0.6 M_{\rm true},\, 2.2 c_{\rm true}$', 'Data'),loc='upper right')
	
  if logscale:
	  plt.xscale('log')
  
  plt.xlabel(r'$r/kpc$')
 
  Sample.clean()
  lib.close()
  if flagsave:
	plt.savefig(lib.rootdir+'/plots/paper/DegeneracyDistributionRad.eps')
	
def fig_DegeneracyDistribution(nbin=50, cumulative=True, flagsave=False):
  npart=1000
  lib.open()
  FullSample=Tracer('Mock')
  Sample=FullSample.copy(5000,npart)
  FullSample.clean()
  plt.figure()
  pars=[1,1]
  Sample.freeze_energy(pars)
  Sample.like_init(pars)
  if cumulative:
	t=Sample.data['theta'][...]
	t.sort()
	plt.plot(t, (np.arange(Sample.nP)+1.)/Sample.nP, 'b-')
  else:
	plt.hist(Sample.data['theta'], nbin, normed=True, histtype='step', color='b')
	
  #pars=[ 0.82423562, 1.28484559]
  pars=[0.6013295 ,  2.22171445]
  Sample.freeze_energy(pars)
  Sample.like_init(pars)
  if cumulative:
	t=Sample.data['theta'][...]
	t.sort()
	plt.plot(t, (np.arange(Sample.nP)+1.)/Sample.nP, 'r-')
  else:
	plt.hist(Sample.data['theta'], nbin, normed=True, histtype='step', color='r')
  
  plt.legend((r'$M_{\rm true},\, c_{\rm true}$', r'$0.6 M_{\rm true},\, 2.2 c_{\rm true}$'),loc='upper left')
  if cumulative:
	plt.plot([0,1],[0,1],'k:')
	plt.ylabel(r'$P(<\theta)$')
	plt.axis([0,1,0,1])
  else:
	plt.plot([0,1],[1,1],'k:')
	plt.ylabel(r'$dP/d\theta$')
  plt.xlabel(r'$\theta$')
  
  Sample.clean()
  lib.close()
  if flagsave:
	plt.savefig(lib.rootdir+'/plots/paper/DegeneracyDistribution.eps')
  
  
def fig_MinDistContour(flagsave=False):
  f1=h5py.File(lib.rootdir+'/plots/paper/raw/ADDist.hdf5','r')
  mm=f1['/logm']
  cc=f1['/logc']
  sig=f1['/sig_dist']
  x=f1['/sig_dist'].attrs['xmin']
  print 10**x
  cs=plt.contour(mm,cc,sig, levels=[1,3], linestyles='solid')
  #plt.clabel(cs, inline=1, fmt={1:'',2:'AD'})
  h1=Ellipse((0,0),0,0,fill=False, linestyle='solid')
  #x1=np.log10(np.array([0.80234718,  1.09236641])) #initial value [1,1]
  #x2=np.log10(np.array([ 0.65940932,  1.49481424])) #initial value [1.5,1.5]
  #plt.plot(x[0],x[1],'ro',markersize=10)
  #plt.plot(x2[0],x2[1],'gs',markersize=10)
  f2=h5py.File(lib.rootdir+'/plots/paper/raw/MeanDist.hdf5','r')
  mm=f2['/logm']
  cc=f2['/logc']
  sig=f2['/sig_dist']
  x=f2['/sig_dist'].attrs['xmin']
  cs=plt.contour(mm,cc,sig, levels=[1,3], linestyles='dashed')
  h2=Ellipse((0,0),0,0,fill=False, linestyle='dashed')
  #x1=np.log10(np.array([ 0.95707194,  1.0434845])) #initial value [1,1]
  x2=np.log10(np.array([ 0.6013295 ,  2.22171445])) #initial value [2,2]
  #plt.plot(x1[0],x1[1],'gx',markersize=10)
  plt.plot(x2[0],x2[1],'gx',markersize=10)
  #plt.clabel(cs, inline=1, fmt={1:'',2:'Mean'})
  #plt.plot(x[0],x[1],'g+',markersize=10)
  #plt.xlim([-0.3,0.8])
  #plt.ylim([-0.8,0.8])
  plt.plot(plt.xlim(),[0,0],'k:',[0,0],plt.ylim(),'k:')
  plt.xlabel(r'$\log(M/M_{\rm true})$')
  plt.ylabel(r'$\log(c/c_{\rm true})$')
  plt.legend((h1,h2),('AD','Mean'))
  if flagsave:
	plt.savefig(lib.rootdir+'/plots/paper/MinDistContour.eps')

def fig_MaxLikeContour(estimator='Mean',flagsave=False):
  '''estimator in 'AD' or 'Mean' '''
  def plot_siglike(name, color, linestyle='solid'):
	marker='o' #{'solid':'s','dashed':'d','dotted':'o','dashdot':'x'}[linestyle]
	f=h5py.File(lib.rootdir+'/plots/paper/raw/'+name+'.hdf5','r')
	cs=plt.contour(f['/logm'],f['/logc'],f['/sig_like'], levels=[1,], colors=color, linestyles=linestyle)
	h=Ellipse((0,0),0,0,fill=False, color=color, linestyle=linestyle)
	x=f['/sig_like'].attrs['xmin']
	print 10**x
	plt.plot(x[0],x[1],marker,color=color,markersize=6)
	f.close()
	return h
  plt.figure()
  f=h5py.File(lib.rootdir+'/plots/paper/raw/MeanDist.hdf5','r')
  plt.contourf(f['/logm'],f['/logc'],f['/sig_dist'], levels=[0,1], colors=((0.8,0.8,0.8),))
  #plt.contour(f['/logm'],f['/logc'],f['/sig_dist'], levels=[3], linestyles='dashed', colors=((0.2,0.2,0.2),))
  #h6=Ellipse((0,0),0,0,fill=False, color='y')
  f.close()
  f=h5py.File(lib.rootdir+'/plots/paper/raw/ADDist.hdf5','r')
  plt.contourf(f['/logm'],f['/logc'],f['/sig_dist'], levels=[0,1], colors=((0.7,0.7,0.7),)) #colors='k', alpha=0.2
  #h5=Ellipse((0,0),0,0,fill=False, color='m')
  f.close()
  #data=np.loadtxt(lib.rootdir+'/plots/paper/raw/f(E,L)_marginal.dat', usecols=[0,1,6])
  #n=sqrt(data.shape[0])
  #m=np.log10(data[:,0].reshape([n,n], order='F')/1.873)
  #c=np.log10(data[:,1].reshape([n,n], order='F')/16.3349)
  #l=data[:,2].reshape([n,n], order='F')
  #sig=P2Sig(chi2.sf(2*(l-l.ravel().min()),2))
  #plt.contour(m,c,sig, levels=[1,], colors='k')
  #h0=Ellipse((0,0),0,0,fill=False, color='k')
  #plt.plot(m.ravel()[l.argmin()], c.ravel()[l.argmin()], 'ko')
  h1=plot_siglike('f(E,L)_condition','r')
  h2=plot_siglike('RBinLog30','g')
  h22=plot_siglike('RBinLog30_L','k')
  name={'Mean':'Mean', 'AD':'ADBN'}[estimator]
  h3=plot_siglike(name+'_L','b') 
  h4=plot_siglike(name+'Iter_E','c','dashed')
  h5=plot_siglike(name+'Iter_LE','m','dashed')
  h6=plot_siglike(name+'Iter_EL','y','dashed')
  
  if estimator=='Mean':
	plt.axis([-0.3,0.3, -0.3,0.3])
  else:
	plt.axis([-0.4,0.5, -0.5, 0.4])
  plt.plot(plt.xlim(),[0,0],'k:',[0,0],plt.ylim(),'k:')
  plt.xlabel(r'$\log(M/M_{\rm true})$')
  plt.ylabel(r'$\log(c/c_{\rm true})$')
  plt.legend((h1,h2,h22,h3,h4,h5,h6),('f(E,L)|','RBin','RBin|L',estimator+'|L',estimator+'|E',estimator+'|LE',estimator+'|EL'))
  if flagsave:
	plt.savefig(lib.rootdir+'/plots/paper/MaxLikeContour'+estimator+'.eps') #rasterize=True, dpi=300
	
def plot_pot(pars, linestyle='-'):
  print pars
  r=np.logspace(0,2,20)
  Halo=NFWHalo()
  Halo.define_halo(pars)
  y=-np.array(map(Halo.pot, r))
  plt.loglog(r/Halo.halo.Rs,y, linestyle)
  
def show_legend():
    a=[x for x in gca().get_children() if isinstance(x,matplotlib.patches.Ellipse)]
    l=[getp(x,'label') for x in a]
    legend(a,l,loc=3)
    
ColorList={0:'k',4:'r',8:'g',9:'b',10:'c',11:'m',12:'y',13:'k', 14:'c'}
def image_newscan(mid, proxy, indir='scanMock10000Zoom', flagsave=False):
  if proxy=='EL':
	proxy='LE'
  name='|'.join([lib.NameList[mid],proxy]).rstrip('|')
  infile=lib.rootdir+'/plots/'+indir+'/'+name.replace('|','_')+'.hdf5'
  f=h5py.File(infile,'r')
  mm=f['/logm'][...]
  cc=f['/logc'][...]
  ts=f['/ts']
  nbin=ts.attrs['nbinE']*ts.attrs['nbinL']
  dm=mm[0,1]-mm[0,0]
  dc=cc[1,0]-cc[0,0]
  extent=[mm.ravel().min()-dm/2,mm.ravel().max()+dm/2,cc.ravel().min()-dc/2,cc.ravel().max()+dc/2]
  plt.figure()
  plt.imshow(ts, extent=extent, cmap=plt.cm.summer)
  #cs=plt.contour(mm,cc,y)
  #plt.clabel(cs, inline=1)
  levels=chi2.ppf([0.683, 0.954, 0.997], nbin)
  cs=plt.contour(mm,cc,ts, levels=levels)
  plt.clabel(cs, inline=1, fmt={levels[0]:r'1$\sigma$', levels[1]: r'2$\sigma$', levels[2]: r'3$\sigma$'})
  plt.plot(mm.ravel()[ts.value.argmin()],cc.ravel()[ts.value.argmin()],'ro')
  plt.title(name)
  plt.xlabel(r'$\log(M/M_0)$')
  plt.ylabel(r'$\log(c/c_0)$')
  if flagsave:
	plt.savefig(infile.replace('.hdf5','.eps'))

def contour_newscan(mid, proxy, indir='scanMock10000Zoom', levels={0.68:r'1$\sigma$'}, likerat=False, flagsave=False, colors=None, clabel=True, **kwargs):
  if proxy=='EL':
	proxy='LE'
  name='|'.join([lib.NameList[mid],proxy]).rstrip('|')
  infile=lib.rootdir+'/plots/'+indir+'/'+name.replace('|','_')+'.hdf5'
  f=h5py.File(infile,'r')
  mm=f['/logm'][...]
  cc=f['/logc'][...]
  ts=f['/ts']
  print ts.value.min()
  nbin=ts.attrs['nbinE']*ts.attrs['nbinL']
  fmt={}
  for k,v in levels.iteritems():
	if likerat: #likelihood ratio to determine confidence level
	  lvl=ts.value.min()+chi2.ppf(k, 2)
	else:
	  lvl=chi2.ppf(k,nbin)
	fmt[lvl]=v
  cs=plt.contour(mm,cc,ts, levels=fmt.keys(), colors=colors, **kwargs)	
  if clabel:
	plt.clabel(cs, inline=1, fmt=fmt)
  h=Ellipse((0,0),0,0,fill=False, color=list(colors)[0], label=name)
  plt.plot(mm.ravel()[ts.value.argmin()],cc.ravel()[ts.value.argmin()],'o',color=list(colors)[0])
  plt.title(name)
  plt.xlabel(r'$\log(M/M_0)$')
  plt.ylabel(r'$\log(c/c_0)$')
  if flagsave:
	plt.savefig(infile.replace('.hdf5','.eps'))
  return h

def plot_newscan(indir='scanMock10000Zoom',percent=0.68, likerat=False, savefig=False):
  nmax=5
  colors=itertools.cycle(plt.cm.jet(np.linspace(0.,1.,nmax)))
  plt.figure()
  hall=[]
  for mid,proxy in [(8,''),(8,'L'),(10,''),(10,'L')]:
	  name='|'.join([lib.NameList[mid],proxy]).rstrip('|')
	  h=contour_newscan(mid, proxy, indir, levels={percent:name}, colors=(colors.next(),), clabel=True, likerat=likerat)
	  hall.append(h)
  plt.legend(hall,[plt.getp(x,'label') for x in hall],loc=3)
  plt.plot(plt.xlim(),[0,0],'k:',[0,0],plt.ylim(),'k:')
  plt.title(indir)
  plt.xlim([-0.25,0.25])
  plt.ylim([-0.25,0.25])
  if savefig:
	if likerat:
	  plt.savefig(lib.rootdir+'/plots/'+indir+'/contours_likerat.eps')
	else:
	  plt.savefig(lib.rootdir+'/plots/'+indir+'/contours.eps')

def plot_scan(mid,T=1,sample=0):
    SigmaDef={4:1.1,8:-1.084,9:-2.3,10:-1.15,11:-1.101,12:-0.714,13:-1.15}
    data=loadtxt(rootdir+'/data/scanZoom%d/model%d.T%d'%(sample,mid,T))
    if int(sqrt(data.shape[0]))**2!=data.shape[0]:
        print data.shape
        return [],[],[]
    n=sqrt(data.shape[0])
    x=data[:,0].reshape([n,n],order='F')
    y=data[:,1].reshape([n,n],order='F')
    z=data[:,2].reshape([n,n],order='F')
    x=log10(x)
    y=log10(y)
    figure()
    #imshow(z,extent=[x.min(),x.max(),y.min(),y.max()],origin='lower')
    if mid==4:
        lvls=[nanmax(z)-SigmaDef[mid]]
    else:
        lvls=[SigmaDef[mid]]
    lvls=nanmax(z)-array([0,0.01,0.1,0.2,1,3,5,10,20])    
    #CS=contour(x,y,z,levels=lvls,origin='lower',colors=ColorList[mid],linestyles='solid')
    CS=contour(x,y,z,levels=lvls,origin='lower',linestyles='solid') 
    #CS=contour(z,levels=lvls,extent=[x.min(),x.max(),y.min(),y.max()],origin='lower',colors=ColorList[mid])
    #clabel(CS,inline=1)
    xlabel(r'$\log(M/M_0)$')
    ylabel(r'$\log(c/c_0)$')
    title(NameList[mid]+' T%d'%T)
    #axis([-1,1,-1,1])
    h=Ellipse((0,0),0,0,fill=False, color=ColorList[mid],label=NameList[mid])
    #savefig('/work/Projects/DynDistr/plots/Scan%d_Model%d_T%d.eps'%(sample,mid,T))
    return x,y,z,h

def check_fit():
	hall=[]
	fminfit={10:[0.511067, 3.418906],11:[0.698747,1.774041],12:[0.620638, 2.495314], 13:[0.528289,3.213150],4:[0.887955,1.211126],8:[0.762275,1.527824],9:[0.862825,1.347704]}
	minuitfit={10:[1.220458,0.821241], 11:[0.484838,4.199952], 12:[0.762192, 1.660001], 13:[0.530566, 3.177460], 4:[1.218284,0.729737], 8:[0.762047,1.527773], 9:[0.862815,1.347850]}
	LowPrecfit={10:[1.015049,1.018265], 11:[1.010887, 1.058054], 4:[1.010851, 0.970722], 8:[0.761972,1.527838]}
	IniOfffit={10:[1.096413,0.925138], 11:[0.908809, 1.200061], 4:[0.907599,1.199348], 8:[0.761999, 1.527716]}
	HighPrecfit={10:[0.468171, 4.413201], 11:[0.769737,1.500002], 4:[0.798530,1.443032], 8:[0.761952, 1.527847]}
	for mid in [4,8,10,11]:
		h=plot_scan(mid,T=1,sample=0)
		plot(log10(fminfit[mid][0]), log10(fminfit[mid][1]), 'o')
		plot(log10(minuitfit[mid][0]), log10(minuitfit[mid][1]), 'd', markersize=10)
		plot(log10(LowPrecfit[mid][0]), log10(LowPrecfit[mid][1]), 'x', markersize=10)
		plot(log10(HighPrecfit[mid][0]), log10(HighPrecfit[mid][1]), '^', markersize=10)
		plot(log10(IniOfffit[mid][0]), log10(IniOfffit[mid][1]), 'v', markersize=10)
		hall.append(h[-1])
	#legend(hall,[getp(x,'label') for x in hall],loc=1)   
	#savefig('/work/Projects/DynDistr/plots/SigmaScan.eps')    
  
class EnsembleFile:
  """ file of the fitted parameters for an ensemble of realizations """
  #below are class-specific variable, static class-scope variable, that are shared by all the instances
  NameList={0:'f(E,L)',4:'RBin',8:'AD',9:'Resultant',10:'Mean',11:'KS',12:'Kuiper',13:'CosMean'}
  ColorList={0:'k',4:'r',8:'g',9:'b',10:'c',11:'m',12:'y',13:'k'} 
  def __init__(self,infile):
    """infile: the file to be examined """
    #variables here are instance-specific, different for different instances
    self.file=infile
    self.rawdata=loadtxt(infile)
    self.pars=self.rawdata[:,0:2]
    #self.pars=self.rawdata[self.rawdata[:,-1]>0,0:2]
    self.id=int(os.path.basename(self.file).split('_')[0][3:])
    self.estimator=self.NameList[self.id]
    self.color=self.ColorList[self.id]
    self.stat()
  
  def stat(self):
    self.median=median(self.pars,axis=0)
    self.mean=self.pars.mean(axis=0)
    self.std=self.pars.std(axis=0)
    self.cov=cov(self.pars.T)
    self.corrcoef=corrcoef(self.pars.T)[0,1]
    self.mean_log=log10(self.pars).mean(axis=0)
    self.std_log=log10(self.pars).std(axis=0)
    self.cov_log=cov(log10(self.pars.T))
    self.corrcoef_log=corrcoef(log10(self.pars.T))[0,1]
    print "Model %d: "%self.id+self.estimator
    print "%d out of %d good fits"%(self.pars.shape[0], self.rawdata.shape[0])
    print "Median=", self.median
    print "-----------------------"
    print "m = %.3f +- %.3f"%(self.mean[0],self.std[0])
    print "c = %.3f +- %.3f"%(self.mean[1],self.std[1])
    print "corr= %.3f"%self.corrcoef
    print "-----------------------"
    print "log(m) = %.3f +- %.3f"%(self.mean_log[0],self.std_log[0])
    print "log(c) = %.3f +- %.3f"%(self.mean_log[1],self.std_log[1])
    print "corr(log)= %.3f"%self.corrcoef_log
    print "********************************************"
          
  def plot_corner(self):    
    #labels=[r"$\log (\rho_s/\rho_{s0})$",r"$\log (r_s/r_{s0})$"]
    labels=[r"$M/M_0$",r"$c/c_0$"]
    fig=triangle.corner(self.pars,labels=labels,truths=self.median,quantiles=[0.5-0.683/2,0.5+0.683/2])
    return fig
  
  def plot_cov(self,logscale=True,**kwargs):
    if logscale:
      h=plot_cov_ellipse(self.cov_log, self.mean_log, color=self.color, label=self.estimator, **kwargs)
      plot(self.mean_log[0],self.mean_log[1], 'x', color=self.color, markersize=10, **kwargs)
      axis([-1,1,-1,1])
    else:
      h=plot_cov_ellipse(self.cov, self.mean, color=self.color, label=self.estimator, **kwargs)
      plot(self.mean[0],self.mean[1], 'x', color=self.color, markersize=10, **kwargs)
      axis([0,2,0,2])
    return h      
  
  def plot_contour(self,nbin=100, percents=0.683, logscale=True, **kwargs):
    """percents can be a list, specify the contour percentile levels"""
    h,h0=percentile_contour(self.pars.T, nbin=nbin, percents=percents, color=self.color, logscale=logscale, **kwargs)
    h.set_label(self.estimator+' %d/%d'%(self.pars.shape[0],self.rawdata.shape[0]))
    #plot(self.median[0],self.median[1],'o',markersize=10,color=self.color)
    return h

def plot_dir_cov(dir='/gpfs/data/jvbq85/DynDistr/data/mockfit_mc/',logscale=True):
    hall=[]
    for f in glob.glob(dir+'/fit*.dat'):
      ff=EnsembleFile(f)
      hall.append(ff.plot_cov(logscale))
    legend(hall,[getp(x,'label') for x in hall],loc=3)
    if logscale:
      plot([-1,1],[0,0],'k:',[0,0],[-1,1],'k:')
      axis([-0.5,0.5,-0.5,0.5])
    else:
      plot([0,2],[1,1],':',[1,1],[0,2],':')
      c=array([[0.90, -1.18], [-1.18, 1.55]])
      pos=array([0.956,0.945])
      plot_cov_ellipse(c, pos,color='k')
      axis([0,2,0,2])
      

def plot_dir_contour(name='HighPrecRand',logscale=True):
      #dir='/gpfs/data/jvbq85/DynDistr/data/mockfit'+name+'_mc/'
      dir='/mnt/charon/DynDistr/data/mockfit'+name+'_mc/'
      hall=[]
      for f in glob.glob(dir+'/fit*.dat'):
	ff=EnsembleFile(f)
	hall.append(ff.plot_contour(percents=[0.683],logscale=logscale))
      legend(hall,[getp(x,'label') for x in hall],loc=1)
      plot([0.1,10],[1,1],'k:',[1,1],[0.1,10],'k:')
      axis([0.3,3,0.3,3])
      if logscale:
	xtickloc=list(arange(xlim()[0],1,0.1))
	xtickloc.extend(arange(1,xlim()[1]+0.1,0.2))
	xticks(xtickloc, map(str, xtickloc))
	ytickloc=list(arange(xlim()[0],1,0.1))
	ytickloc.extend(arange(1,xlim()[1]+0.1,0.2))
	yticks(ytickloc, map(str, ytickloc))
	
def plot_high_prec(T=1,fit='Fmin'):
  """or fit='Fmin' """
  basedir='/mnt/charon/DynDistr/data/'
  dirlist=['mockfit'+fit+'Rand_mc','mockfit'+fit+'SwapT_mc']
  #tlist={0:0,4:0,8:0,9:1,10:0,11:0,12:1,13:0}
  tlist={0:0,4:0,8:0,9:1,10:0,11:0,12:1}
  hall=[]
  for m,t in tlist.items():
    if T==2:
      t=int(not t)
    f=basedir+dirlist[t]+'/fit%d_mc.dat'%m
    ff=EnsembleFile(f)
    hall.append(ff.plot_contour(percents=0.68,logscale=True))
  legend(hall,[getp(x,'label') for x in hall],loc=1)
  plot([0.1,10],[1,1],'k:',[1,1],[0.1,10],'k:')
  #axis([0.5,2,0.5,2])
  #xtickloc=[0.5,0.6,0.7,0.8,0.9,1,1.2,1.4,1.6,1.8,2.0]
  axis([0.3,3,0.3,3])
  xtickloc=[0.3,0.5,1.0,1.5,2.0,2.5,3.0]
  xticks(xtickloc, map(str, xtickloc))
  yticks(xtickloc, map(str, xtickloc))
  xlabel(r'$M/M_0$')
  ylabel(r'$c/c_0$')
  savefig('/work/Projects/DynDistr/plots/Ensemble%s_T%d.eps'%(fit,T))

def plot_AD(mid=8,percents=0.6):
	l=['','InitOff','HighPrec','HighPrecRand','Fmin','FminRand']
	colors=['r','g','b','c','m','k']
	hall=[]
	for i,a in enumerate(l):
		f=EnsembleFile(rootdir+'/data/mockfit'+a+'_mc/fit%d_mc.dat'%mid)
		f.color=colors[i]
		hall.append(f.plot_contour(percents=percents, logscale=True))
	l[0]='LowPrec'	
	legend(hall,l, loc=1)	
	plot([0.1,10],[1,1],'k:',[1,1],[0.1,10],'k:')
	axis([0.3,3,0.3,3])
	xtickloc=[0.3,0.5,1.0,1.5,2.0,2.5,3.0]
	xticks(xtickloc, map(str, xtickloc))
	yticks(xtickloc, map(str, xtickloc))
	xlabel(r'$M/M_0$')
	ylabel(r'$c/c_0$')
	
def chi2sig(ts,dof=1):
  '''convert a chi2 TS to sigma'''
  pval=chi2.sf(ts,dof)
  sigma=norm.ppf(1.0-pval/2)
  return sigma,pval

x0=linspace(0.8,1.2,10)

def plot_sigma_prof(id=0,dim=0,x=x0):
  '''plot the significance on the subsample specified by id'''
  select_particles(id)
  #x=linspace(0.8,1.2,20)
  if dim==0:
    y=array([like2(a,1.) for a in x])
  else:
    y=array([like2(1.,a) for a in x])
  #y,py=chi2sig(-y,1)
  #print zip(x,y)
  plot(x,y)
  if dim==0:
    xlabel(r'$\log(\rho_s/\rho_{s0})$')
  else:
    xlabel(r'$\log(\r_s/\r_{s0})$')
  ylabel(r'Discrepancy Level($\sigma$)')
  #plot(xlim(),[-2.48,-2.48],'--')
  #plot(xlim(),[-2.28,-2.28],'--')
  ylim(0,3)
  #title('RADIAL_BIN_ESTIMATOR')
  #savefig('ORBITAL_ROULETTE_ESTIMATOR.eps') 
  #show()
  return y

def normfit(data,nbin=50):
  h1=hist(data,nbin,normed=1,histtype='step')
  hold('on')
  x=h1[1]
  #x=arange(3*mean(TS))
  h2,=plot(x,norm.eps(x),'r')
  par=norm.fit(data)
  h3,=plot(x,norm.eps(x,loc=par[0],scale=par[1]),'k')
  rc('text', usetex=True)
  legend((h1[2][0],h2,h3),('data','Normal(0,1)',r'Normal$(%.3f,%.3f)$'%(par[0],par[1])))
  #ylim(0,max(h1[0])*1.2)
  return par


def plot_subsamples(dataid):
  '''plot the statistics of subsamples'''
  init(dataid)
  sid=range(100)
  y=[]
  for i in sid:
    select_particles(i)
    s=like2(1.,1.)
    y.append(s)
    print i,s
  subplot(211)  
  plot(sid,y)
  subplot(212)
  normfit(y)
  return y


#figure()
#y=plot_subsamples(0)

#init(0)  
#figure()  
#xmin=[]
#for i in range(0,100,5):
  #y=plot_significance(i)
  #print i,x0[y.argmin()],y.min()
  #xmin.append(x0[y.argmin()])
  #draw()

  
"""
figure()
x0=linspace(-0.5,1,10)
for i,m in enumerate(x0):
	freeze_energy(m,0.)
	x=m+arange(-0.5,1,0.01)
	y=array([like(a,0.) for a in x])
	plot(x,y)
	hold('on')

x0=linspace(-0.01,0.01,50)
y0=zeros_like(x0)
for i,m in enumerate(x0):
	freeze_energy(m,0.)
	y0[i]=(like(m,0.)+0.5*like(m-0.01,0.)+0.5*like(m+0.01,0.))/2.
plot(x0,y0,'k-')
xlabel('log10(Rhos/Rhos0)')
ylabel('ln(Like)')
#savefig('Like_iterative_mixed.eps')
show()
"""
