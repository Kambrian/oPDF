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

class EnsembleFile(object):
  """ file of the fitted parameters for an ensemble of realizations """
  #below are class-specific variable, static class-scope variable, that are shared by all the instances
  NameList={0:'f(E,L)',4:'RBin',8:'AD',9:'Resultant',10:'Mean',11:'KS',12:'Kuiper',13:'CosMean'}
  ColorList={0:'k',4:'r',8:'g',9:'b',10:'c',11:'m',12:'y',13:'k',16:'k'} 
  def __init__(self,infile,name=None,color=None,oldformat=False):
    """infile: the file to be examined """
    #variables here are instance-specific, different for different instances
    self.file=infile
    self.rawdata=np.loadtxt(infile)
    if oldformat:
	  self.pars=self.rawdata[self.rawdata[:,-1]>0,0:2]
	  self.id=int(os.path.basename(self.file).split('_')[0][3:])
    else:
	  self.pars=self.rawdata[self.rawdata[:,-1]==0,0:2]
	  self.id=int(os.path.basename(self.file).split('.')[0].lstrip('fit').rstrip('EL'))
    #self.pars=self.rawdata[:,0:2]
    if name is None:
	  self.name=lib.NameList[self.id]
    else:
	  self.name=name
    #self.estimator=self.NameList[self.id]
    if color is None:
	  self.color=self.ColorList[self.id]
    else:
	  self.color=color
    self.stat()
  
  def stat(self):
    self.median=np.median(self.pars,axis=0)
    self.mean=self.pars.mean(axis=0)
    self.std=self.pars.std(axis=0)
    self.cov=np.cov(self.pars.T)
    self.corrcoef=np.corrcoef(self.pars.T)[0,1]
    self.mean_log=np.log10(self.pars).mean(axis=0)
    self.std_log=np.log10(self.pars).std(axis=0)
    self.cov_log=np.cov(np.log10(self.pars.T))
    self.corrcoef_log=np.corrcoef(np.log10(self.pars.T))[0,1]
    print "Model %d: "%self.id+self.name
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
      h=plot_cov_ellipse(self.cov_log, self.mean_log, color=self.color, label=self.name, **kwargs)
      plt.plot(self.mean_log[0],self.mean_log[1], 'x', color=self.color, markersize=10, **kwargs)
      plt.axis([-1,1,-1,1])
    else:
      h=plot_cov_ellipse(self.cov, self.mean, color=self.color, label=self.name, **kwargs)
      plt.plot(self.mean[0],self.mean[1], 'x', color=self.color, markersize=10, **kwargs)
      plt.axis([0,2,0,2])
    return h      
  
  def plot_contour(self,nbin=100, percents=0.683, logscale=True, fill=False, **kwargs):
    """percents can be a list, specify the contour percentile levels"""
    data=self.pars.T
    if logscale:
	  data=np.log10(data)
    X,Y,Z=density_of_points(data, nbin=nbin)
    h,h0=percentile_contour(X,Y,Z, percents=percents, colors=(self.color,), fill=fill, **kwargs)
    h.set_label(self.name)#+' %d/%d'%(self.pars.shape[0],self.rawdata.shape[0]))
    if logscale:
	  plt.plot(np.log10(self.median[0]),np.log10(self.median[1]),'o',markersize=8,color=self.color)
	  #plt.errorbar(self.mean_log[0],self.mean_log[1],xerr=self.std_log[0]/sqrt(self.pars.shape[0]), 
			#yerr=self.std_log[1]/sqrt(self.pars.shape[0]),
			#marker='o',markersize=6,mfc=self.color)
    else:
	  plt.plot(self.median[0],self.median[1],'o',markersize=8,color=self.color)
    #plt.plot(10**self.mean_log[0],10**self.mean_log[1],'o',markersize=10,color=self.color)
    return h

def fig_EnsembleContour(logscale=True, flagsave=False, init='IniRand'):
  basedir=lib.rootdir+'/data'
  nbin=80
  plt.figure()
  hall=[]
  ff=EnsembleFile(basedir+'/MinDistEnsemble'+init+'/fit8.dat', name='AD', color=(0.,0.,0.,0.3))
  ff.plot_contourf(percents=[0.68,1],nbin=nbin,logscale=logscale)
  ff=EnsembleFile(basedir+'/MinDistEnsemble'+init+'/fit10.dat', name='Mean', color=(0,0,0,0.2))
  ff.plot_contourf(percents=[0.68,1],nbin=nbin,logscale=logscale)
  ff=EnsembleFile(basedir+'/MaxLikeEnsemble/fit0_mc.dat', name='f(E,L)', color='k', oldformat=1)
  hall.append(ff.plot_contour(percents=[0.68],nbin=nbin,logscale=logscale))
  ff=EnsembleFile(basedir+'/MaxLikeEnsemble'+init+'/fit4.dat', name='RBin', color='g')
  hall.append(ff.plot_contour(percents=[0.68],nbin=nbin,logscale=logscale))
  plt.legend(hall,[plt.getp(x,'label') for x in hall],loc=1)
  if logscale:
	#plt.axis([-0.5,0.5,-0.5,0.5])
	plt.axis([-0.3,0.3,-0.3,0.3])
	plt.minorticks_on()
	plt.plot(plt.xlim(),[0,0],'k:',[0,0],plt.ylim(),'k:')
	plt.xlabel(r'$\log(M/M_{\rm true})$')
	plt.ylabel(r'$\log(c/c_{\rm true})$')
  else:
	plt.axis([0,2,0,2])
	plt.plot(plt.xlim(),[1,1],'k:',[1,1],plt.ylim(),'k:')
	plt.xlabel(r'$M/M_{\rm true}$')
	plt.ylabel(r'$c/c_{\rm true}$')
  if flagsave:
	plt.savefig(lib.rootdir+'/plots/paper/Ensemble'+init+'.pdf')# rasterize=True, dpi=300)
	
def fig_EnsembleContourExtra(logscale=False, flagsave=False, init='IniRand'):
  basedir=lib.rootdir+'/data'
  flist=[
		 #('/MinDistEnsemble'+init+'/fit8.dat','AD'), #fail
		 #('/MinDistEnsemble'+init+'/fit10.dat','Mean'), #fail
		 #('/MaxLikeEnsemble'+init+'/fit0.dat','f(E,L)|'), #perfect
		 ('/MaxLikeEnsemble'+init+'/fit4.dat','RBin'), #good
		 #('/MaxLikeEnsemble'+init+'/fit4L.dat','RBin|L'), #worse than RBin
		 #('/MaxLikeEnsemble'+init+'/fit8.dat','AD|L'), #slightly biased at 1sigma
		 #('/MaxLikeEnsemble'+init+'/fit10.dat','Mean|L'), #ok
		 #('/MaxLikeEnsemble'+init+'/fit16.dat','ADBN|L'), #fail
		 #('/IterLikeEnsemble'+init+'/fit8E.dat','AD|E'), #biased
		 #('/IterLikeEnsemble'+init+'/fit8EL.dat','AD|EL'),#biased
		 #('/IterLikeEnsemble'+init+'/fit8LE.dat','AD|LE'),#biased
		 #('/IterLikeEnsemble'+init+'/fit10E.dat','Mean|E'), #ok
		 #('/IterLikeEnsemble'+init+'/fit10EL.dat','Mean|EL'), #ok
		 #('/IterLikeEnsemble'+init+'/fit10LE.dat','Mean|LE'), #ok
		 #('/IterLikeEnsemble'+init+'/fit16E.dat','ADBN|E'),
		 #('/IterLikeEnsemble'+init+'/fit16EL.dat','ADBN|EL'),
		 #('/IterLikeEnsemble'+init+'/fit16LE.dat','ADBN|LE')
		 ]
  plt.figure()
  hall=[]
  ColorList=itertools.cycle(plt.cm.jet(np.linspace(0.,1.,len(flist)+2)))
  #ff=EnsembleFile(basedir+'/MinDistEnsemble'+init+'/fit8.dat', name='AD', color=(0.7,0.7,0.7))
  #ff.plot_contourf(percents=[0.68,1],nbin=80,logscale=logscale)
  #ff=EnsembleFile(basedir+'/MinDistEnsemble'+init+'/fit10.dat', name='Mean', color=(0.8,0.8,0.8, 0.8))
  #ff.plot_contourf(percents=[0.68,1],nbin=80,logscale=logscale)
  ff=EnsembleFile(basedir+'/MaxLikeEnsemble/fit0_mc.dat', name='f(E,L)', color=ColorList.next(), oldformat=1)
  hall.append(ff.plot_contour(percents=[0.68],nbin=80,logscale=logscale))
  for f,name in flist:
	ff=EnsembleFile(basedir+f,name, ColorList.next())
	hall.append(ff.plot_contour(percents=[0.68],nbin=80,logscale=logscale))
  plt.legend(hall,[plt.getp(x,'label') for x in hall],loc=1)
  if logscale:
	#plt.axis([-0.5,0.5,-0.5,0.5])
	plt.axis([-0.3,0.3,-0.3,0.3])
	plt.plot(plt.xlim(),[0,0],'k:',[0,0],plt.ylim(),'k:')
	plt.xlabel(r'$\log(M/M_{\rm true})$')
	plt.ylabel(r'$\log(c/c_{\rm true})$')
  else:
	plt.axis([0,2,0,2])
	plt.plot(plt.xlim(),[1,1],'k:',[1,1],plt.ylim(),'k:')
	plt.xlabel(r'$M/M_{\rm true}$')
	plt.ylabel(r'$c/c_{\rm true}$')
  if flagsave:
	plt.savefig(lib.rootdir+'/plots/paper/Ensemble'+init+'.eps')
	
def fig_DegeneracyDistributionRad(nbin=500, logscale=True, cumulative=False, flagsave=False):
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
  
  def plot_rad_distr(pars, color, label):
	Sample.freeze_energy(pars)
	Sample.like_init(pars)
	pred=Sample.predict_radial_count(nbin, logscale)
	if cumulative:
	  plt.plot(xbin, np.hstack((0.,pred)).cumsum()/Sample.nP, '-', color=color, label=label)
	  f=np.hstack((0.,pred)).cumsum()/Sample.nP
	  return f
	else:
	  plt.step(xbin,hstack((pred,nan)), where='post', color=color)
  
  fig,ax=plt.subplots(2,sharex=True, figsize=(8,8))
  plt.axes(ax[0])
  f1=plot_rad_distr([1,1], 'r', r'$M_{\rm true},\, c_{\rm true}$')
  f2=plot_rad_distr([0.6013295 ,  2.22171445], 'g', r'$0.6M_{\rm true},\, 2.2c_{\rm true}$')
  #f3=plot_theta_distr([2.16523294471, 0.476724073476], 'b', r'$2.2M_{\rm true},\, 0.5c_{\rm true}$')
  
  if cumulative:
	plt.plot(xbin, np.hstack((0.,count)).cumsum()/Sample.nP, 'k--', label='Data')
	f0=np.hstack((0.,count)).cumsum()/Sample.nP
	plt.ylabel(r'$P(<r)$')
	plt.legend(loc='upper left')
	plt.ylim([0,1])
	plt.xlim([1,300])
  else:
	plt.step(xbin,hstack((count,nan)), where='post', color='k', linestyle='--', linewidth=2)
	plt.ylabel('Counts')
	plt.legend((r'$M_{\rm true},\, c_{\rm true}$', r'$0.6 M_{\rm true},\, 2.2 c_{\rm true}$', 'Data'),loc='upper right')
	
  if logscale:
	  plt.xscale('log')
  
  plt.axes(ax[1])
  plt.plot(xbin, f1-f0, 'r-')
  plt.plot(xbin, f2-f0, 'g-')
  #plt.plot(xbin, f3-f0, 'b-')
  plt.plot(plt.xlim(), [0,0], 'k:')
  plt.ylabel(r'$\Delta P$')
  plt.xlabel(r'$r/\mathrm{kpc}$')
  if logscale:
	  plt.xscale('log')

  plt.ylim([-0.04, 0.04])
    
  fig.subplots_adjust(left=0.15, hspace=0)
  plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
  #plt.setp(ax[1].yaxis, major_locator=MaxNLocator(nbins=5, prune='upper'))
  plt.yticks([-0.04, -0.02, 0, 0.02])
  plt.minorticks_on()
  
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
  
  def plot_theta_distr(pars, color, label):
	Sample.freeze_energy(pars)
	Sample.like_init(pars)
	if cumulative:
	  t=Sample.data['theta'][...]
	  t.sort()
	  f=(np.arange(Sample.nP)+1.)/Sample.nP
	  plt.plot(t, f, '-', label=label, color=color)
	  return t.copy(),f
	else:
	  plt.hist(Sample.data['theta'], nbin, normed=True, histtype='step', color=color, label=label)
	  
  fig,ax=plt.subplots(2, sharex=True, figsize=(8,8))
  plt.axes(ax[0])
  t1,f1=plot_theta_distr([1,1], 'r', r'$M_{\rm true},\, c_{\rm true}$')
  t2,f2=plot_theta_distr([0.6013295 ,  2.22171445], 'g', r'$0.6M_{\rm true},\, 2.2c_{\rm true}$')
  #t3,f3=plot_theta_distr([2.16523294471, 0.476724073476], 'b', r'$2.2M_{\rm true},\, 0.5c_{\rm true}$')
  plt.legend(loc='upper left')
  if cumulative:
	plt.plot([0,1],[0,1],'k:')
	#plt.plot([0.5,0.5],plt.ylim(), 'k:', plt.xlim(), [0.5, 0.5], 'k:')
	plt.ylabel(r'$P(<\theta)$')
	plt.axis([0,1,0,1])
  else:
	plt.plot([0,1],[1,1],'k:')
	plt.ylabel(r'$dP/d\theta$')
  
  plt.axes(ax[1])
  
  plt.plot(t1, f1-t1, 'r-')
  plt.plot(t2, f2-t2, 'g-')
  #plt.plot(t3, f3-t3, 'b-')
  plt.plot(plt.xlim(), [0, 0], 'k:')
  plt.ylabel(r'$\Delta P$')
  plt.xlabel(r'$\theta$')
  plt.ylim([-0.04, 0.04])
    
  fig.subplots_adjust(left=0.15, hspace=0)
  plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
  #plt.setp(ax[1].yaxis, major_locator=MaxNLocator(nbins=5, prune='upper'))
  plt.yticks([-0.04, -0.02, 0, 0.02])
  plt.minorticks_on()
  
  Sample.clean()
  lib.close()
  if flagsave:
	plt.savefig(lib.rootdir+'/plots/paper/DegeneracyDistribution.eps')
  

def fig_MinDistContourRhosRs(flagsave=False):
  f2=h5py.File(lib.rootdir+'/plots/paper/raw/MeanDistRhosRs.hdf5','r')
  mm=f2['/logm'][...]
  cc=f2['/logc'][...]
  sig=f2['/ts'][...]
  x=f2['/sig_dist'].attrs['xmin']
  dm=mm[0,1]-mm[0,0]
  dc=cc[1,0]-cc[0,0]
  extent=[mm.ravel().min()-dm/2,mm.ravel().max()+dm/2,cc.ravel().min()-dc/2,cc.ravel().max()+dc/2]
  plt.imshow(sig, extent=extent, cmap=plt.cm.summer)
  #plt.contour(mm,cc, sig,  linestyles='dashed')
  plt.contour(mm,cc, mm+2*cc)
  plt.plot(plt.xlim(),[0,0],'k:',[0,0],plt.ylim(),'k:')
  plt.xlabel(r'$\log(\rho_s/\rho_{s,{\rm true}})$')
  plt.ylabel(r'$\log(r_s/r_{s,{\rm true}})$')
  if flagsave:
	plt.savefig(lib.rootdir+'/plots/paper/MinDistRhosRs.eps')

def fig_MinDistContourPotsRs(flagsave=False):
  f1=h5py.File(lib.rootdir+'/plots/paper/raw/MeanDistPotsRs.hdf5','r')
  mm=f1['/logm'][...]
  cc=f1['/logc'][...]
  ts=f1['/ts'][...]
  sig=f1['/sig_dist'][...]
  x=f1['/sig_dist'].attrs['xmin']
  dm=mm[0,1]-mm[0,0]
  dc=cc[1,0]-cc[0,0]
  extent=[mm.ravel().min()-dm/2,mm.ravel().max()+dm/2,cc.ravel().min()-dc/2,cc.ravel().max()+dc/2]
  plt.imshow(ts, extent=extent, cmap=plt.cm.summer)
  plt.contour(mm,cc, sig, [1,3], linestyles='dashed')
  h1=Ellipse((0,0),0,0,fill=False, linestyle='dashed')
  f1.close()
  
  f2=h5py.File(lib.rootdir+'/plots/paper/raw/RBinLog30PotsRs.hdf5','r')
  mm=f2['/logm'][...]
  cc=f2['/logc'][...]
  ts=f2['/ts'][...]
  sig=f2['/sig_like'][...]
  x=f2['/sig_like'].attrs['xmin']
  dm=mm[0,1]-mm[0,0]
  dc=cc[1,0]-cc[0,0]
  extent=[mm.ravel().min()-dm/2,mm.ravel().max()+dm/2,cc.ravel().min()-dc/2,cc.ravel().max()+dc/2]
  #plt.imshow(ts, extent=extent, cmap=plt.cm.winter)
  plt.contour(mm,cc, sig, [1,3], linestyles='solid')
  h2=Ellipse((0,0),0,0,fill=False, linestyle='solid')
  f2.close()
  
  plt.legend((h1,h2), ('Mean','RBin'))
  plt.plot(plt.xlim(),[0,0],'k:',[0,0],plt.ylim(),'k:')
  plt.xlabel(r'$\log(\psi_s/\psi_{s,{\rm true}})$')
  plt.ylabel(r'$\log(r_s/r_{s,{\rm true}})$')
  if flagsave:
	plt.savefig(lib.rootdir+'/plots/paper/MinDistPotsRs.eps')
	
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
  #x2=np.log10(np.array([ 0.51680255,  3.11238955])) #initial value [3,3]
  #plt.plot(x1[0],x1[1],'gx',markersize=10)
  plt.plot(x2[0],x2[1],'gx',markersize=10)
  #x3=np.log10(np.array([2.16523294471, 0.476724073476]))
  #plt.plot(x3[0],x3[1],'gx',markersize=10)
  #plt.clabel(cs, inline=1, fmt={1:'',2:'Mean'})
  #plt.plot(x[0],x[1],'g+',markersize=10)
  #plt.xlim([-0.3,0.8])
  #plt.ylim([-0.8,0.8])
  plt.minorticks_on()
  plt.plot(plt.xlim(),[0,0],'k:',[0,0],plt.ylim(),'k:')
  plt.xlabel(r'$\log(M/M_{\rm true})$')
  plt.ylabel(r'$\log(c/c_{\rm true})$')
  plt.legend((h1,h2),('AD','Mean'))
  if flagsave:
	plt.savefig(lib.rootdir+'/plots/paper/MinDistContour.eps')

def fig_MaxLikeContour(estimator='',flagsave=False):
  '''estimator in 'AD' or 'Mean' '''
  def plot_siglike(name, color, linestyle='solid', levels=[1,]):
	marker='o' #{'solid':'s','dashed':'d','dotted':'o','dashdot':'x'}[linestyle]
	f=h5py.File(lib.rootdir+'/plots/paper/raw/'+name+'.hdf5','r')
	cs=plt.contour(f['/logm'],f['/logc'],f['/sig_like'], levels=levels, colors=color, linestyles=linestyle)
	h=Ellipse((0,0),0,0,fill=False, color=color, linestyle=linestyle)
	x=f['/sig_like'].attrs['xmin']
	print 10**x
	plt.plot(x[0],x[1],marker,color=color,markersize=6)
	f.close()
	return h
  plt.figure()
   
  f=h5py.File(lib.rootdir+'/plots/paper/raw/ADDist.hdf5','r')
  plt.contourf(f['/logm'],f['/logc'],f['/sig_dist'], levels=[0,1], colors=((0.,0,0, 0.3),)) #colors='k', alpha=0.2
  #h5=Ellipse((0,0),0,0,fill=False, color='m')
  f.close()
  
  f=h5py.File(lib.rootdir+'/plots/paper/raw/MeanDist.hdf5','r')
  plt.contourf(f['/logm'],f['/logc'],f['/sig_dist'], levels=[0,1], colors=((0,0,0, 0.2),))
  #plt.contour(f['/logm'],f['/logc'],f['/sig_dist'], levels=[3], linestyles='dashed', colors=((0.2,0.2,0.2),))
  #h6=Ellipse((0,0),0,0,fill=False, color='y')
  f.close()
 

  
  data=np.loadtxt(lib.rootdir+'/plots/paper/raw/f(E,L)_marginal.dat', usecols=[0,1,6])
  n=sqrt(data.shape[0])
  m=np.log10(data[:,0].reshape([n,n], order='F')/1.873)
  c=np.log10(data[:,1].reshape([n,n], order='F')/16.3349)
  l=data[:,2]
  l[l>1e9]=inf
  l=l.reshape([n,n], order='F')
  sig=P2Sig(chi2.sf(2*(l-l.ravel().min()),2))
  plt.contour(m,c,sig, levels=[1,], colors='k',linestyles='solid')
  h0=Ellipse((0,0),0,0,fill=False, color='k',linestyle='solid')
  plt.plot(m.ravel()[l.argmin()], c.ravel()[l.argmin()], 'ko')
  #h1=plot_siglike('f(E,L)_condition','r', levels=[1,2])
  h2=plot_siglike('RBinLog30','g')
  #h22=plot_siglike('RBinLog30_L','k')
  #name={'Mean':'Mean', 'AD':'ADBN'}[estimator]
  #h3=plot_siglike(name+'_L','b') 
  #h4=plot_siglike(name+'Iter_E','c','dashed')
  #h5=plot_siglike(name+'Iter_LE','m','dashed')
  #h6=plot_siglike(name+'Iter_EL','y','dashed')
  
  if estimator in ['Mean', '']:
	plt.axis([-0.3,0.3, -0.3,0.3])
  else:
	plt.axis([-0.4,0.5, -0.5, 0.4])
  plt.plot(plt.xlim(),[0,0],'k:',[0,0],plt.ylim(),'k:')
  plt.xlabel(r'$\log(M/M_{\rm true})$')
  plt.ylabel(r'$\log(c/c_{\rm true})$')
  plt.legend((h0,h2),('f(E,L)','RBin'))
  #plt.legend((h0,h1,h2,h22,h3,h4,h5,h6),('f(E,L)','f(E,L)|','RBin','RBin|L',estimator+'|L',estimator+'|E',estimator+'|LE',estimator+'|EL'))
  plt.minorticks_on()
  if flagsave:
	plt.savefig(lib.rootdir+'/plots/paper/MaxLikeContour'+estimator+'.pdf') #rasterize=True, dpi=300

def fig_PhaseContours(flagsave=False, nbin=20, scaled=True):
  '''plot the MeanPhase as a function of Pots and Rs. Also plot the Phase of a random particle'''
  lib.open()
  npart=1000
  pid=np.random.randint(0, npart)
  x=np.logspace(-0.5,0.5,nbin)
  z0=np.empty([nbin,nbin])
  z=np.empty([nbin,nbin])
  with Tracer('Mock') as F:
	with F.copy(5000,npart) as Sample:
	  for i in xrange(nbin):
		for j in xrange(nbin):
		  z[j,i]=Sample.freeze_and_like([x[i], x[j]],14) #c: row; m: column
		  z0[j,i]=Sample.P[pid].theta
	  #x,y,Z=Sample.scan_like(14, x,x)
  lib.close()
  x=np.log10(x)

  plt.figure()
  if not scaled:
	z=z/sqrt(12.*npart)+0.5
  cs=plt.contour(x,x,z)
  plt.clabel(cs, inline=1)
  plt.xlabel(r'$\log(\psi_s/\psi_{s0})$')
  plt.ylabel(r'$\log(r_s/r_{s0})$')
  if flagsave:
	  plt.savefig(lib.rootdir+'/plots/paper/PhaseMock_Mean.eps') #rasterize=True, dpi=300
  plt.figure()
  cs2=plt.contour(x,x,z0)
  plt.clabel(cs2, inline=1)
  plt.xlabel(r'$\log(\psi_s/\psi_{s0})$')
  plt.ylabel(r'$\log(r_s/r_{s0})$')
  if flagsave:
	  plt.savefig(lib.rootdir+'/plots/paper/extra/PhaseMock_Single.eps') #rasterize=True, dpi=300
  
def fig_MockTSprof(flagsave=False, estimator=14, nbin=30):
  """TS profile for mocks"""
  lib.open()
  with Tracer('Mock') as FullSample:
	with FullSample.copy(5000, 1000) as Sample:
	  f,ax = plt.subplots(3, sharex=True, sharey=True, figsize=(8,8))
	  for i,proxy in enumerate(['r','E','L']):
		h=[]
		for pars,linestyle in [([1,1],'b-'),([ 0.51680255,  3.11238955],'r--')]:
		  ts,x=Sample.TSprof(pars, estimator, viewtypes=proxy, nbins=nbin)
		  tmp,=ax[i].plot(xrange(nbin+1), ts, linestyle)
		  #ax[i].step(xrange(nbin+1), ts, linestyle, where='post')
		  h.append(tmp)
		  ax[i].plot(plt.xlim(),[0,0],'k:')
		  ax[i].text(nbin*0.05, 2, proxy)
	  ax[1].set_ylabel(r'$\bar{\Theta}$')
	  ax[-1].set_xlabel('Bins')
	  ax[-1].legend(h, (r'$(M_{\rm true},\, c_{\rm true})$', r'$(0.6 M_{\rm true},\, 2.2 c_{\rm true})$'),loc='lower left', frameon=0, ncol=2, fontsize=18)
  f.subplots_adjust(hspace=0)
  plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
  nbins = 7 #len(ax[0].get_yticklabels())
  plt.setp([a.yaxis for a in ax], major_locator=MaxNLocator(nbins=nbins, prune='lower',symmetric=True))
  lib.close()
  if flagsave:
	plt.savefig(lib.rootdir+'/plots/paper/TSprofMockMean.eps') #rasterize=True, dpi=300

def fig_DMTSprofShow(halo, npart=10000, flagsave=False, estimator=14, nbin=30, rmin=1.01, rmax=499):
  """TS profile for DM"""
  lib.open()
  f,ax = plt.subplots(3, sharex=True, sharey=True, figsize=(8,8))
  pars=[1,1]
  h=[]
  for halotype,linestyle,realpot in [('allN', 'k-',1), ('N','r-',1),('subN','g-',1),('strmN','b-',1), ('allN', 'k--',0)]:#,('sublistN', 'c-')]:
	with Tracer(halo+halotype, DynRMAX=rmax, DynRMIN=rmin) as FullSample:
	  with FullSample.copy(0, npart) as Sample:
		if realpot:
		  lib.init_potential_spline()
		for i,proxy in enumerate(['r','E','L']):  
		  ts,x=Sample.TSprof(pars, estimator, viewtypes=proxy, nbins=nbin)
		  tmp,=ax[i].plot(xrange(nbin+1), ts, linestyle)
		  #tmp,=ax[i].plot(x, ts, linestyle)
		  #ax[i].step(xrange(nbin+1), ts, linestyle, where='post')
		  ax[i].plot(plt.xlim(),[0,0],'k:')
		  ax[i].text(nbin*0.05, 2, proxy)
		if realpot:
		  lib.free_potential_spline()
	h.append(tmp)
  ax[1].set_ylabel(r'$\bar{\Theta}$')
  ax[-1].set_xlabel('Bins')
  ax[-1].legend(h, ('All','Halo','MainSub','MainStream','All+NFWpot'),loc='lower left', frameon=0, ncol=3, fontsize=18)
  f.subplots_adjust(hspace=0)
  plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
  nbins = 7 #len(ax[0].get_yticklabels())
  plt.setp([a.yaxis for a in ax], major_locator=MaxNLocator(nbins=nbins, prune='lower',symmetric=True))
  lib.close()
  ax[0].set_title(halo)
  if flagsave:
	plt.savefig(lib.rootdir+'/plots/paper/TSprof'+halo+'MeanCmp.eps') #rasterize=True, dpi=300
	
def fig_DMTSprof(halo, npart=1e6, flagsave=False, estimator=8, nbin=30, rmin=1, rmax=500, pars=[1,1]):
  """TS profile for DM"""
  lib.open()
  f,ax = plt.subplots(3, sharex=False, sharey=False, figsize=(8,12))
  h=[]
  labels={'r':r'$\log10(r)$', 'E': r'$\log10(E)$', 'L': r'$\log10(L^2)$'}
  with Tracer(halo+'N', DynRMIN=rmin, DynRMAX=rmax) as FullSample:
	with FullSample.copy(0, int(npart)) as All:
	  with All.select(All.data['subid']<=0) as Clean:
		for Sample,realpot,linestyle,linewidth in [(All, 1, 'k-', 2), (Clean, 1, 'r--', 2), (All, 0, 'g:', 2)]:
		  if realpot:
			lib.init_potential_spline()
		  for i,proxy in enumerate(['r','E','L']):  
			ts,x=Sample.TSprof(pars, estimator, viewtypes=proxy, nbins=nbin)
			if estimator==8:
			  ts=AD2Sig(-ts)
			  tsref=1
			  #ax[i].set_yscale('log')
			else:
			  tsref=0
			tmp,=ax[i].plot(xrange(nbin+1), ts, linestyle, linewidth=linewidth)
			#tmp,=ax[i].plot(np.log10(x), ts, linestyle, linewidth=linewidth)
			#ax[i].set_xscale('log')
			#tmp,=ax[i].step(xrange(nbin+1), ts, linestyle, where='post', linewidth=linewidth)
			#tmp,=ax[i].step(np.log10(x), ts, linestyle, where='post', linewidth=linewidth)
			ax[i].plot(ax[i].get_xlim(),[tsref,tsref],'k:')
			#ax[i].text(nbin*0.05, 2, proxy)
			plt.axes(ax[i])
			plt.xticks(range(nbin)[1::3],['%.1f'%np.log10(a) for a in x][1::3])
			ax[i].set_xlabel(labels[proxy])
		  if realpot:
			lib.free_potential_spline()
		  h.append(tmp)
		if estimator==8:
		  ax[1].set_ylabel(r'AD Discrepancy [$\sigma$]')
		elif estimator==14:
		  ax[1].set_ylabel(r'$\bar{\Theta}$')
		#ax[-1].set_xlabel('Bins')
		lib.init_potential_spline()
		sigAll=All.freeze_and_like(estimator=estimator)
		sigClean=Clean.freeze_and_like(estimator=estimator)
		lib.free_potential_spline()
		sigAllNFW=All.freeze_and_like(estimator=estimator)
		if estimator==8:
		  sigAll=AD2Sig(-sigAll)
		  sigClean=AD2Sig(-sigClean)
		  sigAllNFW=AD2Sig(-sigNFW)
		l=ax[0].legend(h, ('All(%.1f)'%sigAll,'Smooth(%.1f)'%sigClean,'All+NFW(%.1f)'%sigAllNFW),loc='upper left', frameon=0, ncol=2)
		plt.setp(l.get_texts(), fontsize=18)
		#ax[0].set_xlim(0.9,2.7);ax[1].set_xlim(4,5.2);ax[2].set_xlim(6,10)
		f.subplots_adjust(hspace=0.25, bottom=0.1)
		#plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
		#nbins = 7 #len(ax[0].get_yticklabels())
		#plt.setp([a.yaxis for a in ax], major_locator=MaxNLocator(nbins=nbins, prune='lower',symmetric=True))
		lib.close()
		ax[0].set_title(halo)

  #strpot={True:'RealPot', False:'NFWPot'}
  if flagsave:
	plt.savefig(lib.rootdir+'/plots/paper/TSprof'+halo+lib.NameList[estimator]+'.%d-%d.eps'%(rmin,rmax)) #rasterize=True, dpi=300

def fig_TSProfDM_all():
  '''plot all halos in a single plot'''
  #infile=lib.rootdir+'/plots/paper/DM/TSprofRawMeanN1e+06R300.DMClean.hdf5'
  infile=lib.rootdir+'/plots/paper/star/TSprofRawMean.hdf5'
  f=h5py.File(infile,'r')
  fig=plt.figure(figsize=(8,12))
  for i,proxy in enumerate(['r','E','L2']):
	plt.subplot(3,1,i+1)
	for halo in 'ABCDE':  
	  ts=f['/'+halo+'/'+proxy+'/tsdiff'][...]
	  #ts=ts.reshape(-1,2).mean(axis=1)*sqrt(2)
	  n=len(ts)
	  plt.plot((np.arange(n)+1.)/n*100, ts)
	plt.plot(plt.xlim(),[0,0],'k--')
	plt.xlabel(proxy[0]+' Percentile')
	plt.ylim([-30,30])
	plt.plot(plt.xlim(),[-3,-3],'k:')
	plt.plot(plt.xlim(),[3,3],'k:')
  plt.subplot(312)
  plt.ylabel(r'$\bar{\Theta}$')  
  fig.subplots_adjust(hspace=0.25, bottom=0.1)
  plt.subplot(311)
  #plt.ylim([-10,10])
  plt.subplot(311)
  plt.legend(list('ABCDE'), ncol=3, loc='upper left')
  plt.savefig(infile.replace('hdf5','eps')) #rasterize=True, dpi=300
  f.close()
  
def fig_StarTSprof(halo, npart=1e5, flagsave=True, estimator=8, nbin=30, rmin=1, rmax=500):
  """TS profile for stars"""
  pars=[1,1]
  npart=int(npart)
  lib.open()
  f,ax = plt.subplots(3, sharex=False, sharey=False, figsize=(8,12))
  h=[]
  labels={'r':r'$\log10(r)$', 'E': r'$\log10(E)$', 'L': r'$\log10(L^2)$'}
  with Tracer(halo, DynRMIN=rmin, DynRMAX=rmax, DynDataFile=halo+'star.hdf5') as FullSample:
	  with FullSample.copy(0, npart) as Star:
		with FullSample.select(FullSample.data['subid']<=0) as Clean:
		  with Clean.copy(0, npart) as StarClean:
			lib.init_potential_spline()
			for Sample,weight,linestyle,linewidth in [(StarClean, 1, 'r--', 2),(StarClean, 0, 'g:', 2)]:
			  Sample.FlagUseWeight=weight
			  for i,proxy in enumerate(['r','E','L']):  
				ts,x=Sample.TSprof(pars, estimator, viewtypes=proxy, nbins=nbin)
				if estimator==8:
				  ts=AD2Sig(-ts)
				  tsref=1
				else:
				  tsref=0
				tmp,=ax[i].plot(xrange(nbin+1), ts, linestyle, linewidth=linewidth)
				ax[i].plot(ax[i].get_xlim(),[tsref,tsref],'k:')
				plt.axes(ax[i])
				plt.xticks(range(nbin)[1::3],['%.1f'%np.log10(a) for a in x][1::3])
				ax[i].set_xlabel(labels[proxy])
			  h.append(tmp)
			if estimator==8:
			  ax[1].set_ylabel(r'AD Discrepancy [$\sigma$]')
			elif estimator==14:
			  ax[1].set_ylabel(r'$\bar{\Theta}$')
			#ax[-1].set_xlabel('Bins')
			StarClean.FlagUseWeight=1
			sigClean=StarClean.freeze_and_like(estimator=estimator)
			StarClean.FlagUseWeight=0
			sigDM=StarClean.freeze_and_like(estimator=estimator)
			if estimator==8:
			  sigClean=AD2Sig(-sigClean)
			  sigDM=AD2Sig(-sigDM)
			l=ax[0].legend(h, ('Star(%.1f)'%sigClean,'SmthNoWeight(%.1f)'%sigDM),loc='upper left', frameon=0, ncol=1)
			plt.setp(l.get_texts(), fontsize=18)
			#ax[0].set_xlim(0.9,2.7);ax[1].set_xlim(4,5.2);ax[2].set_xlim(6,10)
			f.subplots_adjust(hspace=0.25, bottom=0.1)
			lib.free_potential_spline()
			lib.close()
			ax[0].set_title(halo)
  if flagsave:
	plt.savefig(lib.rootdir+'/plots/paper/TSprof'+halo+lib.NameList[estimator]+'.%d-%d.star.eps'%(rmin,rmax)) #rasterize=True, dpi=300

def fig_DensityALLvsFoF(halo='AqB4'):
  '''compare density profile of all particles vs FoF particles. NFW is a fit to all, not the FoF.'''
  npart=1000
  halo='AqB4'
  lib.open()
  FoF=Tracer(halo,DynRMIN=0,DynRMAX=500)
  All=Tracer(halo+'all',DynRMIN=0,DynRMAX=500)

  xbin=np.logspace(0, np.log10(500), 100)
  xcen=xbin/np.sqrt(xbin[1]/xbin[0])
  vol=np.diff(np.hstack([0., xbin])**3)*np.pi*4/3
  countFoF,tmp=np.histogram(FoF.data['r'], np.hstack([0., xbin]))

  countAll,tmp=np.histogram(All.data['r'], np.hstack([0., xbin]))
  plt.plot(xcen, countFoF*FoF.mP/vol,'ro')
  plt.plot(xcen, countAll*All.mP/vol,'gx')
  c0=get_config(halo+'N')
  Halo=NFWHalo()
  Halo.M0.value=float(c0['DynM0'])
  Halo.C0.value=float(c0['DynC0'])
  Halo.define_halo([1,1])
  rs=Halo.halo.Rs
  rhos=Halo.halo.Rhos
  rv=Halo.halo.Rv
  y=rhos/(xcen/rs)/(1+xcen/rs)**2
  plt.plot(xcen, y, 'k')
  plt.loglog()

  plt.figure()
  plt.plot(xcen, (countAll-countFoF)/vol,'r')
  plt.xscale('log')
  plt.xlabel('R')
  plt.ylabel(r'$\rho$')
  plt.legend(('All','FoF','NFWfit'))
  plt.savefig(lib.rootdir+'/plots/paper/extra/DensityProf'+halo+'ALLvsFoF.eps') #rasterize=True, dpi=300

def fig_FoFMap(halo='AqA4', npart=1e5, nbin=100, method='kde', percents=[0.1,0.3,0.5,0.7,0.9, 0.95]):
  npart=int(npart)
  lib.open()
  FullSampleAll=Tracer(halo,DynRMIN=1.1,DynRMAX=499)
  All=FullSampleAll.copy(0,npart)
  FoF=All.select(All.data['haloid']==0)
  Bk=All.select(All.data['haloid']!=0) #ungrouped particle dominate over other halo particles.
  #Other=All.select(All.data['haloid']>0)
  bins=np.linspace(-500,500, nbin)
  h,h0,lvls=percentile_contour(*density_of_points(FoF.data['x'][:,:2].T, bins, method), percents=percents, colors='g')
  #plt.contour(*density_of_points(FoF.data['x'][:,:2].T, bins, method), levels=lvls, colors='r')
  plt.plot(Bk.data['x'][:5000,0], Bk.data['x'][:5000,1], '.')
  FullSampleAll.clean()
  All.clean()
  FoF.clean()
  Bk.clean()
  lib.close()
  
  #plot_circle(r=float(os.environ['DynRv']), color='k', linestyle='--')
  plt.xlabel('x[kpc]')
  plt.ylabel('y[kpc]')
  plt.savefig(lib.rootdir+'/plots/paper/extra/FoFvsBkMap'+halo+'.eps')#, rasterize=True, dpi=300)
 
  
def plot_pot(pars, linestyle='-'):
  print pars
  r=np.linspace(0, 100)
  #r=np.logspace(0,2,20)
  Halo=NFWHalo()
  Halo.define_halo(pars)
  y=-np.array(map(Halo.pot, r))
  plt.semilogy(r,y, linestyle)
  
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
  
def plot_dir_cov(dir='/gpfs/data/jvbq85/DynDistr/data/mockfit_mc/',logscale=True):
    hall=[]
    for f in glob.glob(dir+'/fit*.dat'):
      ff=EnsembleFile(f,oldformat=1)
      hall.append(ff.plot_cov(logscale))
    plt.legend(hall,[plt.getp(x,'label') for x in hall],loc=3)
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
	ff=EnsembleFile(f,oldformat=1)
	hall.append(ff.plot_contour(percents=[0.683],logscale=logscale))
      plt.legend(hall,[getp(x,'label') for x in hall],loc=1)
      plt.plot([0.1,10],[1,1],'k:',[1,1],[0.1,10],'k:')
      plt.axis([0.3,3,0.3,3])
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
    ff=EnsembleFile(f,oldformat=1)
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
		f=EnsembleFile(lib.rootdir+'/data/mockfit'+a+'_mc/fit%d_mc.dat'%mid,oldformat=1)
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
