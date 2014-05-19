import sys,os
from dynio import *
os.environ['OMP_NUM_THREADS']='32'

from matplotlib.pyplot import *
from numpy import *
from scipy.stats import chi2,norm
from myutils import *
import triangle
import glob

ion()
        
def show_legend():
  a=[x for x in gca().get_children() if isinstance(x,matplotlib.patches.Ellipse)]
  l=[getp(x,'label') for x in a]
  legend(a,l,loc=3)  
  
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
    self.pars=self.rawdata[self.rawdata[:,-1]>0,0:2]
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
  h2,=plot(x,norm.pdf(x),'r')
  par=norm.fit(data)
  h3,=plot(x,norm.pdf(x,loc=par[0],scale=par[1]),'k')
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
