from pylab import *
from scipy.stats import chi2
import h5py

def fit_TS(TS,npar):
  """fit TS scale"""
  return chi2.fit(TS,f0=npar,floc=0)[2]

def plot_TShist(TS,npar):
  h1=hist(TS,50,normed=1,histtype='step')
  hold('on')
  x=h1[1]
  #x=arange(3*mean(TS))
  h2,=plot(x,chi2.pdf(x,npar),'r')
  nfit,l,TSscale=chi2.fit(TS,floc=0)
  h3,=plot(x,chi2.pdf(x,nfit,scale=TSscale),'k')
  #TSscale=fit_TS(TS,npar)
  #h3,=plot(x,chi2.pdf(x,npar,scale=TSscale),'k')
  rc('text', usetex=True)
  legend((h1[2][0],h2,h3),('data',r'$\chi^2_%d$(TS)'%npar,r'$\chi^2_%d(TS/%.1f)$'%(nfit,TSscale)))
  ylim(0,max(h1[0])*1.2)
  
if __name__=="__main__":
  f=h5py.File('/work/Projects/DynDistr/data/AD-TS.mat','r')
  TS=f['/TS'][:].T
  #scale=fit_TS(TS,1)
  plot_TShist(TS,4)
  show()