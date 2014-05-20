import sys,os,ConfigParser
from dynio import *
#os.environ['OMP_NUM_THREADS']='32'
from matplotlib.pyplot import *
from numpy import *
from scipy.stats import chi2,norm
from myutils import *
import triangle
import glob

def pick_particles(sample=0):
  select_particles(sample)
  P=ParticleData()
  freeze_energy(1,1)
  return P

def TSprof(bintype='E',sample=0,nbin=20,estimator=8):
  P=pick_particles(sample)  
  rmin=P.R_MIN.value
  rmax=P.R_MAX.value
  proxy=copy(P.P[bintype])
  if bintype=='E':
    proxy=-proxy
  xmin,xmax=proxy[proxy>0].min(),proxy.max()
  x=logspace(log10(xmin),log10(xmax),nbin)
  TS=[]
  for i,x0 in enumerate(x[:-1]):
    x1=x[i+1]
    P.P['flag']=0
    P.P['flag'][logical_and(proxy>=x0,proxy<x1)]=1
    if bintype=='r': #update radial limits for radial cut
      P.R_MIN.value=x0
      P.R_MAX.value=x1
    squeeze_data() #radial limits actually get updated inside this func.
    if bintype=='r': #update the bin counts for radial cut
      lib.fill_radial_bin()
    TS.append(like(1,1,estimator))
    print i, P.nP.value, TS[-1]
    P.R_MIN.value=rmin #restore radial limits
    P.R_MAX.value=rmax
    P=pick_particles(sample) 
  
  #figure()
  #subplot(211)
  #hist(proxy,x,log=True)
  #xscale('log')
  #xlim(x.min(),x.max())
  #ylabel('Count')
  #subplot(212)
  plot(x[:-1],AD2Sig(-array(TS)))
  xscale('log')
  xlim(x.min(),x.max())
  xlabel(bintype)
  ylabel(r'Discrepancy/$\sigma$')
  return x,TS

def plot_halo_TS(halo,estimator=8):    
  get_config(halo)
  os.environ['DynRMAX']='500'
  init()
  subplot(311)
  TSprof(bintype='r',estimator=estimator)
  title(" ".join([halo,NameList[estimator]]))
  subplot(312)
  TSprof(bintype='E',estimator=estimator)
  subplot(313)
  TSprof(bintype='L2',estimator=estimator)
  
if __name__=="__main__":
  #ion()
  halo='AqA4'
  estimator=8
  if len(sys.argv)>1:
    halo=sys.argv[1]  
  if len(sys.argv)>2:
    estimator=int(sys.argv[2])
  plot_halo_TS(halo,estimator)
  show()

