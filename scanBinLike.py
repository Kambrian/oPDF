#scan with binned min dist estimator
from dynio import *
import matplotlib.pyplot as pl

def binlike(pars, Sample, estimator=14, proxy='L2', nbin=10):
  '''slice and combine the likelihood'''
  ADpar=[-0.22, 0.66] 
  Sample.freeze_energy(pars)
  x=Sample.gen_bin(proxy, nbin, equalcount=True)[0]
  neglnL=0
  for i in range(nbin):
	with Sample.copy() as sub:
	  sub.data['flag'][:]=(sub.data[proxy]>x[i])*(sub.data[proxy]<x[i+1])
	  sub.squeeze()
	  if proxy=='r':
		sub.rmin,sub.rmax=x[i],x[i+1]
	  if estimator==4:
		sub.radial_count()
		ts=-sub.likelihood(pars, estimator)
	  elif estimator==8:
		ts=sub.likelihood(pars, estimator)
		ts=(np.log(-ts)-ADpar[0])/ADpar[1]  #from lognormal to normal
		ts*=ts/2
	  elif estimator==14:
		ts=sub.likelihood(pars, estimator)
		ts*=ts/2
	  neglnL+=ts 
  return neglnL #negative loglikelihood	
  
FullSample=Tracer('Mock')
Sample=FullSample.copy(0,10000)
Sample.radial_count(30)
lib.open()
x=np.linspace(0.8,1.2,20)
#yL=[binlike([m,1],Sample, proxy='L2') for m in x]
y=np.array([lib.like_to_chi2(Sample.freeze_and_like([m,1],10),10) for m in x])
yL=np.array([Sample.jointLE_like([m,1], nbinL=10, nbinE=1) for m in x])
yE=np.array([Sample.jointE_like([m,1], nbinE=10) for m in x])
yLE=np.array([Sample.jointLE_like([m,1], nbinL=10, nbinE=10) for m in x])

y8=np.array([lib.like_to_chi2(Sample.freeze_and_like([m,1],8),8) for m in x])
yL8=np.array([Sample.jointLE_like([m,1], nbinL=10, nbinE=1, estimator=8) for m in x])
yE8=np.array([Sample.jointE_like([m,1], nbinE=10, estimator=8) for m in x])
yLE8=np.array([Sample.jointLE_like([m,1], nbinL=10, nbinE=10, estimator=8) for m in x])

pl.ion()
y0=0 #ylim[0]
y00=10#ylim[1]
pl.plot(x,y, 'r-')
pl.plot(x,yL/10,'g-')
pl.plot(x,yE/10,'b-')
pl.plot(x,yLE/100,'k-')
pl.plot(x,y8, 'r--')
pl.plot(x,yL8/10,'g--')
pl.plot(x,yE8/10,'b--')
pl.plot(x,yLE8/100,'k--')
l=pl.legend(('Mean','Mean|L','Mean|E','Mean|LE','AD','AD|L','AD|E','AD|LE'),loc='lower left',frameon=False, fontsize=10)
pl.xlabel(r'$M/M_0$')
pl.ylabel(r'Reduced $\chi^2$')
pl.yscale('log')
pl.savefig(lib.rootdir+'/plots/TSchi2_EL.eps')

if __name__=='__main__':
  #---------------------------mass plot------------------------------------
  yr=[binlike([m,1],Sample, proxy='r') for m in x]
  yE=[binlike([m,1],Sample, proxy='E') for m in x]
  yL8=[binlike([m,1],Sample, proxy='L2', estimator=8) for m in x]
  yr8=[binlike([m,1],Sample, proxy='r', estimator=8) for m in x]
  yE8=[binlike([m,1],Sample, proxy='E', estimator=8) for m in x]
  yL4=[binlike([m,1],Sample, proxy='L2', estimator=4) for m in x]
  yr4=[binlike([m,1],Sample, proxy='r', estimator=4) for m in x]
  yE4=[binlike([m,1],Sample, proxy='E', estimator=4) for m in x] #this simply doesn't work, if you allow the energy to freely shift.
  y4=[-Sample.freeze_and_like([m,1], 4) for m in x]
  y8=[-Sample.freeze_and_like([m,1],8) for m in x]
  y10=[-Sample.freeze_and_like([m,1],10) for m in x]

  pl.ion()
  y0=0 #ylim[0]
  y00=10#ylim[1]
  pl.plot(x,yL-np.min(yL)+y0,'r-')
  pl.plot(x,yr-np.min(yr)+y0,'g-')
  pl.plot(x,yE-np.min(yE)+y0,'b-')
  pl.plot(x, y4-np.min(y4)+y0, 'c--')
  pl.plot(x, y8-np.min(y8)+y0, 'm--')
  pl.plot(x, y10-np.min(y10)+y0, 'y--')
  pl.legend(('Mean|L','Mean|r','Mean|E','RBin','AD','Mean'))
  pl.xlabel(r'$M/M_0$')
  pl.ylabel('TS')
  #pl.yscale('log')
  pl.ylim([y0,y00])
  pl.savefig(lib.rootdir+'/plots/TSchi2.eps')

  pl.figure()
  pl.plot(x,yL-np.min(yL)+y0,'r-')
  pl.plot(x,yr-np.min(yr)+y0,'g-')
  pl.plot(x,yE-np.min(yE)+y0,'b-')
  pl.plot(x,yL8-np.min(yL8)+y0,'r--')
  pl.plot(x,yr8-np.min(yr8)+y0,'g--')
  pl.plot(x,yE8-np.min(yE8)+y0,'b--')
  pl.plot(x,yL4-np.min(yL4)+y0,'r:')
  pl.plot(x,yr4-np.min(yr4)+y0,'g:')
  pl.plot(x,yE4-np.min(yE4)+y0,'b:')
  pl.legend(('Mean|L','Mean|r','Mean|E','AD|L','AD|r','AD|E','RBin|L','RBin|r','RBin|E'))
  pl.xlabel(r'$M/M_0$')
  pl.ylabel('TS')
  #pl.yscale('log')
  pl.ylim([y0,y00])
  pl.savefig(lib.rootdir+'/plots/TSchi2_bin.eps')
  #-------------------------c plot----------------
  yL=[binlike([1,c],Sample, proxy='L2') for c in x]
  yr=[binlike([1,c],Sample, proxy='r') for c in x]
  yE=[binlike([1,c],Sample, proxy='E') for c in x]
  yL8=[binlike([1,c],Sample, proxy='L2', estimator=8) for c in x]
  yr8=[binlike([1,c],Sample, proxy='r', estimator=8) for c in x]
  yE8=[binlike([1,c],Sample, proxy='E', estimator=8) for c in x]
  y4=[-Sample.freeze_and_like([1,c], 4) for c in x]
  y8=[-Sample.freeze_and_like([1,c],8) for c in x]
  y10=[-Sample.freeze_and_like([1,c],10) for c in x]

  pl.figure()
  pl.ion()
  pl.plot(x,yL-np.min(yL)+y0,'r-')
  pl.plot(x,yr-np.min(yr)+y0,'g-')
  pl.plot(x,yE-np.min(yE)+y0,'b-')
  pl.plot(x, y4-np.min(y4)+y0, 'c--')
  pl.plot(x, y8-np.min(y8)+y0, 'm--')
  pl.plot(x, y10-np.min(y10)+y0, 'y--')
  pl.legend(('Mean|L','Mean|r','Mean|E','RBin','AD','Mean'))
  pl.xlabel(r'$c/c_0$')
  pl.ylabel('TS')
  #pl.yscale('log')
  pl.ylim([y0,y00])
  pl.savefig(lib.rootdir+'/plots/TSchi2_concentration.eps')

  pl.figure()
  pl.plot(x,yL-np.min(yL)+y0,'r-')
  pl.plot(x,yr-np.min(yr)+y0,'g-')
  pl.plot(x,yE-np.min(yE)+y0,'b-')
  pl.plot(x,yL8-np.min(yL8)+y0,'r--')
  pl.plot(x,yr8-np.min(yr8)+y0,'g--')
  pl.plot(x,yE8-np.min(yE8)+y0,'b--')
  pl.legend(('Mean|L','Mean|r','Mean|E','AD|L','AD|r','AD|E'))
  pl.xlabel(r'$c/c_0$')
  pl.ylabel('TS')
  #pl.yscale('log')
  pl.ylim([y0,y00])
  pl.savefig(lib.rootdir+'/plots/TSchi2_AD_concentration.eps')
  #from scipy.optimize import fmin
  #blike=lambda pars: binlike(pars, Sample)
  #x=fmin(blike, [1,1])
  #like10=lambda x: -Sample.freeze_and_like(x, 10)
  #x10=fmin(like10, [1,1])
  #like8=lambda x: -Sample.freeze_and_like(x, 8)
  #x8=fmin(like8,[1,1])
  #print x,x8,x10
  #print blike(x),blike(x8),blike(x10)

  #lib.free_integration_space()
  #with Tracer('Mock') as FullSample:
	  #with FullSample.sample(0) as Sample:
		#binlike([1,1],Sample)
		
  #lib.close()	  