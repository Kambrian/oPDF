#scan with binned min dist estimator
from dynio import *
from myutils import *

def im_plot(mm,cc,y,name=''):
  extent=[mm.ravel().min(),mm.ravel().max(),cc.ravel().min(),cc.ravel().max()]
  plt.figure()
  plt.imshow(y, extent=extent, cmap=cm.summer)
  cs=plt.contour(mm,cc,y)
  plt.clabel(cs, inline=1)
  plt.contour(mm,cc,y, levels=[1], colors='r')
  plt.title(name)
  
lib.open()

FullSample=Tracer('Mock')
Sample=FullSample.copy(0,1000)
#Sample.radial_count(30)
nx=100
nbin=10
x=np.logspace(-0.5,0.5,nx)
mm=[m for m in x for c in x] #the first for is top layer, the second nested, so c varies first
cc=[c for m in x for c in x]
y=[lib.like_to_chi2(Sample.freeze_and_like([m,c], estimator=10),10) for m in x for c in x]
yE=[Sample.jointE_like([m,c], nbinE=nbin)/nbin for m in x for c in x]
yL=[Sample.jointLE_like([m,c], nbinL=nbin, nbinE=1)/nbin for m in x for c in x]
yLE=[Sample.jointLE_like([m,c], nbinL=nbin, nbinE=nbin)/nbin**2 for m in x for c in x]
data=np.array([mm,cc,yE,yL,yLE]).T

estimator=8
y8=[lib.like_to_chi2(Sample.freeze_and_like([m,c], estimator=estimator),estimator) for m in x for c in x]
yE8=[Sample.jointE_like([m,c], nbinE=nbin, estimator=estimator)/nbin for m in x for c in x]
yL8=[Sample.jointLE_like([m,c], nbinL=nbin, nbinE=1, estimator=estimator)/nbin for m in x for c in x]
yLE8=[Sample.jointLE_like([m,c], nbinL=nbin, nbinE=nbin, estimator=estimator)/nbin**2 for m in x for c in x]
data8=np.array([mm,cc,yE8,yL8,yLE8]).T

mm=np.log10(np.array(mm).reshape([nx,nx], order="F"))
cc=np.log10(np.array(cc).reshape([nx,nx], order="F"))
y=np.array(y).reshape([nx,nx], order="F")
yE=np.array(yE).reshape([nx,nx], order="F")
yL=np.array(yL).reshape([nx,nx], order="F")
yLE=np.array(yLE).reshape([nx,nx], order="F")
y8=np.array(y8).reshape([nx,nx], order="F")
yE8=np.array(yE8).reshape([nx,nx], order="F")
yL8=np.array(yL8).reshape([nx,nx], order="F")
yLE8=np.array(yLE8).reshape([nx,nx], order="F")

im_plot(mm,cc,y,'Mean')
im_plot(mm,cc,yLE,'Mean|LE')

extent=np.log10(x[[0,-1,0,-1]])
plt.ion()
plt.contour(mm,cc,yE, lvls=[1], origin='lower',linestyles='solid', colors='r')
plt.contour(mm,cc,yL, lvls=[1], origin='lower',linestyles='solid', colors='g')
plt.contour(mm,cc,yLE, lvls=[1], origin='lower',linestyles='solid', colors='b')
plt.contour(mm,cc,yE8, lvls=[1], origin='lower',linestyles='dashed', colors='r')
plt.contour(mm,cc,yL8, lvls=[1], origin='lower',linestyles='dashed', colors='g')
plt.contour(mm,cc,yLE8, lvls=[1], origin='lower',linestyles='dashed', colors='b')

FullSample.clean()
Sample.clean()
lib.close()