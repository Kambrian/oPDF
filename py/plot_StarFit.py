import matplotlib
#matplotlib.use('Agg')
from dynio import *
from myutils import *
import h5py,os,sys
from scipy.stats import chi2
plt.ion()

dm=np.loadtxt('/work/Projects/DynDistr/plots/paper/extra/RBin50/fit.dat', delimiter=',');
TMP=np.loadtxt('/work/Projects/DynDistr/plots/paper/extra/Star/Weight0/Clean/fitTMP.dat', delimiter=',');
NFW=np.loadtxt('/work/Projects/DynDistr/plots/paper/extra/Star/Weight0/Clean/fitNFW.dat', delimiter=',');
dm=dm.reshape(5,-1,2)
TMP=TMP.reshape(5,2)
NFW=NFW.reshape(5,2)
refM=np.array([1.842, 0.819, 1.774, 1.774, 1.185])
WentingM=np.array([[1.302, 1.970, 1.616, 1.910, 0.744],
				   [1.150, 0.867, 1.411, 1.406, 0.995]]);
refC=np.array([16.10, 8.16, 12.34, 8.73, 8.67])
WentingC=np.array([[8.616, 3.082, 7.682, 6.722, 8.751],
				   [15.269,  8.186, 15.876, 10.339, 10.297]]);
WentingM/=refM;
WentingC/=refC;
Wenting=np.vstack([WentingM.T.flat,WentingC.T.flat]).T.reshape([5,-1,2])
#data=np.concatenate((TMP, NFW, Wenting[:,[1],:], dm[:,[3],:]), 1);
def plot_fit(mc, color, fill, alpha=1):
  shapes='osd^*'
  facecolor='none'
  if fill:
	facecolor=color
  if mc.shape[0]==2&mc.shape[1]!=2:#change to [n,2]
	mc=mc.T
  h=[]
  for i,x in enumerate(mc):
	htmp,=plt.plot(x[0],x[1], shapes[i], mec=color, mfc=facecolor, markersize=12)
	h.append(htmp)
  cov=np.cov(mc.T)
  print mc.mean(0), mc.std(0,ddof=1), np.sqrt(((mc-1)**2).mean(0))
  #plot_cov_ellipse(cov, [1,1], color=color)
  h2=plot_cov_ellipse(cov, mc.mean(0), color=color, fill=1, alpha=alpha)
  return h,h2

plt.figure()
h,h1=plot_fit(NFW, 'r', 1, 0.3)
_,h2=plot_fit(TMP, 'g', 1, 0.4)
_,h3=plot_fit(Wenting[:,1,:], 'b', 1, 0.5)
_,h4=plot_fit(dm[:,3,:],'grey', 1, 0.7)
plt.axis([0.5,1.5,0.5,1.5]);
plt.plot(plt.xlim(),[1,1],'k--',[1,1],plt.ylim(),'k--');
plt.xlabel(r'$M/M_0$');
plt.ylabel(r'$c/c_0$');
legend1=plt.legend(h, ('A','B','C','D','E'),loc=1,fontsize=15)
plt.legend((h1,h2,h3,h4),('NFW','TMP','f(E,L)','DM-TMP'), loc=3,fontsize=15)
plt.gca().add_artist(legend1)
plt.savefig('/work/Projects/DynDistr/plots/paper/extra/StarFit.Filled.pdf');

