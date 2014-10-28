import matplotlib
#matplotlib.use('Agg')
from dynio import *
from myutils import *
import h5py,os,sys
from scipy.stats import chi2
plt.ion()

dm=np.loadtxt('/work/Projects/DynDistr/plots/paper/extra/RBin50/fit.dat', delimiter=',').reshape(5,-1,2);
TMP=np.loadtxt('/work/Projects/DynDistr/plots/paper/extra/Star/Weight0/Clean/fitTMP.dat', delimiter=',').reshape(5,2);
NFW=np.loadtxt('/work/Projects/DynDistr/plots/paper/extra/Star/Weight0/Clean/fitNFW.dat', delimiter=',').reshape(5,2);
starFromDM=np.loadtxt('/work/Projects/DynDistr/plots/paper/extra/Star/FromDM/Weight0/Clean/fit.dat', delimiter=',').reshape(5,2);
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
  h2=plot_cov_ellipse(cov, mc.mean(0), color=color, fill=0, alpha=alpha)
  return h,h2

plt.figure()
h,h1=plot_fit(TMP, 'r', 1, 1)
_,h3=plot_fit(dm[:,3,:],'grey', 1, 1)
_,h2=plot_fit(starFromDM, 'g', 1, 1)
plt.axis([0.5,1.5,0.5,1.5]);
plt.plot(plt.xlim(),[1,1],'k--',[1,1],plt.ylim(),'k--');
plt.xlabel(r'$M/M_0$');
plt.ylabel(r'$c/c_0$');
legend1=plt.legend(h, ('A','B','C','D','E'),loc=1,fontsize=15)
plt.legend((h1,h2,h3),('Star','dmStar','DM'), loc=3,fontsize=15)
plt.gca().add_artist(legend1)
#plt.savefig('/work/Projects/DynDistr/plots/paper/StarFit.FromDM.eps');

