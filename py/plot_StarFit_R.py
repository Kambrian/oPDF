#TODO: plot in pots rs space.
import matplotlib
#matplotlib.use('Agg')
from dynio import *
from myutils import *
import h5py,os,sys
from scipy.stats import chi2
plt.ion()

#dm=np.loadtxt('/work/Projects/DynDistr/plots/paper/extra/RBin50/fit.dat', delimiter=',');
#dm=dm.reshape(5,-1,2)
TMP=np.loadtxt('/work/Projects/DynDistr/plots/paper/extra/Star/Weight0/Clean/fitTMP.dat.PsiR', delimiter=',');
TMP=TMP.reshape(5,2)
#NFW=np.loadtxt('/work/Projects/DynDistr/plots/paper/extra/Star/Weight0/Clean/fitNFW.dat', delimiter=',');
#NFW=NFW.reshape(5,2)
#DUP=np.loadtxt('/work/Projects/DynDistr/plots/paper/extra/Star/Dup/Weight0/Clean/fit.dat', delimiter=',');
#DUP=DUP.reshape(5,2)

R20=np.loadtxt('/work/Projects/DynDistr/plots/paper/extra/Star/Weight0/Clean/r20/fitTMP.dat.PsiR', delimiter=',');
R20=R20.reshape(5,2)
R40=np.loadtxt('/work/Projects/DynDistr/plots/paper/extra/Star/Weight0/Clean/r40/fitTMP.dat.PsiR', delimiter=',');
R40=R40.reshape(5,2)
R60=np.loadtxt('/work/Projects/DynDistr/plots/paper/extra/Star/Weight0/Clean/r60/fitTMP.dat.PsiR', delimiter=',');
R60=R60.reshape(5,2)

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
  print mc.mean(0), mc.std(0,ddof=1)
  #plot_cov_ellipse(cov, [1,1], color=color)
  h2=h[0]
  h2=plot_cov_ellipse(cov, mc.mean(0), color=color, fill=0, alpha=alpha)
  return h,h2

plt.figure()
h,h1=plot_fit(TMP, 'grey', 1, 1)
_,h2=plot_fit(R20, 'r', 1, 1)
_,h3=plot_fit(R40, 'g', 1, 1)
_,h4=plot_fit(R60, 'k', 1, 1)
#plt.axis([0.5,1.5,0.5,1.5]);
plt.plot(plt.xlim(),[1,1],'k--',[1,1],plt.ylim(),'k--');
#plt.xlabel(r'$M/M_0$');
#plt.ylabel(r'$c/c_0$');
plt.xlabel(r'$\psi_s/\psi_{s0}$');
plt.ylabel(r'$r_s/r_{s0}$');
legend1=plt.legend(h, ('A','B','C','D','E'),loc=2,fontsize=15)
plt.legend((h1,h2,h3,h4),('Full',r'$r>20\%$',r'$r>40\%$',r'$r>60\%$'), loc=4,fontsize=15)
plt.gca().add_artist(legend1)
plt.savefig('/work/Projects/DynDistr/plots/paper/extra/StarFit.R.PotsRs.eps');

