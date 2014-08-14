#using nested_views rather than jointEL
from dynio import *
from myutils import *
import h5py,os,sys
from scipy.stats import chi2
plt.ion()

#estimator=8
#proxy='LE'
#nbins=[10,10]
npart=1000

#lib.open()
#FullSample=Tracer('Mock')
#Sample=FullSample.copy(0,npart)
#FullSample.clean()

def fig_MockTSprof(halo, npart=1000, flagsave=True, estimator=14, nbin=30):
  """TS profile for mocks"""
  lib.open()
  with Tracer(halo) as FullSample:
	with FullSample.copy(0, npart) as Sample:
	  #Sample.FlagUseWeight=useweight
	  f,ax = plt.subplots(3, sharex=True, sharey=True, figsize=(8,8))
	  for i,proxy in enumerate(['r','E','L']):
		h=[]
		pars=[1,1]
		for Sample.FlagUseWeight,linestyle in [(0,'b-'), (1,'r--')]:
		  ts,x=Sample.TSprof(pars, estimator, viewtypes=proxy, nbins=nbin)
		  tmp,=ax[i].plot(xrange(nbin+1), ts, linestyle)
		  #ax[i].step(xrange(nbin+1), ts, linestyle, where='post')
		  h.append(tmp)
		  ax[i].plot(plt.xlim(),[0,0],'k:')
		  ax[i].text(nbin*0.05, 2, proxy)
	  ax[1].set_ylabel(r'$\bar{\Theta}$')
	  ax[-1].set_xlabel('Bins')
	  ax[-1].legend(h, ('No Weight','Weight'),loc='lower left', frameon=0, ncol=2, fontsize=18)
  f.subplots_adjust(hspace=0)
  plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
  nbins = 7 #len(ax[0].get_yticklabels())
  plt.setp([a.yaxis for a in ax], major_locator=MaxNLocator(nbins=nbins, prune='lower',symmetric=True))
  lib.close()
  ax[0].set_title(halo)
  if flagsave:
	plt.savefig(lib.rootdir+'/plots/paper/TSprof'+halo+'Mean.eps') #rasterize=True, dpi=300
	


#fig_MockTSprof('AqA2star')
fig_MockTSprof('AqA2starN')
#fig_MockTSprof('AqA2starN', useweight=1)
#fig_MockTSprof('AqB2star')
fig_MockTSprof('AqB2starN')
#fig_MockTSprof('AqB2starN', useweight=1)
  
#Sample.clean()
#lib.close()