import matplotlib
matplotlib.use('Agg')
from dynio import *
from myutils import *
import h5py,os,sys
from scipy.stats import chi2
#plt.ion()

estimator=4
halo=sys.argv[1] #'A2'
cleansample=1
usetemplate=1
useweight=0
scancontour=0
proxycut=None
rmin=10;rmax=300

npart=int(1e6)
nbin_r=50
ngrid=20
lib.open()

outdir=lib.rootdir+'/plots/paper/extra/Star/Dup/Weight%d/'%useweight
if cleansample==0:
  outdir+='Full/'
elif cleansample==1:
  outdir+='Clean/'
elif cleansample==2:
  outdir+='BranchClean/'
else:
  print "error: cleansample unknown"
  raise 
if proxycut!=None:
  outdir=outdir+proxycut+'%02d/'%(percentcut*100)
if not os.path.exists(outdir):
  try:
	os.makedirs(outdir)  
  except:
	pass
outfile=outdir+halo+lib.NameList[estimator]+'Template%d.R%d-%dN%g'%(usetemplate,rmin,rmax,npart)
sys.stdout = open(outfile+'.log', 'w')

FullSample=Tracer(halo,DynRMIN=rmin,DynRMAX=rmax,DynDataFile='old/'+halo+'star.hdf5')

if cleansample:
  if cleansample==1:
	newsample=FullSample.select(FullSample.data['subid']<=0)
  if cleansample==2:
	newsample=FullSample.select(FullSample.data['haloid']==0)
  FullSample.clean()
  FullSample=newsample

#cut in r or L percentiles
if proxycut!=None:
  if proxycut=='L':
	proxycut='L2'
  x=np.percentile(FullSample.data[proxycut], percentcut*100)
  if proxycut=='E':
	if usetemplate:
	  lib.init_potential_spline()
	FullSample.freeze_energy()
	x=np.percentile(FullSample.data[proxycut], percentcut*100)
	newsample=FullSample.select(FullSample.data[proxycut]>x)
	if usetemplate:
	  lib.free_potential_spline()
  else:
	newsample=FullSample.select(FullSample.data[proxycut]<x)
  FullSample.clean()
  FullSample=newsample
  if proxycut=='r':
	FullSample.rmax=x
	print 'Rcut=', x

if FullSample.nP<npart:
  print "Warning: full sample size=%d, required=%d"%(FullSample.nP,npart)
  npart=0
  
Sample=FullSample.copy(0,npart)
#Nunique=np.sum(Sample.data['haloid']>0)
#if useweight:
  #Sample.data['w']=Sample.data['w']*float(Nunique)/Sample.nP #rescale the sample size to correct for duplication
#print Sample.nP, Nunique
Sample.FlagUseWeight=useweight
Sample.radial_count(nbin_r)

if usetemplate:
  lib.init_potential_spline()

x1=Sample.gfmin_like(estimator, [1.,1.])
print 'x1', x1
#with open(outfile+'.log','w') as logfile:
  #logfile.write(str(x1))

if scancontour:
  delta_x=0.05/sqrt(Sample.nP/1e5) #sigma~0.02 for 1e5 particles. for 1e6 particles, one can get 1 percent accuracy in M,c
  x=np.logspace(log10(x1[0][0])-delta_x,log10(x1[0][0])+delta_x,ngrid)
  y=np.logspace(log10(x1[0][1])-delta_x,log10(x1[0][1])+delta_x,ngrid)
  #x=np.logspace(-delta_x,delta_x,20)
  #y=np.logspace(-delta_x,delta_x,20)
  cont1=Sample.scan_Chi2(estimator, x,y)
  cont1=cont1[0],cont1[1],cont1[2]-x1[1]
  if estimator==8:
	sig=AD2Sig(cont1[2])
  elif estimator==10:
	sig=P2Sig(chi2.sf(cont1[2],1))
  elif estimator==4:
	sig=P2Sig(chi2.sf(cont1[2],2))
  else:
	raise('estimator must be 4, 8 or 10')

  f=h5py.File(outfile+'.hdf5','w')
  grp=f.create_group('/Dynamical')
  grp.create_dataset('logm', data=np.log10(x))
  grp.create_dataset('logc', data=np.log10(y))
  dset=grp.create_dataset('sig_like', data=sig)
  dset.attrs['xmin']=np.log10(x1[0])
  dset=grp.create_dataset('ts', data=cont1[2])

  f.close()

  f=h5py.File(outfile+'.hdf5','r')
  grp=f['/Dynamical']
  x2=grp['sig_like'].attrs['xmin']
  h2,=plt.plot(x2[0],x2[1],'go', markersize=10)
  plt.contour(grp['logm'],grp['logc'],grp['sig_like'],levels=[1,2,3], colors='g', linestyles='solid')
  #h1=contour_handle('r', 'dashed')
  #h2=contour_handle('g', 'solid')
  #plt.legend((h1,h2), ('Density','Dynamics'))
  plt.xlabel(r'$\log(M/M_0)$')
  plt.ylabel(r'$\log(c/c_0)$')
  #plt.axis([x2[0]-delta_x,log10(x2[0])+delta_x,x2[1]-delta_x,x2[1]+delta_x])
  plt.plot(plt.xlim(),[0,0],'k:',[0,0],plt.ylim(),'k:')
  plt.subplots_adjust(left=0.15)
  plt.savefig(outfile+'.eps')
  f.close()

if usetemplate:
  lib.free_potential_spline()  
Sample.clean()
FullSample.clean()
lib.close()