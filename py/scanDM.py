import matplotlib
matplotlib.use('Agg')
from dynio import *
from myutils import *
import h5py,os,sys
from scipy.stats import chi2
#plt.ion()

estimator=4
halo=sys.argv[1] #'A2'
cleansample=int(sys.argv[2]) #0
usetemplate=int(sys.argv[3]) #1
scancontour=int(sys.argv[4]) #1
rmin=float(sys.argv[5])
rmax=float(sys.argv[6])
#rmin=1;rmax=300
#rmin=10;rmax=100 

npart=1e6#1000
nbin_r=100
delta_x=0.05/sqrt(npart/1e5) #sigma~0.02 for 1e5 particles. for 1e6 particles, one can get 1 percent accuracy in M,c
ngrid=10
lib.open()

outfile=lib.rootdir+'/plots/paper/extra/RBin100/'+halo+lib.NameList[estimator]+'Template%dClean%d.R%d-%dN%g'%(usetemplate,cleansample,rmin,rmax,npart)
sys.stdout = open(outfile+'.log', 'w')

FullSample=Tracer(halo,DynRMIN=rmin,DynRMAX=rmax)

Sample=FullSample.copy(0,npart)
print np.sum(Sample.data['subid']<=0)
if cleansample:
  newsample=Sample.select(Sample.data['subid']<=0)
  Sample.clean()
  Sample=newsample
print Sample.nP
Sample.mP*=float(FullSample.nP)/Sample.nP
#Sample.FlagUseWeight=useweight
Sample.radial_count(nbin_r)

if usetemplate:
  lib.init_potential_spline()

x1=Sample.gfmin_like(estimator, [1.,1.])
print 'x1', x1
#with open(outfile+'.log','w') as logfile:
  #logfile.write(str(x1))

if scancontour:
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