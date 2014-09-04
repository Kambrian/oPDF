from dynio import *
from myutils import *
import h5py,os,sys
from scipy.stats import chi2
plt.ion()

estimator=4
halo='AqD2starN'
useweight=1

npart=1000
nbin_r=20
delta_x=0.5
lib.open()

FullSample=Tracer(halo,DynRMAX=100)
outfile=lib.rootdir+'/plots/paper/'+halo+lib.NameList[estimator]+'Weight%d'%useweight

Sample=FullSample.copy(0,npart)
Sample.mP*=float(FullSample.nP)/Sample.nP
Sample.FlagUseWeight=useweight
Sample.radial_count(nbin_r)
x1=Sample.gfmin_like(estimator, [1.,1.])
x=np.logspace(log10(x1[0][0])-delta_x,log10(x1[0][0])+delta_x,20)
y=np.logspace(log10(x1[0][1])-delta_x,log10(x1[0][1])+delta_x,20)
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
  raise('estimator must be 8 or 10')

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
plt.savefig(outfile+'.eps')
f.close()
  
Sample.clean()
FullSample.clean()
lib.close()