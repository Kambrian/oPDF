#using nested_views rather than jointEL
from dynio import *
from myutils import *
import h5py,os,sys
from scipy.stats import chi2
import itertools
plt.ion()

#estimator=8
#proxy='LE'
#nbins=[10,10]
npart=10000

halo='AqA4'
lib.open()

c0=get_config(halo)
Halo=NFWHalo()
Halo.M0.value=float(c0['DynM0'])
Halo.C0.value=float(c0['DynC0'])
Halo.define_halo([1,1])
rv0=Halo.halo.Rv

c1=get_config(halo+'N')
Halo=NFWHalo()
Halo.M0.value=float(c1['DynM0'])
Halo.C0.value=float(c1['DynC0'])
Halo.define_halo([1,1])
rv1=Halo.halo.Rv
m1,c1=float(c1['DynM0'])/float(c0['DynM0']),float(c1['DynC0'])/float(c0['DynC0'])

rmax=[30, 50,100,200,300,rv0,rv1]
ColorList=itertools.cycle(plt.cm.jet(np.linspace(0.,1.,len(rmax)+1)))
hall=[]
mlist=[]
for DynRMAX in rmax:
  with Tracer(halo+'all',DynRMIN=1, DynRMAX=DynRMAX) as FullSample:
	with FullSample.copy(0,npart) as Sample:
	  Sample.mP*=float(FullSample.nP)/Sample.nP
	  #Sample.radial_count()
	  result,m=Sample.minuit_NFWlike()
	  mlist.append((m.values['m'],m.errors['m']))
	  cont=m.contour('m','c',subtract_min=True)
	  color=ColorList.next()
	  plt.contour(cont[0],cont[1],cont[2], levels=[1.15,], colors=[color,])
	  htmp=Ellipse((0,0),0,0,fill=False, color=color)
	  hall.append(htmp)
	
h1,=plt.plot(1,1,'o')
h2,=plt.plot(m1, c1, 'd')
hall.extend([h1,h2])
l=plt.legend(hall,('30', '50','100','200','300','RvS','RvN','FitS','FitN'), loc='lower left', ncol=3, fontsize=15)
plt.xlabel(r'$M/M_0$')
plt.ylabel(r'$c/c_0$')
plt.axis([0.6,1.2,0.3,1.1])
plt.plot([1,1],plt.ylim(),'k:', plt.xlim(), [1,1], 'k:')
plt.plot([m1,m1],plt.ylim(),'k:', plt.xlim(), [c1,c1], 'k:')
plt.show()
plt.savefig(lib.rootdir+'/plots/paper/'+halo+'NFWcontour.eps')
			
plt.figure()
mlist=np.array(mlist)
plt.errorbar(rmax, mlist[:,0], mlist[:,1], fmt='o')
plt.xscale('log')
plt.plot([Halo.Rs0.value,Halo.Rs0.value],plt.ylim(), 'k--')
plt.plot([rv0, rv0], plt.ylim(), 'k--')
plt.plot(plt.xlim(),[1,1],'k:', plt.xlim(), [m1,m1],'k:')
plt.xlabel(r'$R_{max}$[kpc]')
plt.ylabel(r'$M/M_0$')
plt.savefig(lib.rootdir+'/plots/paper/'+halo+'NFWmass.eps')
#Sample.clean()
#FullSample.clean()
lib.close()