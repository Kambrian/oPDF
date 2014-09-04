""" density profile 
NFW should be a better description of the total rather than the FoF profile.
The difference is only in the outer region, and affects little the mass
"""
from dynio import *
from myutils import *
import h5py,os,sys
from scipy.stats import chi2
import itertools
from matplotlib import gridspec
plt.ion()

#estimator=8
#proxy='LE'
#nbins=[10,10]
npart=10000

halo='AqB4'
alpha=2 #plot rho*r^alpha
lib.open()
c0=get_config(halo)
rconv=float(c0['DynRconv'])
Halo=NFWHalo()

gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1]) 
fig=plt.figure(figsize=(8,8))
ax0=plt.subplot(gs[0])
  
xbin=np.logspace(0, log10(400), 100)
xcen=xbin/np.sqrt(xbin[1]/xbin[0])
vol=np.diff(np.hstack([0., xbin])**3)*np.pi*4/3
with Tracer(halo+'all',DynRMIN=0, DynRMAX=400) as FullSample:
  count,tmp=np.histogram(FullSample.data['r'], np.hstack([0., xbin]))
  density=count*FullSample.mP/vol
  plt.plot(xcen, density*xcen**alpha, '-', color=[0.6,0.6,0.6], linewidth=5)  
  
rmax=[30, 50,100,300] #,rv0,rv1]
ColorList=itertools.cycle(plt.cm.jet(np.linspace(0.,1.,len(rmax)+2)))
ColorList.next()
hall=[]
ylist=[]
for DynRMAX in rmax:
  with Tracer(halo+'all',DynRMIN=rconv, DynRMAX=DynRMAX) as FullSample:
	with FullSample.copy(0,npart) as Sample:
	  Sample.mP*=float(FullSample.nP)/Sample.nP
	  #Sample.radial_count()
	  result,m=Sample.minuit_NFWlike()
	  color=ColorList.next()
	  y=Halo.halo.Rhos/(xcen/Halo.halo.Rs)/(1+xcen/Halo.halo.Rs)**2
	  ylist.append(y)
	  y=y*xcen**alpha
	  filt=(xcen>rconv)*(xcen<=DynRMAX)
	  htmp,=plt.plot(xcen[filt], y[filt], '-', color=color, linewidth=1)
	  plt.plot(xcen, y , '--', color=color, linewidth=1)
	  hall.append(htmp)
	  plt.plot(Halo.halo.Rv, Halo.halo.Rhos/Halo.halo.c/(1+Halo.halo.c)**2*Halo.halo.Rv**alpha, 'x', color=color)

c0=get_config(halo)
Halo.M0.value=float(c0['DynM0'])
Halo.C0.value=float(c0['DynC0'])
Halo.define_halo([1,1])
rs=Halo.halo.Rs
rv0=Halo.halo.Rv
m0=Halo.M0.value
y=Halo.halo.Rhos/(xcen/Halo.halo.Rs)/(1+xcen/Halo.halo.Rs)**2
ylist.append(y)
h0,=plt.plot(xcen, y*xcen**alpha, '-', color=ColorList.next(), linewidth=1)

plt.ylabel(r'$\rho r^2[10^{10} M_\odot/\mathrm{kpc}]$')
plt.xscale('log')
plt.yscale('log')
plt.xlim([rconv,300])
if halo=='AqA4':
  c1=get_config(halo+'N')
  Halo.M0.value=float(c1['DynM0'])
  Halo.C0.value=float(c1['DynC0'])
  Halo.define_halo([1,1])
  rs=Halo.halo.Rs
  rv1=Halo.halo.Rv
  m1=float(c1['DynM0'])
  y=Halo.halo.Rhos/(xcen/Halo.halo.Rs)/(1+xcen/Halo.halo.Rs)**2
  ylist.append(y)
  h1,=plt.plot(xcen, y*xcen**alpha, '-', color=ColorList.next(), linewidth=1)
  #x=np.logspace(2,np.log10(300),20)
  #h2,=plt.plot(plist[0], plist[1], 'k--')
  hall.extend([h0,h1])
  plt.legend(hall,('30', '50','100','300','FitS','FitN'), loc='lower center', ncol=2, fontsize=15)
  plt.ylim([0.009,0.2])
if halo=='AqB4':
  #h2,=plt.plot(plist[0], plist[1], 'k--')
  hall.extend([h0])
  plt.legend(hall,('30', '50','100','300','FitS'), loc='lower center', ncol=2, fontsize=15)
  plt.ylim([0.008,0.1])  
#plt.plot([rv0,rv0],[plt.ylim()[0], m0/rv0**alpha],'k:', [plt.xlim()[0], rv0], [m0/rv0**alpha,m0/rv0**alpha], 'k:')
#plt.plot([rv1,rv1],[plt.ylim()[0], m1/rv1**alpha],'k:', [plt.xlim()[0], rv1], [m1/rv1**alpha,m1/rv1**alpha], 'k:')
plt.plot([rs,rs],plt.ylim(),'k:') #,[rconv,rconv],plt.ylim(),'k:')

ax1=plt.subplot(gs[1])
ColorList=itertools.cycle(plt.cm.jet(np.linspace(0.,1.,len(rmax)+2)))
ColorList.next()
for i,y in enumerate(ylist):
  color=ColorList.next()
  y=(y/density)
  if i<len(rmax):
	filt=(xcen>rconv)*(xcen<=rmax[i])
	plt.plot(xcen[filt], y[filt], '-', color=color)
	plt.plot(xcen, y, '--', color=color)
  else:
	plt.plot(xcen, y, '-', color=color)
	
plt.xlabel(r'$R$[kpc]')
plt.ylabel(r'$\rho_{fit}/\rho_{data}$')
plt.xscale('log')
plt.xlim([rconv,300])
if halo=='AqA4':
  plt.ylim([0.5,1.5])
  plt.yticks(np.arange(0.7,1.4,0.3))
if halo=='AqB4':
  plt.ylim([0.9,1.15])
  plt.yticks(np.arange(0.9,1.13,0.1))  
plt.plot(plt.xlim(), [1,1], 'k:')
plt.plot([rs,rs],plt.ylim(),'k:')
#,[rconv,rconv],plt.ylim(),'k:')
fig.subplots_adjust(left=0.15, hspace=0)
plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
#plt.setp(ax[1].yaxis, major_locator=MaxNLocator(prune='upper'))
plt.minorticks_on()

plt.show()
plt.savefig(lib.rootdir+'/plots/paper/'+halo+'NFWDensityProf.eps')
			
#Sample.clean()
#FullSample.clean()
lib.close()