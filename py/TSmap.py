import sys,os,ConfigParser,h5py
from collections import OrderedDict
from numpy import *
import matplotlib
#import __main__ as main
#if hasattr(main, '__file__'): #non-interactive
#matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 15, 'ps.fonttype' : 42 , 'pdf.fonttype' : 42 ,' image.origin': 'lower', 'image.interpolation': None})
from matplotlib.pyplot import *
from myutils import *
from dynio import *
#os.environ['OMP_NUM_THREADS']='32'

def load_TSmap(file):
  f=h5py.File(file,'r')
  bintype=f['/bintype']
  x=OrderedDict()
  for b in bintype:
    x[b]=f[b][:]
  TS=f['/TS'][:]
  return TS,x
 
def pick_particles(sample=0):
  select_particles(sample)
  P=ParticleData()
  freeze_energy(1,1)
  return P

def scanTS2D(bintype=('E','r'), sample=0, nbin=30, estimator=8, logscale=True, plotcount=False):
  P=pick_particles(sample)  
  #p=P.StructP #original copy
  rmin=P.R_MIN.value
  rmax=P.R_MAX.value
  proxy=OrderedDict()
  x=OrderedDict()
  for b in bintype:
    x[b],proxy[b]=P.gen_bin(b,nbin,logscale)
  TS=empty([nbin,nbin])
  for i,x0 in enumerate(x[bintype[0]][:-1]):
    x1=x[bintype[0]][i+1]
    fx=(proxy[bintype[0]]>=x0)*(proxy[bintype[0]]<x1)
    if bintype[0]=='r': #update radial limits for radial cut
		P.R_MIN.value=x0
		P.R_MAX.value=x1
    for j,y0 in enumerate(x[bintype[1]][:-1]):
      y1=x[bintype[1]][j+1]
      fy=(proxy[bintype[1]]>=y0)*(proxy[bintype[1]]<y1)
      if sum(fx*fy)==0:
		TS[i,j]=nan
		continue
      #P.P2P.contents=(Particle*P.nP.value).from_buffer_copy(p) #make a copy to the c-pointer
      #P=ParticleData() #reinitialize the particle data to be the new copy
      P.P['flag']=0
      P.P['flag'][fx*fy]=1
      if bintype[1]=='r': #update radial limits for radial cut
		P.R_MIN.value=y0
		P.R_MAX.value=y1
      P.squeeze_data()
      #if (array(bintype)=='r').any(): #update the bin counts for radial cut
	#lib.fill_radial_bin() #TODO: uncomment this when estimator=4
      TS[i,j]=elike(1,1,estimator)
      print i,j,P.nP.value, TS[i,j]
      P.R_MIN.value=rmin #restore radial limits
      P.R_MAX.value=rmax
      P=pick_particles(sample)  #restore sample
  if logspace:
    extent=[log10(x[b][j]) for b in bintype for j in [0,-1]]
  else:
    extent=[x[b][j] for b in bintype for j in [0,-1]]
  if plotcount:
    figure()
    subplot(1,2,1)
    H=histogram2d(proxy[bintype[0]],proxy[bintype[1]],[x[b] for b in bintype])
    imshow(H[0].T,extent=extent,origin='lower', interpolation='None')
    subplot(1,2,2)
  if estimator==8:
    TS=AD2Sig(-TS)  
  imshow(TS.T,extent=extent,origin='lower', interpolation='None')
  colorbar()
  #CS=contour(AD2Sig(-TS.T),levels=[1,3,5,7],extent=extent,origin='lower')
  #clabel(CS,inline=1)
  if logscale:
    xlabel(r'$\log$(%s)'%bintype[0])
    ylabel(r'$\log$(%s)'%bintype[1])
  else:  
    xlabel(bintype[0])
    ylabel(bintype[1])
  return TS,x,proxy

def TSprof(bintype='E',sample=0,nbin=100,estimator=8,logscale=True,plotcount=True,equalcount=True):
  P=pick_particles(sample)
  rmin=P.R_MIN.value
  rmax=P.R_MAX.value
  try:
    len(nbin)
    x=nbin
    proxy=copy(P.P[bintype])
  except:
    x,proxy=P.gen_bin(bintype,nbin,logscale,equalcount=equalcount)
  TS=[]
  for i,x0 in enumerate(x[:-1]):
    x1=x[i+1]
    fx=(proxy>=x0)*(proxy<x1)
    if sum(fx)==0:
      TS.append(nan)
      continue
    P.P['flag']=0
    P.P['flag'][fx]=1
    if bintype=='r': #update radial limits for radial cut
      P.R_MIN.value=x0
      P.R_MAX.value=x1
    P.squeeze_data() #radial limits actually get updated inside this func.
    #if bintype=='r': #update the bin counts for radial cut
      #lib.fill_radial_bin() #TODO: uncomment this when estimator=4
    TS.append(elike(1,1,estimator))
    print i, P.nP.value, TS[-1]
    P.R_MIN.value=rmin #restore radial limits
    P.R_MAX.value=rmax
    P=pick_particles(sample) 
  
  if plotcount:
    figure()
    subplot(211)
    hist(proxy,x,log=logscale)
    if logscale:
      xscale('log')
    xlim(x.min(),x.max())
    ylabel('Count')
    subplot(212)
  #plot(x[:-1],AD2Sig(-array(TS)))
  #plot(x[:-1],(-array(TS)))
  TS.append(nan)
  TS=array(TS)
  if estimator==8:
    TS=AD2Sig(-TS)
  step(x,TS, where='post')
  if logscale:
    xscale('log')
  xlim(x.min(),x.max())
  xlabel(bintype)
  ylabel('TS')
  #ylabel(r'Discrepancy/$\sigma$')
  if bintype=='r':
    rconv=float(os.environ['DynRconv'])
    rv=float(os.environ['DynRv'])
    plot([rconv,rconv],ylim(),'k:')
    plot([rv,rv],ylim(),'k:')
  return TS,x

def plot_halo_TS(halo,estimator=8, nbin=100, rmin=None, rmax=None, equalcount=True, flagsave=False):    
  get_config(halo)
  if rmin!=None:
    os.environ['DynRMIN']=format(rmin)
  else:
    rmin=float(os.environ['DynRMIN'])
  if rmax!=None:
    os.environ['DynRMAX']=format(rmax)
  else:
    rmax=float(os.environ['DynRMAX'])
  init()
  x=dict()
  TS=dict()
  figure()
  subplot(311)
  TS['r'],x['r']=TSprof(bintype='r',estimator=estimator,nbin=nbin, equalcount=equalcount, plotcount=False)
  title(" ".join([halo,NameList[estimator]]))
  xlim(rmin,rmax)
  subplot(312)
  TS['E'],x['E']=TSprof(bintype='E',estimator=estimator,nbin=nbin, equalcount=equalcount, plotcount=False)
  subplot(313)
  TS['L2'],x['L2']=TSprof(bintype='L2',estimator=estimator,nbin=nbin, equalcount=equalcount, plotcount=False)
  tight_layout()
  if flagsave:
    outfile=rootdir+'plots/'+'_'.join(['TS',halo,NameList[estimator],'R%g'%rmin, format(rmax,'.0f')])
    savefig(outfile+'.eps')
    f=h5py.File(outfile+'.hdf5', 'w')
    for b in ['r','E','L2']:
      group=f.create_group(b)
      group.create_dataset('x',data=x[b])
      group.create_dataset('TS',data=TS[b])
    f.close()  
  free_data()
  return TS,x
    
def plot_halo_scan(halo,bintype,rmax=None, estimator=10, nbin=30, flagsave=False):  
  get_config(halo)
  if rmax!=None:
    os.environ['DynRMAX']=format(rmax)
  else:
    rmax=float(os.environ['DynRMAX'])
  init()
  TS,x,proxy=scanTS2D(bintype,nbin=nbin,estimator=estimator)
  title(halo+' '+NameList[estimator])
  axis('tight')
  if flagsave:
    outfile=rootdir+'plots/'+'_'.join(['Scan',halo,NameList[estimator]]+list(bintype))+'_'+format(rmax,'.0f')
    savefig(outfile+'.eps')
    f=h5py.File(outfile+'.hdf5', 'w')
    f.create_dataset('bintype',data=bintype)
    for b in bintype:
      f.create_dataset(b,data=x[b])
    f.create_dataset('TS',data=TS)
    f.close()   
  free_data()
  return TS,x

def SubOnTSmap(halo, bintype=('E','L2'), mmin=100, scmap=cm.gist_rainbow, shiftcmap=False, markercolor='k', interpolation='None'):
  "plot subhaloes on top of TS contour"
  get_config(halo)
  init()
  P=pick_particles(0)
  free_data()
  S=SubData(rootdir+'/data/'+halo[0:4]+'sublist.hdf5')
  S.eval_energy()
  TS,x=load_TSmap(rootdir+'/plots/full_halo/R1_200/100x100/Scan_'+halo+'_Mean_'+bintype[0]+'_'+bintype[1]+'_%.0f.hdf5'%P.R_MAX.value)
  extent=[log10(b[j]) for b in x.values() for j in [0,-1]]
  figure()
  if shiftcmap: #so that TS=0 is at the center of the colormap
	tsmin=nanmin(TS)
	tsmax=nanmax(TS)
	if tsmin<0 and tsmax>0:
	  scmap=shiftedColorMap(scmap, midpoint=-tsmin/(tsmax-tsmin))
  imshow(TS.T, extent=extent, origin='lower', cmap=scmap, interpolation=interpolation)
  #contour(TS.T, extent=extent, origin='lower', cmap=cm.jet)
  colorbar()
  xlabel('log('+bintype[0]+')')
  ylabel('log('+bintype[1]+')')
  name2data={'E': -S.E, 'L2': S.L2, 'r': S.r}
  for i in range(sum(S.m>mmin)):
	  if S.r[i]<P.R_MAX and S.r[i]>P.R_MIN:
		plot(np.log10(name2data[bintype[0]][i]), np.log10(name2data[bintype[1]][i]), 'o', markeredgecolor=markercolor, markeredgewidth=1, markerfacecolor='none', markersize=S.m[i]**(1./3)/1.5)
  axis('tight')
  #xlim(1,xlim()[1])
  outfile=rootdir+'plots/'+'_'.join(['Scan',halo,'Mean']+list(bintype))+'_'+format(P.R_MAX.value,'.0f')+'_sub'
  savefig(outfile+'.eps')
  
if __name__=="__main__":
 
  #ion()
  halo='AqA4'
  estimator=10
  if len(sys.argv)>1:
    halo=sys.argv[1]  
  if len(sys.argv)>2:
    estimator=int(sys.argv[2])
 
  #get_config(halo)
  #os.environ['DynRMAX']='100'
  #init()  
  #TSprof('r', estimator=estimator, nbin=200)
  #show()
  #TSprof('r', estimator=estimator, nbin=200, equalcount=True)
  #for halo in ['AqA4Fit4','AqA4Fit8','AqA4Fit10','AqA4LargeRFit4','AqA4LargeRFit8','AqA4LargeRFit10']:  
  #plot_halo_TS(halo,estimator=10,flagsave=True)
  #plot_halo_TS(halo,estimator=8,rmin=0.01, rmax=1000, flagsave=True)
  #plot_halo_TS(halo,estimator=8,rmin=1, rmax=200, flagsave=True)
  #plot_halo_TS(halo,estimator=10,rmin=0.01, rmax=1000, flagsave=True)
  #plot_halo_TS(halo,estimator=10,rmin=1, rmax=200, flagsave=True)
  #rmax=100
  #tsmap=dict()
  #for b in [('E','L2'),('L2','r'),('E','r')]:
    #figure();
    #tsmap[b]=plot_halo_scan(halo,b,estimator=10,flagsave=True)
  
  #tsmap=dict()
  for halo in ['AqA4Fit4','AqA4Fit8','AqA4subN','AqA4subFit8','AqB4Fit8','AqB4subFit8']:
	  plot_halo_TS(halo, estimator=10, rmin=1, rmax=200, flagsave=True)
    #figure();
    #tsmap[halo]=plot_halo_scan(halo,('E','L2'),estimator=10,flagsave=True)
  #show()
  #init()
  #P=pick_particles(0)
  #x=dict()
  #p=dict()
  #for b in ['E','L2','r']:
    #x[b],p[b]=P.gen_bin(b)  
