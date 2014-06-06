from TSmap import *

def filt_TS(lim=[-inf, -40], nskip=50):
	P.filter_TSmap(TS,x, lim)
	f=P.P['flag']>0
	#P.filter_TSmap(TS,x, [16, inf])
	#f=(f+P.P['flag'])>0
	figure()
	#plot(P.P['r'][f][0:-1:nskip], P.P['vr'][f][0:-1:nskip],'.')
	data=P.P[['r','vr']].view(float).reshape(P.nP.value,-1).T
	percentile_contour(data[:,f][:,0:-1:nskip], color=None, percents=[0.05,0.1,0.3,0.6,0.9])
	for i in range(500):
		plot(S.r[i], S.vr[i], 'o', markerfacecolor='w', markersize=40./(i+1.)**0.6)
	return f	
#where((S.x[:,0]>68)*(S.x[:,0]<69)*(S.x[:,1]>-112.2)*(S.x[:,1]<-111.8))

def plot_contour(f, dims=[0,1], nskip=50):
	figure();
	fs=S.m>100
	percents=[0.1,0.3,0.6,0.9]
	subplot(121)
	percentile_contour(P.P['x'][0:-1:800,dims].T, percents=percents, color='k', alpha=0.3)
	percentile_contour(P.P['x'][f][0:-1:nskip,dims].T, color='r', percents=percents)
	plot(S.x[fs,dims[0]], S.x[fs,dims[1]], 'ro', alpha=0.1)
	plot_circle(r=float(os.environ['DynRv']),linestyle=':', color='k')
	xlabel('x[kpc]');ylabel('y[kpc]')
	axis('equal')
	ax=subplot(122)
	data=P.P[['r','vr']].view(float).reshape(P.nP.value,-1).T
	percentile_contour(data[:,0:-1:800], percents=percents, color='k', alpha=0.3)
	percentile_contour(data[:,f][:,0:-1:nskip], color='r', percents=percents)
	plot(S.r[fs], S.vr[fs], 'ro', alpha=0.1)
	xlabel('r[kpc]');ylabel(r'$v_r$[km/s]')
	ax.yaxis.tick_right()
	ax.yaxis.set_ticks_position('both')
	ax.yaxis.set_label_position('right')
	subplots_adjust(wspace=0.1, right=0.88)
	savefig(rootdir+'/plots/full_halo/R1_200/100x100/Check'+halo+'_contour.png')
	
def plot_selection(f, nskip=50, dims=[0,1]):
	figure();
	subid=[415,20]
	fs=S.m>100
	subplot(121)
	plot(P.P['x'][0:-1:1000,dims[0]], P.P['x'][0:-1:1000,dims[1]],'k.', alpha=0.1)
	plot(S.x[fs,dims[0]], S.x[fs,dims[1]], 'ro', alpha=0.1)
	plot(P.P['x'][f,dims[0]][0:-1:nskip], P.P['x'][f,dims[1]][0:-1:nskip],'.')
	#plot(S.x[subid,0], S.x[subid,1], 'ro')
	plot_circle(r=float(os.environ['DynRv']),linestyle='--', color='k')
	xlabel('x[kpc]');ylabel('y[kpc]')
	axis('equal')
	ax=subplot(122)
	plot(P.P['r'][0:-1:1000], P.P['vr'][0:-1:1000],'k.',alpha=0.1)
	plot(S.r[fs], S.vr[fs], 'ro', alpha=0.1)
	plot(P.P['r'][f][0:-1:nskip], P.P['vr'][f][0:-1:nskip],'.')
	#plot(S.r[subid], S.vr[subid], 'ro')
	xlabel('r[kpc]');ylabel(r'$v_r$[km/s]')
	ax.yaxis.tick_right()
	ax.yaxis.set_ticks_position('both')
	ax.yaxis.set_label_position('right')
	subplots_adjust(wspace=0.1, right=0.88)
	savefig(rootdir+'/plots/full_halo/R1_200/100x100/Check'+halo+'.png')

def plot_3d(f, nskip=50):  
	from mayavi import mlab
	mlab.figure()
	mlab.points3d(P.P['x'][f,0], P.P['x'][f,1], P.P['x'][f,2], mask_points=nskip, color=(1,0,0), scale_factor=10)
	#mlab.points3d(S.x[415][0], S.x[415][1], S.x[415][2], color=(0,1,0), scale_factor=20)
	mlab.figure()
	mlab.points3d(P.P['v'][f,0], P.P['v'][f,1], P.P['v'][f,2], mask_points=nskip, color=(1,0,0), scale_factor=10)
	#mlab.points3d(S.v[415][0], S.v[415][1], S.v[415][2], color=(0,1,0), scale_factor=20)

def plot_phase_hist(nbin=100, name=''):
	P.squeeze_data()
	elike(1,1,10)
	figure()
	hist(P.P['theta'],nbin,normed=True, histtype='step')
	plot([0,1],[1,1],'k--')
	xlabel(r'$\theta$')
	ylabel(r'$dP/d\theta$')
	savefig(rootdir+'/plots/full_halo/R1_200/100x100/Check'+halo+'_phase'+name+'.png')

def plot_radial_hist(nbin=100, name=''):
	'''either call plot_phase_hist before this func,
	or call P.squeeze_data() and elike(1,1,4) first.'''
	pred=predict_radial_count(nbin)
	count,bins=histogram(P.P['r'],bins=linspace(P.R_MIN.value,P.R_MAX.value,nbin+1))
	step(bins,hstack((count,nan)), where='post', color='r')
	step(bins,hstack((pred,nan)), where='post', color='k', linestyle='--')
	legend(('Data','Model'))
	xlabel('r[kpc]')
	ylabel('Counts')
	savefig(rootdir+'/plots/full_halo/R1_200/100x100/Check'+halo+'_radial'+name+'.png')
	
if __name__=='__main__':	
	ion()
	halo='AqA4subFit8'
	if len(sys.argv)>1:
		halo=sys.argv[1]  
	get_config(halo)
	#os.environ['DynSIZE']='10000'
	os.environ['DynRMAX']='200'
	os.environ['DynRMIN']='1'
	init()
	P=pick_particles(0)
	try:
		S=SubData(rootdir+'/data/'+halo[0:4]+'sublist.hdf5')
	except:
		print "Error: no subhalolist available"
		raise
	TS,x=load_TSmap(rootdir+'/plots/full_halo/R1_200/100x100/Scan_'+halo+'_Mean_E_L2_200.hdf5')
	#hist(TS[~np.isnan(TS)],100)

	tsmin=nanmin(TS)
	tsmax=nanmax(TS)
	if -tsmin>tsmax:
		lim=[tsmin-10,tsmin+0.01]
	else:
		lim=[tsmax-0.01,tsmax+10]
	#lim=[-2,2]
	f=filt_TS(lim)
	plot_contour(f,nskip=1)
	plot_selection(f,nskip=1)
	#plot_3d(f)	
	plot_phase_hist()
	plot_radial_hist()
	P=pick_particles()
	P.filter_TSmap(TS,x,[-1,1])
	f=P.P['flag']>0
	plot_contour(f,nskip=100)
	plot_selection(f,nskip=100)
	plot_phase_hist(name='_good')
	plot_radial_hist(name='_good')
	
	
	

