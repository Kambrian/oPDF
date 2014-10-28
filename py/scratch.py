# todo: stream-cleaned version
# star 500kpc
from dynio import *
from myutils import *
import h5py,os,sys
from scipy.stats import chi2
plt.ion()

with Tracer('A2', DynDataFile='A2star.hdf5', DynRMIN=10, DynRMAX=300) as F:
  with Tracer('A2', DynDataFile='A2starCleanFromDM.hdf5', DynRMIN=10, DynRMAX=300) as F2:
	lib.init_potential_spline()
	F.freeze_energy()
	F2.freeze_energy()
	F.like_init()
	F2.like_init()

	S=F.select(F.data['E']>0)
	S2=F2.select(F2.data['E']>0)	

SS=S.select(S.data['subid']<=0)
SS2=S2.select(S2.data['subid']<=0)
x,y,z=density_of_points(np.log10(np.array([S.data['E'], S.data['L2']])), method='hist')
x2,y2,z2=density_of_points(np.log10(np.array([S2.data['E'], S2.data['L2']])), method='hist')
plt.contour(x,y,z, colors='r')
plt.contour(x2,y2,z2, colors='g')
plt.figure()
x,y,z=density_of_points(np.log10(np.array([SS.data['E'], SS.data['L2']])), method='hist')
x2,y2,z2=density_of_points(np.log10(np.array([SS2.data['E'], SS2.data['L2']])), method='hist')
plt.contour(x,y,z, colors='r')
plt.contour(x2,y2,z2, colors='g')
