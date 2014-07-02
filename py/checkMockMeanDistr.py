#check the distribution of mock means
from dynio import *
from myutils import *
import matplotlib.pyplot as plt

npart=1000

lib.open()
mean=[]
AD=[]
with Tracer('Mock') as FullSample:
  for sampleid in range(750):
	with FullSample.copy(sampleid*npart,npart) as Sample:
	  mean.append(Sample.freeze_and_like([1,1], 14))
	  AD.append(-Sample.likelihood([1,1], 8))
lib.close()
plt.figure()
x=plt.hist(mean,30, normed=1, histtype='step')
h,=plt.plot(x[1], norm.pdf(x[1], 0, 1), 'r')
plt.xlabel(r'$\bar{\Theta}$')
plt.ylabel(r'$dP/d\bar{\Theta}$')
plt.legend((x[2][0],h),('Mock','Model'))
plt.savefig(lib.rootdir+'/plots/paper/extra/MockMeanPhaseDistr.eps')

plt.figure()
x=plt.hist(np.log(AD),30, normed=1, histtype='step')
h,=plt.plot(x[1], lnADPDF(x[1]), 'r')
plt.xlabel(r'$\ln(D)$')
plt.ylabel(r'$dP/d\ln(D)$')
plt.legend((x[2][0],h),('Mock','Model'))
plt.savefig(lib.rootdir+'/plots/paper/extra/MockADDistr.eps')
