""" To generate mock halos 
this is never good. the data quality is always poor, no matter how many walkers or how many steps your run!
abaddoned!
"""
import sys,os
from dynio import *
#os.environ['OMP_NUM_THREADS']='32'

import emcee
from numpy import *
import h5py

def save_mock(flatchain, outfile="/work/Projects/DynDistr/data/mcmock.hdf5", sampleID=0):
  #random.shuffle(flatchain)
  x=flatchain[:,0:3]
  v=flatchain[:,3:]
  if sampleID==0:
    f=h5py.File(outfile,"w")
  else:
    f=h5py.File(outfile,"a")
  grp = f.create_group("/sample%d"%sampleID)
  grp.attrs['nwalker']=nwalkers
  grp.attrs['nsteps']=nsteps
  pos = grp.create_dataset("x", x.shape, dtype='f')
  vel = grp.create_dataset("v", v.shape, dtype='f')
  pos[...]=x
  vel[...]=v
  f.close()
  
#loglike=lambda x: dataprob(x[0],x[1],x[2])*x[0]*x[0]*x[2]
loglike=lambda x: dataprob6d(PARTYPEXV(x[0],x[1],x[2],x[3],x[4],x[5]))

#def loglike(x):
  #''' to be used with multiprocessing, since lambda function is not pickelable'''
  #return dataprob6d(PARTYPEXV(x[0],x[1],x[2],x[3],x[4],x[5]))

nwalkers=20
if len(sys.argv)>1:
  nwalkers=int(sys.argv[1])
nsteps=50000/nwalkers  
nburn=0
ndim=6
x00=zeros(ndim)
labels=[r"$x$",r"$y$",r"$z$",r"$v_x$",r"$v_y$",r"$v_z$"]

x0=kron(ones([nwalkers,1]),x00)#repmat to nwalkers rows
x0+=random.rand(ndim*nwalkers).reshape(nwalkers,ndim)-0.5 #random offset, [-0.5,0.5]
#x000=random.normal(size=(nwalkers,3))
#x0[:,0:3]=x000/kron(sqrt(sum(x000**2,axis=1)),ones([3,1])).T*1.1 #place x on a shell

#sampler=emcee.EnsembleSampler(nwalkers,ndim,loglike, threads=4)
sampler=emcee.EnsembleSampler(nwalkers,ndim,loglike)   
for i in range(5000):
  print i
  sampler.run_mcmc(x0,nsteps)
  x0=sampler.chain[:,-1,:]
  flatchain=swapaxes(sampler.chain, 0,1).reshape([-1,ndim])
  save_mock(flatchain,"/gpfs/data/jvbq85/DynDistr/data/mcmock_%dwalkers.hdf5"%nwalkers,i)
  sampler.clear_chain()

#pos,prob,state=sampler.run_mcmc(x0,2)
#raw_input("==========burn in finished. press any key to continue========")
#sampler.reset()
#sampler.run_mcmc(pos,nsteps)
"""
from matplotlib.pyplot import *
ion()
plot(range(0,nsteps), sampler.chain.mean(axis=0)[:,0:2])
figure()
for i in [0,3]:
  #subplot(ndim,1,i)
  subplot(2,1,i/3)
  for j in range(0,nwalkers,2): #range(nwalkers):
    plot(range(nsteps),sampler.chain[j,:,i],'-')
  ylabel(labels[i])
xlabel('Step') 
draw()
#show()

sample=sampler.chain[:,nburn:,:]
sample=swapaxes(sample, 0,1) #swap walker with steps, so sample is [nsteps, nwalkers, ndim], 
			      #and different walkers are in vicinity at the same step after reshaped
			      #but not necessary; different walkers are expected to be independent
flatchain=sample.reshape([-1,ndim])
save_mock(flatchain,"/work/Projects/DynDistr/data/mcmock_tmp.hdf5")

##now plot the TS
from scipy.stats import chi2,norm
def chi2sig(ts,dof=1):
  '''convert a chi2 TS to sigma'''
  pval=chi2.sf(ts,dof)
  sigma=norm.ppf(1.0-pval/2)
  return sigma

init()
figure()
samplesize=1000
samplestep=samplesize/nwalkers
sid=range(0,nsteps/samplestep,2)
y=[]
print "step, sigma"
for i in sid:
  select_particles(i)
  s=like2(1.,1.)
  y.append(s)
  print (i+1)*samplestep,s
stepid=(array(sid)+1)*samplestep
y=array(y)
#plot(stepid,chi2sig(y**2))
#xlabel('steps');ylabel('significance')
figure()
plot(stepid,abs(y))
xlabel('steps');ylabel('significance')
yscale('log')
figure()
plot(stepid,y)

#figure()
#plot(flatchain[:,0],flatchain[:,1],'.')
#xlabel('x')
#ylabel('y')

#for i in [1,3]:
    #figure()
    #hist(flatchain[:,i], 50, color="k", histtype="step")
    #title("Dimension {0:d}".format(i))
    #xlabel(labels[i])
#ylabel('Counts')    

#show()    
#import triangle
#fig=triangle.corner(flatchain,labels=labels,truths=median(flatchain,axis=0),quantiles=[0.5-0.683/2,0.5+0.683/2])
#fig.savefig("triangle.png")

#save_mock(flatchain,"/work/Projects/DynDistr/data/mcmock.hdf5")
"""