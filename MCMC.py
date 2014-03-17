import sys,os
from dynio import *
os.environ['OMP_NUM_THREADS']='32'

init()

loglike=lambda par: likefunc2(PARTYPE(par[0],par[1]))

import emcee
from numpy import *

nwalkers=4
nsteps=500
nburn=300
ndim=2
x00=array([0.,0.])
labels=[r"$\log (\rho_s/\rho_{s0})$",r"$\log (r_s/r_{s0})$"]

x0=kron(ones([nwalkers,1]),x00)#repmat to nwalkers rows
x0+=random.rand(ndim*nwalkers).reshape(nwalkers,ndim)-0.5 #random offset, [-0.5,0.5]

sampler=emcee.EnsembleSampler(nwalkers,ndim,loglike)
sampler.run_mcmc(x0,nsteps)
#pos,prob,state=sampler.run_mcmc(x0,2)
#raw_input("==========burn in finished. press any key to continue========")
#sampler.reset()
#sampler.run_mcmc(pos,nsteps)

from matplotlib.pyplot import *
ion()
figure()
for i in range(ndim):
  subplot(ndim,1,i)
  for j in range(nwalkers):
    plot(range(nsteps),sampler.chain[j,:,i],'.')
  ylabel(labels[i])
xlabel('Step')  

sample=sampler.chain[:,nburn:,:]
flatchain=sample.reshape([-1,ndim])
    	
for i in range(ndim):
    figure()
    hist(flatchain[:,i], 50, color="k", histtype="step")
    title("Dimension {0:d}".format(i))
    #plt.xlim(-1000,1000)
    xlabel(labels[i])
ylabel('Counts')    

show()    
import triangle
fig=triangle.corner(flatchain,labels=labels,truths=median(flatchain,axis=0),quantiles=[0.5-0.683/2,0.5+0.683/2])
#fig.savefig("triangle.png")

