import sys,os
from gama_io import *
os.environ['OMP_NUM_THREADS']='32'

paramfile=sys.argv[1]
p=load_params(paramfile)
srcfile='shearmap_ZEBRA.mat'

GAMA_init(srcfile)

loglike=gen_func2(p['modelid'],p['npar']) #logProb

import emcee
from numpy import *

nwalkers=10
nsteps=1000
nburn=100
ndim=p['npar']
x00=array(p['par'])

x0=kron(ones([nwalkers,1]),x00)#repmat to nwalkders rows
x0+=random.rand(ndim*nwalkers).reshape(nwalkers,ndim)-0.5 #random offset, [-0.5,0.5]

sampler=emcee.EnsembleSampler(nwalkers,ndim,loglike)
sampler.run_mcmc(x0,nsteps)
#pos,prob,state=sampler.run_mcmc(x0,2)
#raw_input("==========burn in finished. press any key to continue========")
#sampler.reset()
#sampler.run_mcmc(pos,nsteps)

from matplotlib.pyplot import *
figure()
for i in range(ndim):
  subplot(ndim,1,i)
  for j in range(nwalkers):
    plot(range(nsteps),sampler.chain[j,:,i],'.')

sample=sampler.chain[:,nburn:,:]
flatchain=sample.reshape([-1,ndim])
    	
for i in range(ndim):
    figure()
    hist(flatchain[:,i], 50, color="k", histtype="step")
    title("Dimension {0:d}".format(i))
    #plt.xlim(-1000,1000)

show()    
#import triangle
fig=triangle.corner(flatchain,labels=["$\log M_p$",r"$\alpha$"],truths=median(flatchain,axis=0),quantiles=[0.5-0.683/2,0.5+0.683/2])
#fig.savefig("triangle.png")

