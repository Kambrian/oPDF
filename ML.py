import sys,os
from dynio import *
os.environ['OMP_NUM_THREADS']='32'


fullpars='mc'
if len(sys.argv)>1:
  freepars=sys.argv[1]
else:
  freepars=fullpars

init(-1)
select_particles(0)

from scipy.optimize import *
like= lambda x:-likefunc2(PARTYPE(x[0],x[1]))
x=fmin(like, [2,1.5])
print x
#numpy.gradient(f, *varargs)[source]
 
from iminuit import Minuit
from iminuit.ConsoleFrontend import ConsoleFrontend

#freepars=''
#freepars=fullpars 
fixnames=fullpars.translate(None,freepars)
fixdict={} #gen a dict for fixed pars
for x in fixnames: 
  fixdict['fix_'+x]=True
#p['par'][6]=1
defaults={'m':2,'c':1.5}
defaults.update(fixdict)
m=Minuit(neglike2,print_level=3,errordef=0.5,limit_m=[0.1,10],limit_c=[0.1,10],frontend=ConsoleFrontend(),**defaults)
#m=Minuit(like,print_level=3,errordef=0.5,frontend=ConsoleFrontend(),limit_a=[8,11],limit_b=[10,13],limit_c=[-4,-1],limit_d=[-2,0],**defaults) #for HOD

m.set_strategy(2);
m.tol=1e-3 #default convergence edm<1e-4*tol*errordef, but we do not need that high accuracy
result=m.migrad()
#m.print_param()
m.print_matrix()
#print m.values
#print m.covariance
#print m.matrix(correlation=True)
#print m.hesse
#print m.fitarg

L=-m.fval
L1=-neglike2(defaults['m'],defaults['c'])
print "current Like: ", L
print "Initial Like: ", L1
print "TS %.1f"%(2*(L-L1))
from scipy.stats import chi2,norm
pval=chi2.sf(2*(L-L1),len(m.list_of_vary_param()))
sigma=norm.ppf(1.0-pval/2)
print "Pval= %g, Significance= %.1f-sigma, DoF=%d"%(pval,sigma,len(m.list_of_vary_param()))

from matplotlib import *
#use('agg')
from matplotlib.pyplot import *

figure()
m.draw_profile('m')
show()

#figure()
#m.errordef=1.15 #for two parameter CI.
#m.draw_contour('m','c')
#show()
#savefig('M-c-contour.eps');
#sys.stdout = open(outdir+'log.plt', 'w')    
"""
from matplotlib import *
use('agg')
from matplotlib.pyplot import *

outdir='plots/'+p['mockname']+'/'
if not os.path.exists(outdir):
    os.makedirs(outdir)
figure()
m.draw_profile('a') #the error on the plot cannot be trusted at all, if there are strong covariance among the parameters
		    # but m.errors and m.covariance can be trusted, and are always quite close to m.minos() estimation
savefig(outdir+"MLprof_a.eps")

savefig(outdir+"MLprof_b.eps")
figure()
m.draw_profile('c')
savefig(outdir+"MLprof_c.eps")

savefig(outdir+"MLcontour_a_c.eps")
figure()
m.draw_contour('a','b')
savefig(outdir+"MLcontour_a_b.eps")
figure()
m.draw_contour('b','c')
savefig(outdir+"MLcontour_b_c.eps")
"""
