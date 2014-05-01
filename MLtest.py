#to do statistical test rather than max-like
#differ in errordef and dof, ts, sigma etc.
import sys,os
from dynio import *
os.environ['OMP_NUM_THREADS']='32'


fullpars='mc'
if len(sys.argv)>1:
  freepars=sys.argv[1]
else:
  freepars=fullpars

init()
select_particles(0)

from iminuit import Minuit
from iminuit.ConsoleFrontend import ConsoleFrontend

#freepars=''
#freepars=fullpars 
fixnames=fullpars.translate(None,freepars)
fixdict={} #gen a dict for fixed pars
for x in fixnames: 
  fixdict['fix_'+x]=True
#p['par'][6]=1
defaults={'m':1,'c':1}
defaults.update(fixdict)
m=Minuit(neglike2,print_level=3,errordef=1.0, frontend=ConsoleFrontend(),**defaults)
#m=Minuit(like,print_level=3,errordef=0.5,frontend=ConsoleFrontend(),limit_a=[8,11],limit_b=[10,13],limit_c=[-4,-1],limit_d=[-2,0],**defaults) #for HOD

m.tol=10. #default convergence mdm<1e-4*tol*errordef, but we do not need that high accuracy
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
pval=chi2.sf(-L1,1)
sigma=norm.ppf(1.0-pval/2)
print "Pval= %g, Significance= %.1f-sigma, DoF=%d"%(pval,sigma,1)

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
