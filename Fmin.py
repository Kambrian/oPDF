import sys,os
from dynio import *
from iminuit import Minuit
from iminuit.ConsoleFrontend import ConsoleFrontend
#os.environ['OMP_NUM_THREADS']='32'

halo='AqA4N'
if len(sys.argv)>1:
  halo=sys.argv[1]

get_config(halo)
#os.environ['DynSIZE']='10000'
os.environ['DynRMIN']='1'
os.environ['DynRMAX']='200'
init()
select_particles(0)

from scipy.optimize import *
for estimator in [4,8,10,11]:
  print '------------- estimator ', estimator, '----------------'
  like= lambda x:-likefunc2(PARTYPE(x[0],x[1]), estimator)
  #result=fmin(like, 10**(random.rand(2)-0.5), xtol=0.001, ftol=1e-4, maxiter=1000, maxfun=5000, full_output=True)
  x=fmin(like, [1,1], xtol=0.01, ftol=1e-3)
  print x
  print '================================='
 
def minuit_like(estimator=8):
  neglike2=lambda m,c: -likefunc2(PARTYPE(m,c),estimator)
  defaults={'m':1,'c':1}
  m=Minuit(neglike2, print_level=3,errordef=0.5,limit_m=[0.1,10],limit_c=[0.1,10],frontend=ConsoleFrontend(),**defaults)
  #m.set_strategy(2);
  m.tol=1 #default convergence edm<1e-4*tol*errordef, but we do not need that high accuracy
  result=m.migrad()
  m.print_matrix()
  return result,m