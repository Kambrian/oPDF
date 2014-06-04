import sys,os
from dynio import *
from iminuit import Minuit
from iminuit.ConsoleFrontend import ConsoleFrontend
#os.environ['OMP_NUM_THREADS']='32'

def minuit_like(estimator=8,x0=[1,1]):
	#too difficult for minuit to work.
	#just use this to estimate the error matrix
	#set a huge tol to avoid walking away
	neglike2=lambda m,c: -likefunc2(PARTYPE(m,c),estimator)
	defaults={'m':x0[0],'c':x0[1]}
	m=Minuit(neglike2, print_level=0,pedantic=False, error_m=0.1, error_c=0.1, errordef=0.5,limit_m=[0.1,10],limit_c=[0.1,10],frontend=ConsoleFrontend(),**defaults)
	#m.set_strategy(2);
	m.tol=1e10    #default convergence edm<1e-4*tol*errordef, but we do not need that high accuracy
	result=m.migrad()
	m.print_param()
	m.print_matrix()
	return result,m

halo='AqA4N'
if len(sys.argv)>1:
  halo=sys.argv[1]

get_config(halo)
os.environ['DynSIZE']='10000'
os.environ['DynRMIN']='1'
os.environ['DynRMAX']='200'
init()
select_particles(0)

from scipy.optimize import *

estimator=8
for i in range(5):
	x0=10**(np.random.rand(2)-0.5)
	print '------------- x0= ', x0, '----------------'
	like= lambda x:-likefunc2(PARTYPE(x[0],x[1]), estimator)
	#result=fmin(like, 10**(random.rand(2)-0.5), xtol=0.001, ftol=1e-4, maxiter=1000, maxfun=5000, full_output=True)
	x=fmin(like, x0, xtol=0.01, ftol=1e-3) #fmin works much better than minuit, probably due to the bad hessian calc in minuit
	print x
	print '---------------------------------'
	minuit_like(x0=x0)
	print '================================='

#for estimator in [4,8,10,11]:
  #print '------------- estimator ', estimator, '----------------'
  #like= lambda x:-likefunc2(PARTYPE(x[0],x[1]), estimator)
  ##result=fmin(like, 10**(random.rand(2)-0.5), xtol=0.001, ftol=1e-4, maxiter=1000, maxfun=5000, full_output=True)
  #x=fmin(like, [1,1], xtol=0.01, ftol=1e-3)
  #print x
  #print '================================='
 
