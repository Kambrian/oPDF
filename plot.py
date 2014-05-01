import sys,os
from dynio import *
os.environ['OMP_NUM_THREADS']='32'

init()

from matplotlib.pyplot import *
from numpy import *
from scipy.stats import chi2,norm

def chi2sig(ts,dof=1):
  '''convert a chi2 TS to sigma'''
  pval=chi2.sf(ts,dof)
  sigma=norm.ppf(1.0-pval/2)
  return sigma,pval

x0=linspace(0.8,1.2,10)

def plot_significance(id=0,dim=0,x=x0):
  '''plot the significance on the subsample specified by id'''
  select_particles(id)
  #x=linspace(0.8,1.2,20)
  if dim==0:
    y=array([like2(a,1.) for a in x])
  else:
    y=array([like2(1.,a) for a in x])
  y,py=chi2sig(-y,1)
  #print zip(x,y)
  plot(x,y)
  if dim==0:
    xlabel(r'$\log(\rho_s/\rho_{s0})$')
  else:
    xlabel(r'$\log(\r_s/\r_{s0})$')
  ylabel(r'Discrepancy Level($\sigma$)')
  #plot(xlim(),[-2.48,-2.48],'--')
  #plot(xlim(),[-2.28,-2.28],'--')
  ylim(0,3)
  #title('RADIAL_BIN_ESTIMATOR')
  #savefig('ORBITAL_ROULETTE_ESTIMATOR.eps') 
  #show()
  return y

ion()
figure()
sid=range(100)
y=[]
for i in sid:
  select_particles(i)
  s=like2(1.,1.)
  y.append(s)
  print i,s
plot(sid,y)  

figure()  
xmin=[]
for i in range(0,100,5):
  y=plot_significance(i)
  print i,x0[y.argmin()],y.min()
  xmin.append(x0[y.argmin()])
  draw()

  
"""
figure()
x0=linspace(-0.5,1,10)
for i,m in enumerate(x0):
	freeze_energy(m,0.)
	x=m+arange(-0.5,1,0.01)
	y=array([like(a,0.) for a in x])
	plot(x,y)
	hold('on')

x0=linspace(-0.01,0.01,50)
y0=zeros_like(x0)
for i,m in enumerate(x0):
	freeze_energy(m,0.)
	y0[i]=(like(m,0.)+0.5*like(m-0.01,0.)+0.5*like(m+0.01,0.))/2.
plot(x0,y0,'k-')
xlabel('log10(Rhos/Rhos0)')
ylabel('ln(Like)')
#savefig('Like_iterative_mixed.eps')
show()
"""
