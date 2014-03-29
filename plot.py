import sys,os
from dynio import *
os.environ['OMP_NUM_THREADS']='32'

init()
m=-0.
freeze_energy(m,0.)

from matplotlib.pyplot import *
from numpy import *

x=m+linspace(-0.01,0.01,100)
y=array([like2(a,0.) for a in x])
#print zip(x,y)
plot(x,y)
xlabel(r'$\log(\rho_s/\rho_{s0})$')
ylabel(r'$\ln$(Like)')
#savefig('Like_iterative_step.eps') 
show()

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