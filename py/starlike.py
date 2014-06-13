from ctypes import *
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as pl
import numpy as np
from matplotlib.patches import Ellipse

star=CDLL("./likelihood.so")
prototype=CFUNCTYPE(c_double, c_double, c_double, c_double, c_double, c_double, c_double)
likelihood=prototype(("likelihood",star))

like=lambda beta,mass,c200,gamma1,gamma2,CC: -likelihood(beta,mass,c200,gamma1,gamma2,CC)


#print like(0.696837,2.31556,12.7983,1.8534,9.444,55.72)
#print like(0.7,2.31556,12.7983,1.8534,8.24,55.)
#print like(0.7,1.84,16.1,1.82,8.24,55.)

#print like(0.7,14.,7.227,1.82,8.24,55.)

from iminuit import Minuit
from iminuit.ConsoleFrontend import ConsoleFrontend

m=Minuit(like,beta=0.715,error_beta=0.1, limit_beta=(-0.5,1.),mass=1.842,error_mass=0.1, limit_mass=(0.5,3.),c200=16.2,error_c200=0.1, limit_c200=(3.,30.),gamma1=2.82,error_gamma1=0.1, limit_gamma1=(0.,5.),gamma2=8.24,error_gamma2=0.1, limit_gamma2=(5.,15.),CC=269.,error_CC=0.1, limit_CC=(7.,400.),print_level=2,errordef=0.5, frontend=ConsoleFrontend())
m.print_param()
m.tol=10.
m.migrad();
#m.minos();

c200=m.fitarg['c200']
sigmac200=m.errors['c200']
mass=m.fitarg['mass']
sigmamass=m.errors['mass']

AC=(np.log(1.+c200)-c200/(1.+c200))
rou0=np.log10(c200**3*200.*147.692/3./AC)
sigmarou02=200.*147.692/3.*(3.*c200**2*AC-c200**4/(1.+c200)**2)/AC**2*sigmac200
rou02=10.**rou0

rs=mass*10.**12/4./np.pi/rou02/AC
rs=rs**(1./3.)
sigmars=1./3.*rs**(-2)*np.sqrt((1./4./np.pi/rou02/AC)**2*(sigmamass*10.**12)**2+(rs**3/rou02)**2*sigmarou02**2+(rs**3/AC*c200/(1.+c200)**2)**2*sigmac200**2)

R200=c200*rs
sigmaR200=np.sqrt(rs**2*sigmac200**2+c200**2*sigmars**2)

print 'R200:',R200,sigmaR200
print 'rs:',rs,sigmars
print 'rou0:',rou0,1./rou02/np.log(10.)*sigmarou02

cov=m.matrix()
print cov
#minv=linalg.inv(conv[1:3,1:3])
#print minv
cov2=np.array(cov)
cov3=cov2[1:3,1:3]

fig = pl.figure(figsize=(7, 7))
x=1.842
y=16.19
pl.plot(x,y,'or')
pl.plot(mass,c200,'ok')
pl.plot([mass,mass],[c200-sigmac200,c200+sigmac200],'-k')
pl.plot([mass-sigmamass,mass+sigmamass],[c200,c200],'-k')

ll, v = np.linalg.eig(cov3)
ll= np.sqrt(ll)
print ll
print v
print cov3
print np.rad2deg(np.arccos(v[0, 0]))
ax=pl.gca()
for i in np.arange(1,3,1):
    ell = Ellipse(xy=(mass, c200),width=ll[0]*2.*i*2.25, height=ll[1]*2.*i*2.25,angle=np.rad2deg(np.arccos(v[0, 0])))
    ell.set_facecolor('none')
    ax.add_artist(ell)

pl.ylim(10.,22.)
pl.xlim(1.2,2.3)
pl.xlabel(r'$\mathrm{M_{200}[10^{12}M_\odot]}$',fontsize=15)
pl.ylabel(r'$c_{200}$',fontsize=15)
pl.savefig('contour.ps', dpi=150, format='ps')

pl.show()





