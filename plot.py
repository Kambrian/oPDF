import sys,os
from dynio import *
os.environ['OMP_NUM_THREADS']='32'

init()

from matplotlib.pyplot import *
from numpy import *

x=linspace(-0.02,0.02,20)
y=array([like(a,0.) for a in x])
plot(x,y)
show()
