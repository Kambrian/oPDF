cd py
from pylab import *
ion()
from dynio import *
lib.open()
FullSample=Tracer('Mock')
Sample=FullSample.copy(0,1000)
lib.choose_integral_routines(2, 0, 0)
x=Sample.fmin_jointLE(estimator=16, nbinL=100,nbinE=1)
print x
lib.MODEL_TOL_REL.value=1e-4
x=Sample.fmin_jointLE(estimator=16, nbinL=100,nbinE=1)
print x
lib.MODEL_TOL_REL.value=1e-2
x=Sample.fmin_jointLE(estimator=16, nbinL=100,nbinE=1)
print x

bin=linspace(-0.05,0.05,100)+x[0][0]
def dprof(linestyle):
  c=x[0][1]
  y=[Sample.jointLE_Flike([m, c], estimator=16, nbinL=100, nbinE=1) for m in bin]
  plot(bin, np.array(y)-min(y), linestyle)

figure()
lib.choose_integral_routines(2, 0, 0)#best. 
dprof('rx')
lib.MODEL_TOL_REL.value=1e-4
dprof('go')
lib.MODEL_TOL_REL.value=1e-5
dprof('b')

#figure()
#lib.choose_integral_routines(1, 1, 1) #bad
#dprof('r')
#lib.choose_integral_routines(1, 0, 1) #seems best? innerbin
#dprof('g')
#lib.choose_integral_routines(1, 1, 0) #not good
#dprof('b')
#lib.choose_integral_routines(1, 0, 0) #ok
#dprof('cx')
#lib.choose_integral_routines(1, 0, 0) #ok
#lib.MODEL_TOL_BIN.value=1e-6
#lib.MODEL_TOL_BIN_ABS.value=1e-8
#lib.MODEL_TOL_REL.value=1e-7
#dprof('k')
#legend(('11','01','10','00','000'))

#figure()
#lib.choose_integral_routines(2, 1, 1)#wrong
#dprof('r')
#lib.choose_integral_routines(2, 0, 1)#can be good with higher bin accuracy
#dprof('g')
#lib.choose_integral_routines(2, 1, 0)#good
#dprof('b--')
#lib.choose_integral_routines(2, 0, 0)#good
#dprof('cx')
#lib.choose_integral_routines(2, 0, 0)#best. 
#lib.MODEL_TOL_BIN.value=1e-6
#lib.MODEL_TOL_BIN_ABS.value=1e-8
#lib.MODEL_TOL_REL.value=1e-6
#dprof('k')
#legend(('11','01','10','00','000')) 

