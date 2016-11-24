import pylab
import numpy
from scipy.optimize import curve_fit

avg=pylab.array([78.1,144.0,206.3,283.7,365.1,449.2,532.5,576.0,663.2,788.5,860.9]) #average values for arduino measures
V=pylab.array([0.396,0.716,1.020,1.398,1.795,2.20,2.61,2.82,3.24,3.85,4.2]) #tester corresponding potential measures
dV=pylab.array([0.001,0.015,0.005,0.006,0.011,0.02,0.02,0.02,0.02,0.03,0.03])

def f(X,a,b):
    return a*X+b
    
init=pylab.array([0.005,0.01]) #expected values for calibration coefficient (5/1023) and offset parameter (0)

popt,pcov= curve_fit(f,avg,V,init,dV)
a,b     =popt
da,db   =pylab.sqrt(pcov.diagonal())

print('a=',a,'+-',da,'b=',b,'+-',db,'MCov=',pcov)

x=numpy.linspace(0.0,1023.0,500)
pylab.title('controllo della linearità di calibrazione')
pylab.xlabel('media di arduino [digit]')
pylab.ylabel('ddp [V]')
pylab.plot(x,f(x,a,b))
pylab.plot(avg,V,'o')
pylab.errorbar(avg, V, dV, linestyle='',ecolor='black')

pylab.figure()
pylab.title('residui normalizzati')
pylab.ylabel('(ddp-fit)/dV')
pylab.xlabel('digit arduino')
pylab.plot(avg,(V-f(avg,a,b))/dV)

chi_quadro=0.
for t in range(len(avg)):
    chi_quadro+=((V[t]-f(avg[t], a, b))/dV[t])**2

print('Chi=',chi_quadro, 'dof=',len(avg))
