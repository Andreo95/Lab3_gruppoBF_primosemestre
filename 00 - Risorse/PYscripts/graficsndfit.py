import pylab
import numpy
from scipy.optimize import curve_fit

DV,dV,I,dI = pylab.loadtxt('C:\\Users\\user\\Desktop\\laboratorio\\L2\\data12.txt', unpack = True)

pylab.figure()
pylab.plot(DV,I,'o',color='blue')
pylab.errorbar(DV,I,dV,dI, linestyle='',ecolor='black')

def f(x,a,b):
    return a*x+b

init=pylab.array([0,0])
    
popt,pcov= curve_fit(f,DV,I,init,dI)
a,b     =popt
da,db   =pylab.sqrt(pcov.diagonal())

print('a=',a,'b=',b,'da=',da,'db=',db,pcov)
x=numpy.linspace(0,5,500)
pylab.plot(x,f(x,a,b))

pylab.figure()
pylab.plot(DV,(I-f(DV,a,b))/dI)


Pchi=scipy.special.chdtrc(len(p)-2,chi_quadro)
print('Chi=',chi_quadro,'Pchi=',Pchi, 'dof=',len(Dq))