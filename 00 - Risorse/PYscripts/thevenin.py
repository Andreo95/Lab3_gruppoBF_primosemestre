import pylab
import numpy
from scipy.optimize import curve_fit

R,dR,I,dI,ddp = pylab.loadtxt('C:\\Users\\user\\Desktop\\laboratorio\\L2\\2\\dataTH.txt', unpack = True)

ra=ddp/I
Y=1/(I)
dY=dI/(I**2)
X=R
dX=dR

def rpar(r1,r2):
    return r1*r2/(r1+r2)
#resistenza di thevenin    
def rth(rl,vth,i,ra):
    return vth/i-rl-ra
#errore associato
def drth(drl,dvth,di,vth,i):
    return (di*vth/i**2+drl+dvth/i)
    
def f(X,a,b):
    return a*X+b
#mi aspetto che la differenza del potenziale sia sui 5 volt e la resistenza interna sui 20 ohm, scrivo in kilo-ohm

init=pylab.array([1/5,0.003])
    
popt,pcov= curve_fit(f,X,Y,init,dY)
a,b     =popt
da,db   =pylab.sqrt(pcov.diagonal())

print('a=',a,'+-',da,'b=',b,'+-',db,'MCov=',pcov)

V=1/a
dV=da/(a**2)
rg=b/a
drg=((db/b)+(da/a))*rg
print('rg=',rg,'+-',drg,'V(dal fit)=',V,'+-',dV)

x=numpy.logspace(-2,4,500)
pylab.title('dati e fit numerico')
pylab.xlabel('(R+ra) [kohm]')
pylab.ylabel('1/I [1/mA]')
pylab.xscale('log')
pylab.yscale('log')
pylab.plot(x,f(x,a,b))
pylab.plot(X,Y,'o')
pylab.errorbar(X, Y, dY, linestyle='',ecolor='black')

pylab.figure()
pylab.plot(X,(Y-f(X,a,b))/dY)

chi_quadro=0.
for t in range(len(I)):
    chi_quadro+=((Y[t]-f(X[t], a, b))/dY[t])**2

print('Chi=',chi_quadro, 'dof=',len(I))
