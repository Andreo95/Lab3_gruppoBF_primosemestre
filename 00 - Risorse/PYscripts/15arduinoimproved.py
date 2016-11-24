import scipy
from scipy.optimize import curve_fit
import pylab
import numpy

cap = 0.47*10**(-6)

## indirizzo e nome file
indirizzo_dati    = '/Scrivania/'
file_origine      = 'data15.txt'

# importiamo i dati
sxa,sya = pylab.loadtxt( r'%s%s' %(indirizzo_dati,file_origine), unpack = True ) 

# Funzione di fit
def f(x,a,b,c,d,e):
    return (a*(pylab.exp(-1/b*x)*pylab.cos(2*pylab.pi/c*x+d))+e)

val =[500,20000,2000,0,500]

# Fit nell'intervallo selezionato (modificare gli spacchettamenti e le funzioni da printare)
popt_a, pcov_a = curve_fit(f, sxa, sya,p0=val,sigma=1,maxfev=6000)
b_fit,c_fit= popt_a[1:3]
db_fit,dc_fit= pylab.sqrt(pcov_a.diagonal()[1:3])
print('')
print('tau = %f +- %f micro s' % (b_fit,  db_fit))
print('T   = %f +- %f micro s' % (c_fit,  dc_fit))
print(pcov_a)
print(numpy.corrcoef(pcov_a))
print('')

# grafico fit
pylab.figure(2)
pylab.clf()
pylab.ylabel('tensione [u.a.]')
pylab.xlabel('tempo[$\mu$s]')
pylab.grid(color = 'gray')
pylab.errorbar(sxa, sya, 1, fmt = '.')
curva = numpy.linspace(numpy.amin(sxa), numpy.amax(sxa),2000) 
pylab.plot(curva, f(curva,*popt_a), label='fit')
pylab.legend()

#analisi errori
me=sum(-(f(sxa,*popt_a)-sya))/len(sxa)

#grafico degli errori (inserire nome asse x)
pylab.figure(3)
pylab.clf()
pylab.xlabel('tempo[$\mu$s]') 
pylab.ylabel('residui normalizzati')
pylab.grid(color = 'gray')
pylab.plot(sxa,-(f(sxa,*popt_a)-sya), '.', label='data')
pylab.plot(sxa,scipy.zeros(len(sxa)) , label='rif') 
pylab.plot(sxa,scipy.ones(len(sxa))*me, label='media')
pylab.legend()

pylab.show()

#chiquadro 
za = sum((sya-f(sxa,*popt_a))**2)
print('chi2  = %g' %(za))
print('dof   = %g' %(len(sya)-5))

print('omega = %g +- %g [Hz]' %(2*pylab.pi/c_fit*10**6,2*pylab.pi*10**6*dc_fit/c_fit**2))
ind = (c_fit/10**6/(2*pylab.pi))**2/cap
print('L = %g[H] ' %(ind))
print('R = %g [OHM}' %(2*ind/(b_fit*10**(-6))))
