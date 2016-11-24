import scipy
from scipy.optimize import curve_fit
import pylab
import numpy

## indirizzo e nome file
indirizzo_dati    = '/afs/uz.sns.it/user/albord95/Scrivania/'
file_origine      = 'data05carica.txt'

## importiamo i dati
sxa, sya = pylab.loadtxt( r'%s%s' %(indirizzo_dati,file_origine), unpack = True ) 

## grafici dei dati iniziali
pylab.figure(1)
pylab.clf()
pylab.xlabel('tempo[$\mu$s]')
pylab.ylabel('valori arduino[u.a.]')
pylab.grid(color = 'gray')
pylab.errorbar(sxa, sya, 1, fmt = '.')

## Funzione di fit
def f(x, a, b):
    return a*(1-scipy.exp(-x/b))

## best-fit
popt_a, pcov_a = curve_fit(f, sxa, sya)
a_fit,   b_fit = popt_a
da_fit, db_fit = pylab.sqrt(pcov_a.diagonal())
print(' ')
print('V_0 = %f +- %f u.a.' % (a_fit,  da_fit))
print('tau = %f +- %f micro s'      % (b_fit,  db_fit))
#chiquadro e coefficiente di correlazione
chi2 = ((sya-f(sxa,*popt_a))**2).sum()
cov_norm = pcov_a[0,1]/(numpy.sqrt(pcov_a[0,0]*pcov_a[1,1]))
print('chiquadro  = %g' %(chi2))
print('dof        = %g' %(len(sxa)-2))
print('cov norm   = %g' %(cov_norm))
print(' ')
print(pcov_a)
print(numpy.corrcoef(pcov_a))

## grafico fit
pylab.figure(2)
pylab.clf()
pylab.title('carica condensatore')
pylab.ylabel('valore arduino[u.a.]')
pylab.xlabel('tempo[$\mu$s]')
pylab.grid(color = 'gray')
pylab.errorbar(sxa, sya, 1, fmt = '.') 
pylab.plot(sxa, f(sxa,*popt_a), label='fit')
pylab.legend()

## grafico degli errori
pylab.figure(3)
pylab.clf()
pylab.xlabel('tempo[$\mu$s]')
pylab.title('carica condensatore') 
pylab.ylabel('residui normalizzati')
pylab.grid(color = 'gray')
pylab.plot(sxa,sya-(f(sxa,*popt_a)), '.', label='data')
pylab.plot(sxa,scipy.zeros(len(sxa)) , label='rif') 
media=sum((sya-f(sxa,*popt_a)))/len(sxa) #calcolo media residui
pylab.plot(sxa,scipy.ones(len(sxa))*media, label='media')
pylab.legend()

pylab.show()
