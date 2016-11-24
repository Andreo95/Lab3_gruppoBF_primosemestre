import scipy
from scipy.optimize import curve_fit
import pylab
import numpy
import math

## indirizzo e nome file
indirizzo_dati    = '/afs/uz.sns.it/user/albord95/Scrivania/'
file_origine      = 'data09.txt'

# importiamo i dati 
sya1,sxa = pylab.loadtxt( r'%s%s' %(indirizzo_dati,file_origine), unpack = True ) 
sya = ((sya1-sxa)*4.92/(3.28))*10**(-3) #MODIFICARE  mA
sxa = sxa*4.92*10**(-3) #MODIFICARE
dy  = (2*4.92/3.28)*10**(-3)  #mA
dx  = 4.92*10**(-3)

# Funzione di fit (modificare la funzione e le variabili)
def f(x, a,b):
    return a*(pylab.exp((1/b)*x)-1)

# best-fit
val =[10,1/(0.026*1.5)]
popt_a, pcov_a = curve_fit(f, sxa, sya,p0=val)
a_fit , b_fit= popt_a
da_fit ,db_fit= pylab.sqrt(pcov_a.diagonal())
print(' ')
print('I_0    = %f +- %f micro A' % (a_fit*10**6,  da_fit*10**6))
print('etaV_T = %f +- %f mV' % (b_fit*1000,  db_fit*1000))
#chiquadro
sigma = math.sqrt(dy**2)
chi2 = (((sya-f(sxa,*popt_a))/sigma)**2).sum()
print('chi2   = %g' %(chi2))
print('dof    = %g' %(len(sya)-2))
cov_norm = pcov_a[0,1]/(numpy.sqrt(pcov_a[0,0]*pcov_a[1,1]))
print('cov norm   = %g' %(cov_norm))
print(' ')
print(pcov_a)

# Realizzazione e salvataggio del grafico (inserire nomi assi)
pylab.figure(2)
pylab.clf()
pylab.title('curva caratteristica del diodo')
pylab.ylabel('corrente [$\mu$A]')
pylab.xlabel('tensione[V]')
pylab.grid(color = 'gray')
pylab.errorbar(sxa, sya, dy, dx, fmt = '.') 
pylab.plot(sxa, f(sxa,*popt_a), label='fit')
pylab.legend()

#grafico degli errori (inserire nome asse x)
pylab.figure(3)
pylab.clf()
pylab.xlabel('tensione[mV]') 
pylab.ylabel('residui normalizzati')
pylab.grid(color = 'gray')
pylab.plot(sxa,(sya-f(sxa,*popt_a)), '.', label='data')
pylab.plot(sxa,scipy.zeros(len(sxa)) , label='rif')
media=sum((sya-f(sxa,*popt_a)))/len(sxa)
pylab.plot(sxa,scipy.ones(len(sxa))*media, label='media')
pylab.legend()

pylab.show()
