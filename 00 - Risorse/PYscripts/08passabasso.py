import scipy
from scipy.optimize import curve_fit
import pylab
import numpy
import math

## indirizzo e nome file
indirizzo_dati    = '/afs/uz.sns.it/user/albord95/Scrivania/'
file_origine      = 'data08.txt'

# importiamo i dati (actng!! il file deve contenere solo una serie di dati altrimenti raddoppiare)
sxa,sya1,sya2,sfa = pylab.loadtxt( r'%s%s' %(indirizzo_dati,file_origine), unpack = True ) 
sya = sya2/sya1
dy = 0.02*sya

# Funzione di fit
def f(x, a):
    return 1/(scipy.sqrt(1+(x/a)**2))

# Fit nell'intervallo selezionato (modificare gli spacchettamenti e le funzioni da printare)
popt_a, pcov_a = curve_fit(f, sxa, sya)
a_fit = popt_a
da_fit = pylab.sqrt(pcov_a.diagonal())
print(' ')
print('f_T   = %f +- %f u.a.' % (a_fit,  da_fit))
#chiquadro 
chi2 = (((sya-f(sxa,*popt_a))/dy)**2).sum()
print('chi2  = %g' %(chi2))
print('dof   = %g' %(len(sya)-1))
print(' ')
print(pcov_a)

# grafico fit logaritmico
pylab.figure(2)
pylab.clf()
pylab.title('diagramma di Bode del guadagno (attenuazione)')
pylab.ylabel('attenuazione $V_{out} / V_{in}$ [dB]')
pylab.xscale('log')
pylab.yscale('log')
pylab.xlabel('frequenza[Hz]')
pylab.grid(color = 'gray')
pylab.errorbar(sxa, sya, dy, fmt = '.')
curva = numpy.logspace(math.log10(numpy.amin(sxa)), math.log10(numpy.amax(sxa)), 1000)
pylab.plot(curva, f(curva,*popt_a), label='fit')
pylab.legend()

#grafico degli errori
pylab.figure(3)
pylab.clf()
pylab.xlabel('frequenza[Hz]') 
pylab.xscale('log')
pylab.ylabel('residui normalizzati')
pylab.grid(color = 'gray')
pylab.plot(sxa,sya-(f(sxa,*popt_a)), '.', label='data')
pylab.plot(sxa,scipy.zeros(len(sxa)) , label='rif') 
media=sum((sya-f(sxa,*popt_a)))/len(sxa)
pylab.plot(sxa,scipy.ones(len(sxa))*media, label='media')
pylab.legend()

pylab.show()
