import os
import inspect
import scipy
from scipy.optimize import curve_fit
import pylab
import numpy as np


filename = inspect.getframeinfo(inspect.currentframe()).filename    #ricerca della directory dello script e impostazione come directory di lavoro
path = os.path.dirname(os.path.abspath(filename))
os.chdir(path) 

# indirizzo della cartella dati e cartela di salvataggio (completare con il nome della cartella i due indirizzi)
indirizzo_dati        = ''
indirizzo_salvataggio = ''

#nomi file (inserire i nomi dei file)
file_origine     = 'dati2.txt'
file_salvataggio = 'save.txt'

# importiamo i dati (actng!! il file deve contenere solo una serie di dati altrimenti raddoppiare)
sxa1,sya = pylab.loadtxt( r'%s%s' %(indirizzo_dati,file_origine), unpack = True ) 

sxa = sxa1*py.pi/180

# grafici dei dati iniziali
pylab.figure(1)
pylab.clf()

#grafico a (completare con il nome degliassi)
pylab.xlabel('tempo[us]')
pylab.ylabel('valori arduino[u.a.]')
pylab.grid(color = 'gray')
pylab.errorbar(sxa, sya, fmt = '.')
pylab.savefig('%sgrafico dati2.png' %(indirizzo_salvataggio))
pylab.show

# Funzione di fit (modificare la funzione e le variabili)
def f(x,a,b,c):
    return a*py.sin(2*x-2*b)**2/2+c

V =3  #numero di variabili (modificare il numero di variabili)
print('numero di variabili= %g' %(V))

val =[53.2,15*py.pi/180,0.09]

# Fit nell'intervallo selezionato (modificare gli spacchettamenti e le funzioni da printare)
popt_a, pcov_a = curve_fit(f, sxa, sya,p0=val,maxfev=6000)
a_fit,b_fit,c_fit= popt_a
da_fit,db_fit,dc_fit= pylab.sqrt(pcov_a.diagonal())
print('a = %f +- %f u.a.' % (a_fit,  da_fit))
print('b = %f +- %f u.a.' % (b_fit,  db_fit))
print('c = %f +- %f u.a.' % (c_fit,  dc_fit))
print(pcov_a)
print(np.corrcoef(pcov_a))

xx = pylab.linspace(min(sxa),max(sxa),200000,endpoint=True)
# Realizzazione e salvataggio del grafico (inserire nomi assi)
pylab.figure(2)
pylab.clf()
pylab.ylabel('corrente[uA]')
pylab.xlabel('angolo[rad]')
pylab.grid(color = 'gray')
pylab.errorbar(sxa, sya,[0.1]*len(sxa),[0.5*py.pi/180]*len(sxa), fmt = '.') 
pylab.plot(xx, f(xx,*popt_a), label='fit')
pylab.legend()
pylab.savefig('%sgrafico fit2.png' %(indirizzo_salvataggio))


#analisi errori
me=sum((f(sxa,*popt_a)-sya))/len(sxa)/100

#grafico degli errori (inserire nome asse x)
pylab.figure(3)
pylab.clf()
pylab.xlabel('angolo[rad]') 
pylab.ylabel('residui normalizzati')
pylab.grid(color = 'gray')
pylab.plot(sxa,(f(sxa,*popt_a)-sya)/100, '.', label='data')
pylab.plot(sxa,scipy.zeros(len(sxa)) , label='rif') 
pylab.plot(sxa,scipy.ones(len(sxa))*me, label='media')
pylab.legend()
pylab.savefig('%sgrafico errori2.png' %(indirizzo_salvataggio,))

#chiquadro 
za = sum((sya-f(sxa,*popt_a))**2)/0.01**2
errore_a = scipy.sqrt(za/(len(sxa)))
print('chiquadro a      = %g' %(za))
print('numero di dati A = %g' %(len(sya)))
print('error a          = %g' %(errore_a))

print(pylab.sqrt(c_fit))