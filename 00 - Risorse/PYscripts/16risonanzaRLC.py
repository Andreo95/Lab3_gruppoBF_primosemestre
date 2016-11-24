import scipy
from scipy.optimize import curve_fit
import pylab
import numpy

## indirizzo e nome file
indirizzo_dati    = '/afs/uz.sns.it/user/albord95/Scrivania/'
file_origine      = 'data16.txt'

# importiamo i dati
sxa,dxa,sya1,dsya1,sya2,dsya2 = pylab.loadtxt( r'%s%s' %(indirizzo_dati,file_origine), unpack = True ) 
sya = sya1/sya2
dsya = sya*(dsya2/sya2+dsya1/sya1)

# Funzione di fit (modificare la funzione e le variabili)
def f(x,a,b,c):
    return 2*pylab.pi*a*x/pylab.sqrt((2*pylab.pi*x*b)**2+(1-(x/c)**2)**2)

val =[35.5*10**(-6),35.5*10**(-6),728]

# Fit nell'intervallo selezionato
popt_a, pcov_a = curve_fit(f, sxa, sya,p0=val,sigma=dsya,maxfev=6000)
a_fit,b_fit,c_fit= popt_a
da_fit,db_fit,dc_fit= pylab.sqrt(pcov_a.diagonal())
print('a = %g +- %g u.a.' % (a_fit,  da_fit))
print('b = %g +- %g u.a.' % (b_fit,  db_fit))
print('c = %g +- %g u.a.' % (c_fit,  dc_fit))
print(pcov_a)
print(numpy.corrcoef(pcov_a))
print('cov norm ab = %g' %(pcov_a[0,1]/(numpy.sqrt(pcov_a[0,0]*pcov_a[1,1]))))
print('cov norm ac = %g' %(pcov_a[0,2]/(numpy.sqrt(pcov_a[0,0]*pcov_a[2,2]))))
print('cov norm bc = %g' %(pcov_a[1,2]/(numpy.sqrt(pcov_a[1,1]*pcov_a[2,2]))))

xx = pylab.linspace(302,1418,2000,endpoint=True)

# Realizzazione e salvataggio del grafico (inserire nomi assi)
pylab.figure(2)
pylab.clf()
pylab.title('circuito risonante RLC')
pylab.ylabel('guadagno')
pylab.xlabel('frequenza[Hz]')
pylab.grid(color = 'gray')
pylab.errorbar(sxa, sya,dsya,dxa, fmt = '.') 
pylab.plot(xx, f(xx,*popt_a), label='fit')
pylab.legend()

#analisi errori
me=sum((f(sxa,*popt_a)-sya))/len(sxa)/100

#grafico degli errori (inserire nome asse x)
pylab.figure(3)
pylab.clf()
pylab.xlabel('tempo[ms]') 
pylab.ylabel('residui normalizzati')
pylab.grid(color = 'gray')
pylab.plot(sxa,(f(sxa,*popt_a)-sya)/dsya, '.', label='data')
pylab.plot(sxa,scipy.zeros(len(sxa)) , label='rif') 
pylab.plot(sxa,scipy.ones(len(sxa))*me, label='media')
pylab.legend()

pylab.show()

#chiquadro 
za = sum(((sya-f(sxa,*popt_a))**2)/dsya**2)
print('chi2  = %g' %(za))
print('dof   = %g' %(len(sya)-3))

R = 355
dR = 4
C =a_fit/R
dC =(da_fit/a_fit+dR/R)*a_fit/R
L = 1/c_fit**2/4/pylab.pi**2/C
dL = (2*dc_fit/c_fit+dC/C)*L
r = (b_fit/C)-R
dr = (db_fit/b_fit+dC/C)*(b_fit/C)+dR
print('C = %g +- %g' %(C,dC ))
print('L = %g +- %g' %(L,dL))
print('r = %g +- %g' %(r,dr))
print('da = %g'   %(da_fit))

HIV = f(c_fit,*popt_a)/2
#print(HIV)
'''
for i in range(2000):
    stop = HIV - f(xx[i], *popt_a)
    if (stop<0):
        print(xx[i])
    i=i+1
'''
