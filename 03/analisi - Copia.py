#o fai delle funzioni usabili o questi mostri ricopiati sono proprio inguardabili...se vogliamo uno script che sia uno scrip
import sys, os
sys.path.append( os.path.join(os.path.realpath('..'), '00 - Risorse', 'Python') )

import numpy as np, pylab, matplotlib.pyplot as plt, matplotlib as mpl, scipy.stats.distributions as dists
from lab import *
from PyTeX import *

folder = os.path.realpath('.')

subplotGS = mpl.gridspec.GridSpec(5,1)

plt.close('all')
greencol = (0,1,0)
flatline = np.vectorize(lambda x, c: c)


def cross(m, q, c):
    return  (c-q)/m

def dcross(linepars, linecov, c, dc):
    '''errore sul intesezione retta costante
    linepars:   parametri retta linepars[0]*x + linepars[1]
    linecov:    covaianza su di questo
    c: costante
    dc: varianza c
    '''
    d=[-(c-linepars[1])/linepars[0]**2, -1/linepars[0]]
    return linecov[0][0]*d[0]**2+linecov[1,1]*d[1]**2+2*linecov[1, 0]*d[0]*d[1]+dc*d[0]**2

def FindCross(linepars, linecov, c, dc):
    '''intesezione retta costante
    
    linepars : parametri retta linepars[0]*x + linepars[1]
    linecov:    covaianza su di questo
    c: costante
    dc: varianza c
    Returns
    ---------
    Coppia (valore , varianza) 
    '''
    return np.array([cross(linepars[0], linepars[1], c), dcross(linepars, linecov, c, dc)])

#### Parte 05c

datafile = '{0}{1}_dati.txt'.format(5,'c')

rawdata = np.loadtxt( os.path.join(folder, 'Dati', datafile) ).T


Vbe, Vce, Vrb = rawdata
RL    = 981		#ohm 
RB    = 46.2e3	#kohm
C1    = 10.55e-9	#nF
Vcc   = 10.46	#V
dVcc =  0.04	#V

dVbe, dVce = mme([Vbe, Vce], 'volt', 'oscil')
dVrb = mme(Vrb, 'volt')
dRB, dRL = mme([RB, RL], 'ohm')

Ib = Vrb/RB
Ic = (Vcc-Vce)/RL

dIb = Ib * np.sqrt( (dVrb/Vrb)**2 )
dIc = Ic * np.sqrt( (dVcc**2)/(Vcc-Vce)**2 )




XX  = Vbe
YY  = Ic
dXX = dVbe
dYY = dIc


points = np.linspace( np.amin(XX) - (np.amax(XX)-np.amin(XX))/20 ,np.amax(XX) + (np.amax(XX)-np.amin(XX))/20,num=max(len(XX)*10, 200) )

rescaleX = 1e0
rescaleY = 1e3
plt.figure(1)
plt.title('SONO TITO')
plt.ylabel('Vbe [V]')
plt.xlabel('Ic [mA]')
plt.errorbar(XX*rescaleX, YY*rescaleY, dYY*rescaleY, dXX*rescaleX, fmt='none', ecolor='black', label='data')


f = lambda x, b, a: a+b*x
df = lambda x, a, b: b

line = f

XX  = Ib
YY  = Ic
dXX = dIb
dYY = dIc

linedata = XX<55e-6

pars, pcov = fit_generic_xyerr(f, df, XX[linedata], YY[linedata], dXX[linedata], dYY[linedata], p0=None)

resd = (YY[linedata] - f(XX[linedata], *pars)) /dYY[linedata]
print( 'ChiSquare = {0} ({1} DoF)'.format(np.sum(resd**2), len(XX[linedata])-len(pars)) )
print(pars, pcov)

linepars, linecov = pars, pcov

flatdata = XX>70e-6

flatpar, flatvar = fit_const_yerr( YY[flatdata], dYY[flatdata])

resd = (YY[flatdata] - flatline(XX[flatdata], flatpar)) /dYY[flatdata]
print( 'ChiSquare = {0} ({1} DoF)'.format(np.sum(resd**2), len(XX[flatdata])-1) )
print(flatpar, flatvar)


points = np.linspace( np.amin(XX) - (np.amax(XX)-np.amin(XX))/20 ,np.amax(XX) + (np.amax(XX)-np.amin(XX))/20,num=max(len(XX)*10, 200) )

rescaleX = 1e6
rescaleY = 1e3
plt.figure(2)
plt.title('SONO TITO')
plt.ylabel('Ib [$\mu$A]')
plt.xlabel('Ic [mA]')
plt.errorbar(XX*rescaleX, YY*rescaleY, dYY*rescaleY, dXX*rescaleX, fmt='none', ecolor='black', label='data')
plt.plot(points[points<80e-6]*rescaleX, line(points[points<80e-6], *linepars)*rescaleY, color=greencol, label='regime attivo')
plt.plot(points*rescaleX, flatline(points, flatpar)*rescaleY, color='r', label='saturazione')
plt.legend()



#### Parte 05e

rawdata1=rawdata2=np.asarray([[], []])



datafile = '{0}_{1}_dati.txt'.format('5e', 2)
rawdata2 = np.loadtxt( os.path.join(folder, 'Dati', datafile) ).T

rawdata = np.append(rawdata1, rawdata2, axis=1)

Vcc, Vce = rawdata

RL    = 981     #ohm 
RB    = 46.2e3	#kohm
C1    = 10.55e-9	#nF
Vrb  = 1.564	#V

Ib = Vrb/RB
Ic = (Vcc-Vce)/RL
dVce, dVcc = mme([Vce,Vcc], 'volt', 'oscil')
dIc = Ic * np.sqrt( (dVce**2 )/(Vcc-Vce)**2 )


XX  = Vce
YY  = Ic
dXX = dVce
dYY = dIc



f = line
df = df

pars, pcov = fit_generic_xyerr(f, df, XX, YY, dXX, dYY, p0=None)

linepars = pars
linecov = pcov

print('\n\nVearly = {} \pm {}'.format(*FindCross(linepars, linecov, 0,0)))

print(linepars, linecov)

resd = (YY - f(XX, *pars)) /dYY
print( 'ChiSquare = {0} ({1} DoF)\n'.format(np.sum(resd**2), len(XX)-len(pars)) )

points = np.linspace( np.amin(XX) - (np.amax(XX)-np.amin(XX))/20 , np.amax(XX) + (np.amax(XX)-np.amin(XX))/20, num=max(len(XX)*10, 200) )



rescaleX = 1e0
rescaleY = 1e3

plt.figure(3)
plt.clf()
plt.title("Titolo")
plt.ylabel('Ic [mA]')
plt.xlabel("Vce [V]")
plt.errorbar(XX*rescaleX, YY*rescaleY, dYY*rescaleY, dXX*rescaleX, fmt='none', ecolor='black', label='data')
plt.plot(points*rescaleX, line(points, *linepars)*rescaleY, color=greencol, label='fit')
plt.ylim((3, 7))
plt.legend()


#### Parte 06



R1 = 15.03e3	#kOhm
R2 = 98.4e3		#kOhm
RL = 2.26e3		#kOhm



XX  = [1]
YY  = [1]
dXX = [1]
dYY = [1]


f = lambda x, params: params
df = lambda x, params: params

pars, pcov = fit_generic_xyerr(f, df, XX, YY, dXX, dYY, p0=None)


fignum = 000
titolo = 'SONO TITO'
labelx = 'SONO LA X'
labely = 'SONO LA Y'
rescaleX = 1e0
rescaleY = 1e0

resd = (YY - f(XX, *pars)) /dYY
DoFs = len(XX)-len(pars)
print( 'ChiSquare = {0} ({1} DoF, p = {2})'.format( np.sum(resd**2), DoFs, dists.chi2.sf(np.sum(resd**2), DoFs) ) )


plt.show()