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





fignum = 1
titolo = '$V_{be}$ vs $I_c$'
labelx = 'Vbe [V]'
labely = 'Ic [mA]'
xtype = 'linear'
ytype = 'linear'
rescaleX = 1e0
rescaleY = 1e3
whatplot = ('data','fit') # 'data', 'fit', 'resid', 'line', 'flatline'


f=lambda x, i_0, v_0: i_0*np.exp(x/v_0)
df=lambda x, i_0, v_0: f(x, i_0, v_0)/v_0
AA, BB= fit_generic_xyerr(f, df, XX[YY<9e3], YY[YY<9e3], dXX[YY<9e3], dYY[YY<9e3])
print(AA)
print(BB)
pars=AA


if xtype == 'linear':
	points = np.linspace( np.amin(XX) - (np.amax(XX)-np.amin(XX))/20 ,
	                      np.amax(XX) + (np.amax(XX)-np.amin(XX))/20,
	                      num=max(len(XX)*10, 200) )
elif xtype == 'log':
    points = np.logspace( np.log10( np.amin(XX)/(np.amax(XX)/np.amin(XX))**(1/20) ),
                          np.log10( np.amax(XX)*(np.amax(XX)/np.amin(XX))**(1/20) ),
                          num=max(len(XX)*10, 200) )

if whatplot:
	plt.figure(fignum)
	plt.clf()
	if 'resid' in whatplot:
		plt.subplot(subplotGS[:4])
	plt.title(titolo)
	plt.ylabel(labely)
	if 'resid' not in whatplot:
		plt.xlabel(labelx)
	plt.yscale(ytype)
	plt.xscale(xtype)
	if 'data' in whatplot:
		plt.errorbar(XX*rescaleX, YY*rescaleY, dYY*rescaleY, dXX*rescaleX, fmt='none', ecolor='black', label='data')
	if 'fit' in whatplot:
		plt.plot(points*rescaleX, f(points, *pars)*rescaleY, color='red', label='fit')
	if 'line' in whatplot:
		plt.plot(points[points<80]*rescaleX, line(points[points<80], *linepars)*rescaleY, color=greencol)
	if 'flatline' in whatplot:
		plt.plot(points*rescaleX, flatline(points, flatpar)*rescaleY, color=greencol)
    plt.xlim((0, 1))
	plt.legend()
	if 'resid' in whatplot:
		plt.subplot(subplotGS[4:])
		plt.ylabel('Norm. res.')
		plt.xlabel(labelx)
		plt.xscale(xtype)
		plt.plot(XX*rescaleX, resd , 'ob')
		plt.axhline(y=0, color='k')
		

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



fignum = 2
titolo = ''
labelx = 'Ib [$\mu$A]'
labely = 'Ic [mA]'
xtype = 'linear'
ytype = 'linear'
rescaleX = 1e6
rescaleY = 1e3
whatplot = ('data', 'line', 'flatline') # 'data', 'fit', 'resid', 'line', 'flatline'

if xtype == 'linear':
	points = np.linspace( np.amin(XX) - (np.amax(XX)-np.amin(XX))/20 ,
	                      np.amax(XX) + (np.amax(XX)-np.amin(XX))/20,
	                      num=max(len(XX)*10, 200) )
elif xtype == 'log':
    points = np.logspace( np.log10( np.amin(XX)/(np.amax(XX)/np.amin(XX))**(1/20) ),
                          np.log10( np.amax(XX)*(np.amax(XX)/np.amin(XX))**(1/20) ),
                          num=max(len(XX)*10, 200) )
if whatplot:
    plt.figure(fignum)
    plt.clf()
    if 'resid' in whatplot:
        plt.subplot(subplotGS[:4])
    plt.title(titolo)
    plt.ylabel(labely)
    if 'resid' not in whatplot:
        plt.xlabel(labelx)
    plt.yscale(ytype)
    plt.xscale(xtype)
    if 'data' in whatplot:
		plt.errorbar(XX*rescaleX, YY*rescaleY, dYY*rescaleY, dXX*rescaleX, fmt='none', ecolor='black', label='data')
	if 'fit' in whatplot:
		plt.plot(points*rescaleX, f(points, *pars)*rescaleY, color='red', label='fit')
	if 'line' in whatplot:
		plt.plot(points[points<80e-6]*rescaleX, line(points[points<80e-6], *linepars)*rescaleY, color=greencol, label='regime attivo')
	if 'flatline' in whatplot:
		plt.plot(points*rescaleX, flatline(points, flatpar)*rescaleY, color='r', label='saturazione')
	plt.legend()
    if 'resid' in whatplot:
        plt.subplot(subplotGS[4:])
        plt.ylabel('Norm. res.')
        plt.xlabel(labelx)
        plt.xscale(xtype)
        plt.plot(XX*rescaleX, resd , 'ob')
        plt.axhline(y=0, color='k')

# plt.savefig( os.path.join(folder, 'Figs-Tabs', 'fig_{}.png'.format(fignum)), bbox_inches='tight', dpi = 1000) # format='ps', papertype='a4', orientation='landscape'



#### Parte 05e

rawdata1=rawdata2=np.asarray([[], []])

# datafile = '{0}_{1}_dati.txt'.format('5e', 1)
# rawdata1 = np.loadtxt( os.path.join(folder, 'Dati', datafile) ).T

datafile = '{0}_{1}_dati.txt'.format('5e', 2)
rawdata2 = np.loadtxt( os.path.join(folder, 'Dati', datafile) ).T

rawdata = np.append(rawdata1, rawdata2, axis=1)

Vcc, Vce = rawdata

RL   = 981		#ohm 
RB   = 46.2e3	#kohm
C1   = 10.55e-9	#nF
Vrb  = 1.564	#V

dRB = mme(RB, 'ohm')
dVrb = mme(Vrb, 'volt')
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

fignum = 3
titolo = ''
labelx = 'Vce [V]'
labely = 'Ic [mA]'
xtype = 'linear'
ytype = 'linear'
rescaleX = 1e0
rescaleY = 1e3
whatplot = ('data', 'line') # 'data', 'fit', 'resid', 'line', 'flatline'


resd = (YY - f(XX, *pars)) /dYY
print( 'ChiSquare = {0} ({1} DoF)\n'.format(np.sum(resd**2), len(XX)-len(pars)) )

if xtype == 'linear':
    points = np.linspace( np.amin(XX) - (np.amax(XX)-np.amin(XX))/20 ,
                          np.amax(XX) + (np.amax(XX)-np.amin(XX))/20,
                          num=max(len(XX)*10, 200) )
elif xtype == 'log':
    points = np.logspace( np.log10( np.amin(XX)/(np.amax(XX)/np.amin(XX))**(1/20) ),
                          np.log10( np.amax(XX)*(np.amax(XX)/np.amin(XX))**(1/20) ),
                          num=max(len(XX)*10, 200) )
                          
if any(whatplot):
    plt.figure(fignum)
    plt.clf()
    if 'resid' in whatplot:
        plt.subplot(subplotGS[:4])
    plt.title(titolo)
    plt.ylabel(labely)
    if 'resid' not in whatplot:
        plt.xlabel(labelx)
    plt.yscale(ytype)
    plt.xscale(xtype)
    if 'data' in whatplot:
        plt.errorbar(XX*rescaleX, YY*rescaleY, dYY*rescaleY, dXX*rescaleX, fmt='none', ecolor='black', label='data')
    if 'fit' in whatplot:
        plt.plot(points*rescaleX, f(points, *pars)*rescaleY, color='red', label='fit')
    if 'line' in whatplot:
        plt.plot(points*rescaleX, line(points, *linepars)*rescaleY, color=greencol, label='fit')
    if 'flatline' in whatplot:
        plt.plot(points*rescaleX, flatline(points, flatpars)*rescaleY, color=greencol)
    plt.ylim((3, 7))
    plt.legend()
    if 'resid' in whatplot:
        plt.subplot(subplotGS[4:])
        plt.ylabel('Norm. res.')
        plt.xlabel(labelx)
        plt.xscale(xtype)
        plt.plot(XX*rescaleX, resd , 'ob')
        plt.axhline(y=0, color='k')

# plt.savefig( os.path.join(folder, 'Figs-Tabs', 'fig_{}.png'.format(fignum)), bbox_inches='tight', dpi = 1000) # format='ps', papertype='a4', orientation='landscape'


#### Parte 06



R1 = 15.03e3	#kOhm
R2 = 98.4e3		#kOhm
RL = 2.26e3		#kOhm

dR1, dR2, dRL = mme([R1, R2, RL], 'ohm')

Vcc = 5			#V
dVcc = mme(5, 'volt')

##ingresso_alto
Vin  = 5		#volt
dVin = mme(Vin, 'volt', 'oscil') + .02*Vin
Vbe  = 644e-3		#millivolt
dVbe = mme(Vbe, 'volt', 'oscil') + .02*Vbe
Vout = 25e-3	#millivolt
dVout = 1.4e-3	#mV

Ib = (Vin - Vbe)/R1
dIb = Ib*np.sqrt( (dVin**2+dVbe**2)/(Vin-Vbe)**2 + dR1**2/R1**2)
Ic = (Vcc - Vout)/RL
dIc = Ic*np.sqrt( (dVcc**2 + dVout**2)/(Vcc-Vout)**2 + dRL**2/RL**2)

##ingresso_basso
Vin= 0.8e-3		#millivolt
dVin = 0.8e-3	#millivolt
Vbe = 0			#microvolt
dVbe = 720e-6	#uV
Vout= 4.96		#volt
dVout = mme(Vout, 'volt', 'oscil') + .02*Vout

Ib = (Vin - Vbe)/R1
dIb = Ib*np.sqrt( (dVin**2+dVbe**2)/(Vin-Vbe)**2 + dR1**2/R1**2)
Ic = (Vcc - Vout)/RL
dIc = Ic*np.sqrt( (dVcc**2 + dVout**2)/(Vcc-Vout)**2 + dRL**2/RL**2)

## 100 kOhm

# input rising edge:
rise_trd = 210e-9	#ns
rise_td  = 250e-9	#ns

#input falling edge:
fall_trs = 4.76e-6	#us
fall_ts  = 1.68e-6	#us


## 3.3 kOhm (approx.)

# input rising edge:
rise_trd = 250e-9	#ns
rise_td  = 750e-9	#ns

#input falling edge:
fall_trs = 900e-9	#ns
fall_ts  = 1.0e-6	#us

## 2.2 kOhm (approx.)

# input rising edge:
rise_trd = 300e-9	#ns
rise_td  = 3.0e-6	#us

#input falling edge:
fall_trs = 200e-9	#ns
fall_ts  = 800e-9	#ns


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
xtype = 'linear'
ytype = 'linear'
rescaleX = 1e0
rescaleY = 1e0
whatplot = () # 'data', 'fit', 'resid', 'line', 'flatline'

resd = (YY - f(XX, *pars)) /dYY
DoFs = len(XX)-len(pars)
print( 'ChiSquare = {0} ({1} DoF, p = {2})'.format( np.sum(resd**2), DoFs, dists.chi2.sf(np.sum(resd**2), DoFs) ) )

if xtype == 'linear':
    points = np.linspace( np.amin(XX) - (np.amax(XX)-np.amin(XX))/20 ,
                          np.amax(XX) + (np.amax(XX)-np.amin(XX))/20,
                          num=max(len(XX)*10, 200) )
elif xtype == 'log':
    points = np.logspace( np.log10( np.amin(XX)/(np.amax(XX)/np.amin(XX))**(1/20) ),
                          np.log10( np.amax(XX)*(np.amax(XX)/np.amin(XX))**(1/20) ),
                          num=max(len(XX)*10, 200) )

if any(whatplot):
    plt.figure(fignum)
    plt.clf()
    if 'resid' in whatplot:
        plt.subplot(subplotGS[:4])
    plt.title(titolo)
    plt.ylabel(labely)
    if 'resid' not in whatplot:
        plt.xlabel(labelx)
    plt.yscale(ytype)
    plt.xscale(xtype)
    if 'data' in whatplot:
        plt.errorbar(XX*rescaleX, YY*rescaleY, dYY*rescaleY, dXX*rescaleX, fmt='none', ecolor='black', label='data')
    if 'fit' in whatplot:
        plt.plot(points*rescaleX, f(points, *pars)*rescaleY, color='red', label='fit')
    if 'line' in whatplot:
        plt.plot(points*rescaleX, line(points, *linepars)*rescaleY, color=greencol)
    if 'flatline' in whatplot:
        plt.plot(points*rescaleX, flatline(points, flatpars)*rescaleY, color=greencol)
    plt.legend()
    if 'resid' in whatplot:
        plt.subplot(subplotGS[4:])
        plt.ylabel('Norm. res.')
        plt.xlabel(labelx)
        plt.xscale(xtype)
        plt.plot(XX*rescaleX, resd , 'ob')
        plt.axhline(y=0, color='k')

# plt.savefig( os.path.join(folder, 'Figs-Tabs', 'fig_{}.png'.format(fignum)), bbox_inches='tight', dpi = 1000) # format='ps', papertype='a4', orientation='landscape'

plt.show()