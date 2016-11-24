


import sys, os
sys.path.append( os.path.join(os.path.realpath('..'), '00 - Risorse', 'Python') )

import numpy as np, pylab, matplotlib.pyplot as plt, matplotlib as mpl
from lab import *
from PyTeX import *

pylab.close('all')

folder = os.path.realpath('.')

subplotGS = mpl.gridspec.GridSpec(5,1)
 
line = lambda x, m, q: m*x+q
dline = np.vectorize(lambda x, m, q : m)
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




#### Parte 1

datafile = '{0}{1}_dati.txt'.format('',2)

rawdata = np.loadtxt( os.path.join(folder, 'Dati', datafile) ).T

linedata = rawdata.T[rawdata[0]>14].T


Vin = linedata[2]
Vout = linedata[4]
frq = linedata[0]*1000


XX  = np.log10(frq)
YY  = 20*np.log10(Vout/Vin)
dXX = linedata[1]/frq
dYY = np.sqrt(2*linedata[3]**2/Vin**2 + 2*linedata[5]**2/Vout**2)*20*np.log10(np.e)

linepars, linecov = fit_generic_xyerr(line, dline, XX, YY, dXX, dYY, p0=[-20, 60])

linedata = rawdata.T[rawdata[0]<.100].T


Vin = linedata[2]
Vout = linedata[4]
frq = linedata[0]*1000


XX  = np.log10(frq)
YY  = 20*np.log10(Vout/Vin)
dXX = linedata[1]/frq
dYY = np.sqrt(2*linedata[3]**2/Vin**2 + 2*linedata[5]**2/Vout**2)*20*np.log10(np.e)

# flatpars, flatcov = fit_generic_xyerr(line, dline, XX, YY, dXX, dYY, p0=[0, 0])
flatpars, flatcov = fit_const_yerr(YY, dYY)

ff=FindCross(linepars, linecov, flatpars, flatcov)
pt=10**ff[0]
errpt=(10**ff[0])*np.log(10)*ff[1]**0.5
print(10**ff[0], (10**ff[0])*np.log(10)*ff[1]**0.5)

Vin = rawdata[2]
Vout = rawdata[4]
frq = rawdata[0]*1000


XX  = frq
YY  = 20*np.log10(Vout/Vin)
dXX = rawdata[1]
dYY = np.sqrt(2*rawdata[3]**2/Vin**2 + 2*rawdata[5]**2/Vout**2)*20*np.log10(np.e)


f = lambda x, ft, c: c -10*np.log10(1+(x/ft)**2)
df = lambda x, ft, c: -10*log10(np.e)*2*x/(x**2+ft**2)

pars, pcov = fit_generic_xyerr(f, df, XX, YY, dXX, dYY, p0=[2000, 0])
print(pars, pcov)


print(f(1e100, *pars) - f(1e99, *pars))


fignum = 1
titolo = 'Passa basso'
labelx = 'Frequanza [Hz]'
labely = 'Attenuazione [dB]'
rescaleX = 1e0
rescaleY = 1e0
xtype = 'log'


resd = (YY - f(XX, *pars)) /dYY
print( 'ChiSquare = {0} ({1} DoF)'.format(np.sum(resd**2), len(XX)-len(pars)) )

if xtype == 'linear':
    points = np.linspace( np.amin(XX) - (np.amax(XX)-np.amin(XX))/10 , np.amax(XX) + (np.amax(XX)-np.amin(XX))/10, num=max(len(XX)*10, 200) )
elif xtype == 'log':
    points = np.logspace( np.log10(np.amin(XX) / (np.amax(XX)/np.amin(XX))**0.1 ), np.log10(np.amax(XX) * (np.amax(XX)/np.amin(XX))**0.1), num=max(len(XX)*10, 200) )

plt.figure(fignum)
plt.clf()
plt.subplot(subplotGS[:4])
plt.title(titolo)
plt.ylabel(labely)
plt.yscale('linear')
plt.xscale(xtype)
plt.errorbar(XX*rescaleX, YY*rescaleY, dYY*rescaleY, dXX*rescaleX, fmt='none', ecolor='black', label='data')
plt.plot(points*rescaleX, f(points, *pars)*rescaleY, color='red', label='fit')
plt.plot(points[points>500], line(np.log10(points[points>500]),linepars[0],linepars[1]), color=(0, 1, 0.2))
# plt.plot(points, line(np.log10(points),flatpars[0],flatpars[1]), color='b')
plt.plot(points, flatline(np.log10(points),flatpars), color=(0, 1, 0.2))
plt.errorbar([pt], [flatpars], flatcov, errpt, fmt="none", color="b", label="intersezione")
plt.legend()
plt.subplot(subplotGS[4:])
plt.ylabel('Norm. res.')
plt.xlabel(labelx)
plt.xscale(xtype)
plt.plot(XX*rescaleX, resd , 'ob')
plt.axhline(y=0, color='k')

# plt.savefig( os.path.join(folder, 'Figs-Tabs', 'fig_{}.png'.format(fignum)), bbox_inches='tight', dpi = 1000) # format='ps', papertype='a4', orientation='landscape'

#### Parte 1-guadagno non dB

datafile = '{0}{1}_dati.txt'.format('',2)

rawdata = np.loadtxt( os.path.join(folder, 'Dati', datafile) ).T



Vin = rawdata[2]
Vout = rawdata[4]
frq = rawdata[0]*1000


XX  = frq
YY  = Vout/Vin
dXX = rawdata[1]
dYY =  np.sqrt((2*rawdata[3]**2/Vin**2 + 2*rawdata[5]**2/Vout**2))*YY


f = lambda x, ft, c: c/np.sqrt(1+(x/ft)**2)
df = lambda x, ft, c: -c/ft/np.sqrt(1+(x/ft)**2)**3

pars, pcov = fit_generic_xyerr(f, df, XX, YY, dXX, dYY, p0=[2000, 1])





fignum = 2
titolo = 'SONO TITO'
labelx = 'SONO LA X'
labely = 'SONO LA Y'
rescaleX = 1e0
rescaleY = 1e0
xtype = 'log'


resd = (YY - f(XX, *pars)) /dYY
print( 'ChiSquare = {0} ({1} DoF)'.format(np.sum(resd**2), len(XX)-len(pars)) )

if xtype == 'linear':
    points = np.linspace( np.amin(XX) - (np.amax(XX)-np.amin(XX))/10 , np.amax(XX) + (np.amax(XX)-np.amin(XX))/10, num=max(len(XX)*10, 200) )
elif xtype == 'log':
    points = np.logspace( np.log10(np.amin(XX) / (np.amax(XX)/np.amin(XX))**0.1 ), np.log10(np.amax(XX) * (np.amax(XX)/np.amin(XX))**0.1), num=max(len(XX)*10, 200) )

plt.figure(fignum)
plt.clf()
plt.subplot(subplotGS[:4])
plt.title(titolo)
plt.ylabel(labely)
plt.yscale('linear')
plt.xscale(xtype)
plt.errorbar(XX*rescaleX, YY*rescaleY, dYY*rescaleY, dXX*rescaleX, fmt='none', ecolor='black', label='data')
plt.plot(points*rescaleX, f(points, *pars)*rescaleY, color='red', label='fit')
plt.legend()
plt.subplot(subplotGS[4:])
plt.ylabel('Norm. res.')
plt.xlabel(labelx)
plt.xscale(xtype)
plt.plot(XX*rescaleX, resd , 'ob')
plt.axhline(y=0, color='k')

# plt.savefig( os.path.join(folder, 'Figs-Tabs', 'fig_{}.png'.format(fignum)), bbox_inches='tight', dpi = 1000) # format='ps', papertype='a4', orientation='landscape'

#### Parte 2

datafile = '{0}{1}_dati.txt'.format('p2-',4)

rawdata = np.loadtxt( os.path.join(folder, 'Dati', datafile) ).T


R1 = 3.22   *1e3    #kOhm
R2 = 3.28   *1e3    #kOhm
C1 = 10.70  *1e-9   #nF
C2 = 110.5  *1e-9   #nF

tftb = 1/(2*np.pi*R1*C1)
tfta = 1/(2*np.pi*R2*C2)
RR = R1/R2

linedata = rawdata.T[rawdata[0]<200].T


Vin = linedata[1]-linedata[2]
Vout = linedata[4]-linedata[5]
frq = linedata[0]


XX  = np.log10(frq)
YY  = 20*np.log10(Vout/Vin)
dXX = 1e-3
dYY = np.sqrt(2*linedata[3]**2/Vin**2 + 2*linedata[6]**2/Vout**2)*20*np.log10(np.e)

line1, line1cov = fit_generic_xyerr(line, dline, XX, YY, dXX, dYY, p0=[-20, 60])

linedata = rawdata.T[rawdata[0]>27000].T


Vin = linedata[1]-linedata[2]
Vout = linedata[4]-linedata[5]
frq = linedata[0]


XX  = np.log10(frq)
YY  = 20*np.log10(Vout/Vin)
dXX = 1e-3
dYY = np.sqrt(2*linedata[3]**2/Vin**2 + 2*linedata[6]**2/Vout**2)*20*np.log10(np.e)

line2, line2cov = fit_generic_xyerr(line, dline, XX, YY, dXX, dYY, p0=[-20, 60])

linedata = rawdata.T[(rawdata[0]>600) & (rawdata[0]<3000)].T


Vin = linedata[1]-linedata[2]
Vout = linedata[4]-linedata[5]
frq = linedata[0]


XX  = np.log10(frq)
YY  = 20*np.log10(Vout/Vin)
dXX = 1e-3
dYY = np.sqrt(2*linedata[3]**2/Vin**2 + 2*linedata[6]**2/Vout**2)*20*np.log10(np.e)

flatpars, flatcov = fit_const_yerr(YY, dYY)

pt1, err1 = FindCross(line1, line1cov, flatpars, flatcov)
pt2, err2 = FindCross(line2, line2cov, flatpars, flatcov)
pt1 = 10**pt1
pt2 = 10**pt2
err1 = pt1*np.log(10)*np.sqrt(err1)
err2 = pt2*np.log(10)*np.sqrt(err2)
print(pt1, err1)
print(pt2, err2)

Vin = rawdata[1] - rawdata[2]
Vout = rawdata[4] - rawdata[5]
frq = rawdata[0]*1


XX  = frq
YY  = 20*np.log10(Vout/Vin)
dXX = XX/1000
dYY = np.sqrt(2*rawdata[3]**2/Vin**2 + 2*rawdata[6]**2/Vout**2)*20*np.log10(np.e)


pb = lambda x, ft, c: c/(1+1j*(x/ft))
dpb = lambda x, ft, c: -c/ft/np.sqrt(1+(x/ft)**2)**3


pa = lambda x, ft, c: c/(1+1j*(ft/x))
dpa = lambda x, ft, c: c*ft**2/np.sqrt(1+(ft/x)**2)**3/x**3

g=lambda x, fta, ftb, r: 20*np.log10(np.abs(pb(x, ftb, 1)*pa(x, fta, 1)/(1+r*pa(x,fta, 1)*pb(x, ftb, 1))))
dg=lambda x, fta, ftb, r: np.abs(pa(x, fta, 1)*dpb(x, ftb, 1)+(pb(x, ftb, 1)*dpa(x, fta, 1))*(1+r*pb(x, ftb, 1)*pa(x, fta, 1))/(1+r*pb(x, ftb, 1)*pa(x, fta, 1)))/f(x, fta, ftb, r)

f=g; df =dg
# f = lambda x, fta, ftb : g(x, fta, ftb, RR)
# df = lambda x, fta, ftb : dg(x, fta, ftb, RR)

pars, pcov = fit_generic_xyerr(f, df, XX, YY, 0, dYY, p0=[400, 4000, 1])


def rlcSolve(w, R1, C1, R2, C2):
    Vin  = 1
    Zout = R2 + 1/(1j*w*C2)
    Z2   = 1/(1j*w*C1 + 1/Zout)
    Ztot = Z2 + R1
    Itot = Vin/Ztot
    Iout = Itot*Z2/Zout
    Vout = Iout*R2
    return Vout

lastcat = lambda f: 20*np.log10( np.abs(rlcSolve(f*2*np.pi, R1, C1, R2, C2)) )
    


fignum = 3
titolo = 'Filtro passa-banda'
labelx = 'Frequenza [Hz]'
labely = 'Attenuazione [dB]'
rescaleX = 1e0
rescaleY = 1e0
xtype = 'log'


resd = (YY - f(XX, *pars)) /dYY
# resd = (YY - lastcat(XX)) /dYY
print( 'ChiSquare = {0} ({1} DoF)'.format(np.sum(resd**2), len(XX)-len(pars)) )

if xtype == 'linear':
    points = np.linspace( np.amin(XX) - (np.amax(XX)-np.amin(XX))/20 , np.amax(XX) + (np.amax(XX)-np.amin(XX))/20, num=max(len(XX)*10, 200) )
elif xtype == 'log':
    points = np.logspace( np.log10(np.amin(XX) / (np.amax(XX)/np.amin(XX))**0.05 ), np.log10(np.amax(XX) * (np.amax(XX)/np.amin(XX))**0.05), num=max(len(XX)*10, 200) )
    
plt.figure(fignum)
plt.clf()
plt.subplot(subplotGS[:4])
plt.title(titolo)
plt.ylabel(labely)
plt.yscale('linear')
plt.xscale(xtype)
plt.errorbar(XX*rescaleX, YY*rescaleY, dYY*rescaleY, dXX*rescaleX, fmt='none', ecolor='black', label='data')
plt.plot(points*rescaleX, f(points, *pars)*rescaleY, color='red', label='fit')

# RR=3.22/3.28
# tfta=1/(2*np.pi*3.22e3*10.7e-9)
# tftb=1/(2*np.pi*3.28e3*110.5e-9)


plt.plot(points[points>2800]*rescaleX, line(np.log10(points[points>2800]), *line2)*rescaleY, color=(0, 1, 0.2))
plt.plot(points[points<1000]*rescaleX, line(np.log10(points[points<1000]), *line1)*rescaleY,color=(0, 1, 0.2))
plt.plot(points*rescaleX, flatline(np.log10(points), flatpars)*rescaleY, color=(0, 1, 0.2) )
plt.errorbar([pt1, pt2], [flatpars,flatpars],np.sqrt([flatcov, flatcov]),[err1, err2], fmt='none', color='b', label="intersezioni")


plt.legend()
plt.subplot(subplotGS[4:])
plt.ylabel('Norm. res.')
plt.xlabel(labelx)
plt.xscale(xtype)
plt.plot(XX*rescaleX, resd , 'ob')
plt.axhline(y=0, color='k')

#plt.savefig( os.path.join(folder, 'Figs-Tabs', 'fig_{}.png'.format(fignum)), bbox_inches='tight', dpi = 1000) # format='ps', papertype='a4', orientation='landscape'



fignum = 4
titolo = 'Filtro passa-banda, curve teoriche'
labelx = 'Frequenza [Hz]'
labely = 'Attenuazione [dB]'
xtype = 'log'
ytype = 'linear'
rescaleX = 1e0
rescaleY = 1e0


resd = (YY - lastcat(XX)) /dYY
print( 'ChiSquare = {0} ({1} DoF)'.format(np.sum(resd**2), len(XX)-len(pars)) )

if xtype == 'linear':
    points = np.linspace( np.amin(XX) - (np.amax(XX)-np.amin(XX))/20 , np.amax(XX) + (np.amax(XX)-np.amin(XX))/20, num=max(len(XX)*10, 200) )
elif xtype == 'log':
    points = np.logspace( np.log10( np.amin(XX)/(np.amax(XX)/np.amin(XX))**0.05 ), np.log10( np.amax(XX)*(np.amax(XX)/np.amin(XX))**0.05), num=max(len(XX)*10, 200) )

plt.figure(fignum)
plt.clf()
plt.subplot(subplotGS[:4])
plt.title(titolo)
plt.ylabel(labely)
plt.yscale(ytype)
plt.xscale(xtype)
plt.errorbar(XX*rescaleX, YY*rescaleY, dYY*rescaleY, dXX*rescaleX, fmt='none', ecolor='black', label='data')

plt.plot(points*rescaleX, lastcat(points)*rescaleY, color='blue', label='predicted')
plt.plot(points*rescaleX, g(points, tfta, tftb, RR)*rescaleY, color='green', label='blackbox')

plt.plot(points[points>8600]*rescaleX, line(np.log10(points[points>8600]), -20, 20*np.log10(tftb*2)-6 )*rescaleY, color=(0, 1, 0.2))
plt.plot(points[points<220]*rescaleX, line(np.log10(points[points<220]), 20, -20*np.log10(tfta/2)-6 )*rescaleY,color=(0, 1, 0.2))

plt.legend()
plt.subplot(subplotGS[4:])
plt.ylabel('Norm. res.')
plt.xlabel(labelx)
plt.xscale(xtype)
plt.plot(XX*rescaleX, resd , 'ob')
plt.axhline(y=0, color='k')

# plt.savefig( os.path.join(folder, 'Figs-Tabs', 'fig_{}.png'.format(fignum)), bbox_inches='tight', dpi = 1000) # format='ps', papertype='a4', orientation='landscape'

plt.show()

plt.show()
