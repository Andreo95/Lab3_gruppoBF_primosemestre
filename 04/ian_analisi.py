import sys, os

sys.path.append( os.path.join(os.path.realpath('..'), '00 - Risorse', 'Python') )
folder = os.path.realpath('.')

import numpy as np, pylab, matplotlib.pyplot as plt, matplotlib as mpl, scipy.stats.distributions as dists
from lab import *
from ian import *
from PyTeX import *

plt.close('all')
rawdata = None


#### Parte 2

datafile = 'dati_{0}{1}.txt'.format(2,'')

try:
	rawdata = np.loadtxt( os.path.join(folder, 'Dati', datafile) ).T
except FileNotFoundError:
	print('Hai sbagliato file')

Vin = (rawdata[0] - rawdata[1])
dVin = rawdata[2]
Vout = rawdata[3] - rawdata[4]
dVout = rawdata[5]

used = Vin<1.5
XX  = Vin
YY  = Vout
dXX = dVin/2*.7
dYY = dVout/2*.7

linepars, linecov = fit_generic_xyerr(line, dline, XX[used], YY[used], dXX[used], dYY[used], p0=[10, 0])

print(linepars, np.sqrt(np.diag(linecov)))

resd = (YY[used] - line(XX[used], *linepars)) /dYY[used]
DoFs = len(XX[used])-len(linepars)
tellChi2(resd, DoFs)

plotLines = {
	line: dict(
		pars = linepars,
		linetype = dict(color=greencol),
		mask = used
	)
}


first = Graph(XX, YY, dXX, dYY, funFacts = plotLines)
first.compute()
first.whatplot |= {'data', line}
first.draw(resid=True)
		
# plt.savefig( os.path.join(folder, 'Figs-Tabs', 'fig_{}.png'.format(fignum)), bbox_inches='tight', dpi = 1000) # format='ps', papertype='a4', orientation='landscape'


#### Parte 3

datafile = 'dati_{0}{1}.txt'.format(3,'')

try:
	rawdata = np.loadtxt( os.path.join(folder, 'Dati', datafile) ).T
except FileNotFoundError:
	print('Hai sbagliato file')

freq = rawdata[0]
dfreq = freq/1000
Vout = rawdata[1] - rawdata[2]
dVout = rawdata[3]

Vin = .352 + .440	#Volt
dVin = mme(Vin, 'volt', 'oscil')

Av = 20*np.log10(Vout/Vin)
dAv = 20*np.log10(np.e)*np.sqrt(dVin**2/Vin**2, dVout**2/Vout**2)

used = (freq>500) & (freq<11000)
XX  = freq
YY  = Av
dXX = dfreq
dYY = dAv

flatpar, flatvar = fit_const_yerr(YY[used], dYY[used])

resd = (YY[used] - const(XX[used], flatpar)) /dYY[used]
DoFs = len(XX[used])-1
tellChi2(resd, DoFs)


plotLines.update({
	logline: dict(
		pars = linepars,
		linetype = dict(color=greencol)
	),
	const: dict(
		pars = (flatpar,),
		linetype = dict(color='green'),
		mask = used
	)
})

third = Graph(XX, YY, dXX, dYY, funFacts = plotLines)
third.typeX = 'log'
third.compute()
third.whatplot |= {'data', const}
third.draw(resid=True)
		
# plt.savefig( os.path.join(folder, 'Figs-Tabs', 'fig_{}.png'.format(fignum)), bbox_inches='tight', dpi = 1000) # format='ps', papertype='a4', orientation='landscape'


#### Parte 4

datafile = 'dati_{0}{1}.txt'.format(4,'')

try:
	rawdata = np.loadtxt( os.path.join(folder, 'Dati', datafile) ).T
except FileNotFoundError:
	print('Hai sbagliato file')

Vin = (rawdata[0] - rawdata[1])
dVin = rawdata[2]
Vout = rawdata[3] - rawdata[4]
dVout = rawdata[5]

used = Vin<1000
XX  = Vin
YY  = Vout
dXX = dVin/2*.7
dYY = dVout/2*.7

linepars, linecov = fit_generic_xyerr(line, dline, XX[used], YY[used], dXX[used], dYY[used], p0=[10, 0])

resd = (YY[used] - line(XX[used], *linepars)) /dYY[used]
DoFs = len(XX[used])-len(linepars)
tellChi2(resd, DoFs)

plotLines = {
	line: dict(
		pars = linepars,
		linetype = dict(color=greencol),
		mask = used
	)
}

second = Graph(XX, YY, dXX, dYY, funFacts = plotLines)
second.compute()
second.whatplot |= {'data', line}
second.draw(resid=True)

###### ending
plt.draw_all()
plt.show()
Graph.countfigs.send(0)