import sys, os
sys.path.append( os.path.join(os.path.realpath('..'), '00 - Risorse', 'Python') )
folder = os.path.realpath('.')
import numpy as np, pylab, matplotlib.pyplot as plt, matplotlib as mpl, scipy.stats.distributions as dists
from lab import *
from ian import *
from PyTeX import *

plt.close('all')
rawdata = None

#### Parte 000

datafile = 'dati_{0}{1}.txt'.format(1,'')

try:
	rawdata = np.loadtxt( os.path.join(folder, 'Dati', datafile) ).T
except FileNotFoundError:
	print('Hai sbagliato file')

#freq = rawdata[0]
Vout = rawdata[2] - rawdata[3]
Vin = rawdata[0] - rawdata[1]
#Av = 20*np.log10(Vout/Vin)

 fitfun = f
 dfitfun = df

 fitpars, fitpcov = fit_generic_xyerr(fitfun, dfitfun, XX, YY, dXX, dYY, p0=[1])
 linepars, linecov = fit_generic_xyerr(line, dline, XX, YY, dXX, dYY, p0=[1, 0])
 
 resd = (YY - fitfun(XX, *fitpars)) /dYY
 DoFs = len(XX)-len(fitpars)
 tellChi2(resd, DoFs)
 
 usualLines = {
 	fitfun: dict(
 		pars = fitpars,
 		linetype = dict(color='red', label='fit')
 	),
 	line: dict(
 		pars = linepars,
 		linetype = dict(color=greencol)
 	)
 }


first = Graph(Vin, Vout, funFacts = {})
#first.typeX = 'log'
first.compute()
first.draw('data')
	

###### ending
# plt.draw_all()
plt.show()
Graph.countfigs.send(0)
