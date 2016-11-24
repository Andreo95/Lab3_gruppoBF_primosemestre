import sys, os
if __name__ == '__main__':
	sys.path.append( os.path.join(os.path.realpath('..'), '00 - Risorse', 'Python') )
	folder = os.path.realpath('.')

import numpy as np, pylab, matplotlib.pyplot as plt, matplotlib as mpl, scipy.stats.distributions as dists
from lab import *
from PyTeX import *
#from ian import *

subplotGS = mpl.gridspec.GridSpec(5,1)
plt.close('all')
greencol = (0,1,0)
rawdata = None


#### Parte 000

datafile = 'dati_{0}{1}.txt'.format(4,'')

try:
	rawdata = np.loadtxt( os.path.join(folder, 'Dati', datafile) ).T
except FileNotFoundError:
	print('Hai sbagliato file')

Vin = (rawdata[0] - rawdata[1])
dVin = rawdata[2]
Vout = rawdata[3] - rawdata[4]
dVout = rawdata[5]

Av = Vout/Vin


dmax=0.1
XX  = Vin[Vin<dmax]
YY  = Vout[Vin<dmax]
dXX = dVin[Vin<dmax]/2
dYY = dVout[Vin<dmax]/2



pars, pcov = fit_generic_xyerr(line, dline ,XX, YY, dXX, dYY)
print(pars, pcov)
resd = (YY - line(XX, *pars)) /dYY
DoFs = len(XX)-len(pars)
print('ChiSquare = {0} ({1} DoF, p = {2})\n'.format( np.sum(resd**2), DoFs, dists.chi2.sf(np.sum(resd**2), DoFs) ))

f=line
fignum = 000
titolo = 'SONO TITO'
labelx = 'SONO LA X'
labely = 'SONO LA Y'
xtype = 'linear'
ytype = 'linear'
rescaleX = 1e0
rescaleY = 1e0
whatplot = ("data", "fit") # 'data', 'fit', 'resid', 'line', 'flatline'

if xtype == 'linear':
	xlows  = np.amin(XX) - (np.amax(XX)-np.amin(XX))/20
	xhighs = np.amax(XX) + (np.amax(XX)-np.amin(XX))/20
	points = np.linspace(xlows, xhighs, num=max(len(XX)*10, 200) )
elif xtype == 'log':
	xlows  = np.log10( np.amin(XX)/(np.amax(XX)/np.amin(XX))**(1/20) )
	xhighs = np.log10( np.amax(XX)*(np.amax(XX)/np.amin(XX))**(1/20) )
	points = np.logspace(xlows, xhighs, num=max(len(XX)*10, 200) )
	
if ytype == 'linear':
	ylows  = np.amin(YY) - (np.amax(YY)-np.amin(YY))/20
	yhighs = np.amax(YY) + (np.amax(YY)-np.amin(YY))/20
elif ytype == 'log':
	ylows  = np.log10( np.amin(YY)/(np.amax(YY)/np.amin(YY))**(1/20) )
	yhighs = np.log10( np.amax(YY)*(np.amax(YY)/np.amin(YY))**(1/20) )

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
	plt.xlim(xlows, xhighs)
	plt.ylim(ylows, yhighs)
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
		plt.xlim( (xlows, xhighs) )
		plt.ylim( (-4, 4) )
		plt.plot(XX*rescaleX, resd , 'ob')
		plt.axhline(y=0, color='k')
		
# plt.savefig( os.path.join(folder, 'Figs-Tabs', 'fig_{}.png'.format(fignum)), bbox_inches='tight', dpi = 1000) # format='ps', papertype='a4', orientation='landscape'

plt.show()