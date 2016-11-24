import sys, os
# if __name__ == '__main__':
sys.path.append( os.path.join(os.path.realpath('..'), '00 - Risorse', 'Python') )
folder = os.path.realpath('.')
 
import numpy as np, pylab, matplotlib.pyplot as plt, matplotlib as mpl, scipy.stats.distributions as dists
from lab import *
from ian import *
from PyTeX import *

plt.close('all')
rawdata = None


#### Parte 1

datafile = 'dati_{0}{1}{2}.txt'.format(1,'b', 1)

try:
	filename =os.path.join(folder, 'Dati', datafile)
	rawdata = np.loadtxt( filename ).T
except FileNotFoundError:
	print("{0} non si può aprire".format(filename))

DV, VGS = rawdata
R_1 = 551	# ohm
dDV = mme(DV, "volt")

XX  = VGS
YY  = DV/R_1
dXX = mme(XX, "volt") 
dYY = dDV/R_1

def fitfun(x, a, b):
	return b * (x/a - 1)**2
	
def dfitfun(x, a, b):
	return b/a * 2 * (x/a - 1)
	
mask = XX > -4.5

fitpars, fitpcov = fit_generic_xyerr(fitfun, dfitfun, XX[mask], YY[mask], dXX[mask], dYY[mask], p0=[-4.4, 14e-3])

resd = (YY[mask] - fitfun(XX[mask], *fitpars)) /dYY[mask]
DoFs = len(XX[mask])-len(fitpars)
tellChi2(resd, DoFs)

usualLines = {
	fitfun: dict(
		pars = fitpars,
		linetype = dict(color='red', label='fit'),
		mask = mask
	)
}

first = Graph(XX, YY, dXX, dYY, funFacts = usualLines)
first.compute()
first.reY = 1e3
first.labelY = "$I_D [mA]$"
first.labelX = "$V_{GS} [V]$"
first.titolo = ''
first.draw('data', resid=False)
		
# plt.savefig( os.path.join(folder, 'Figs-Tabs', 'fig_{}.png'.format(fignum)), bbox_inches='tight', dpi = 1000) # format='ps', papertype='a4', orientation='landscape'





############## parte 2
datafile = 'dati_{0}{1}{2}.txt'.format(1,'b', 1)

try:
	filename =os.path.join(folder, 'Dati', datafile)
	rawdata = np.loadtxt( filename ).T
except FileNotFoundError:
	print("{0} non si può aprire".format(filename))

DV, VGS = rawdata
R_1 = 551	# ohm
dDV = mme(DV, "volt")

XX  = VGS
YY  = DV/R_1
dXX = mme(XX, "volt") 
dYY = dDV/R_1

def fitfun(x, a, b):
	return (a/(b**2))*(x-b)**2

def dfitfun(x, a, b):
	return 2*(a/(b**2))*(x-b)

mask = XX > -4.5

fitpars, fitpcov = fit_generic_xyerr(fitfun, dfitfun, XX[mask], YY[mask], dXX[mask], dYY[mask], p0=[0.013,-4.4])

resd = (YY[mask] - fitfun(XX[mask], *fitpars)) /dYY[mask]
DoFs = len(XX[mask])-len(fitpars)
tellChi2(resd, DoFs)

usualLines = {
	fitfun: dict(
		pars = fitpars,
		linetype = dict(color='red', label='fit'),
		mask = mask
	)
}

first = Graph(XX, YY, dXX, dYY, funFacts = usualLines)
first.compute()
first.reY = 1e3
first.labelY = "$I_D [mA]$"
first.labelX = "$V_{GS} [V]$"
first.titolo = ''
first.draw('data', resid=False)

first.funPars = usualLines
first.compute()
first.whatplot |= {'data', fitfun}
first.draw(resid=[])



##############parte 3
datafile = 'dati_{0}.txt'.format(3)

try:
	filename =os.path.join(folder, 'Dati', datafile)
	rawdata = np.loadtxt( filename ).T
except FileNotFoundError:
	print("{0} non si può aprire".format(filename))

rawdata /= 1000
Vin = rawdata[1] - rawdata[0]
Vs = rawdata[3] - rawdata[2]
Vd = rawdata[5] - rawdata[4]
R_1 = 551



XX  = Vin
YY  = Vd/Vin
dXX = mme(XX, "volt", "oscil") 
dYY = YY*((mme(Vin, "volt", "oscil")/Vin)**2 + (mme(Vd, "volt", "oscil")/Vd)**2)**0.5


second = Graph(XX, YY, dXX, dYY, funFacts = usualLines)
second.compute()
second.draw('data')



XX  = Vin
YY  = Vd 
dXX = mme(XX, "volt", "oscil")
dYY = mme(YY, "volt", "oscil")

mask = XX < 4.5
mmx=np.mean(dXX[mask])
mmy=np.mean(dYY[mask])
linepars, linecov = fit_generic_xyerr(line, dline, XX[mask], YY[mask], dXX[mask], dYY[mask], p0=[1, 0])

print("Common source")
print(linepars)
print(linecov[0, 0]**0.5)
print(linecov[1, 1]**0.5)
print(linecov[0, 1]/(linecov[0, 0]*linecov[1, 1])**0.5)


resd = (YY[mask] - line(XX[mask], *linepars)) /dYY[mask]
DoFs = len(XX[mask])-len(linepars)
tellChi2(resd, DoFs)

usualLines = {
	line: dict(
		pars = linepars,
		linetype = dict(color=greencol),
		mask = mask
	)
}


third = Graph(XX, YY, dXX, dYY, funFacts = usualLines)
third.titolo="$V_{in}$ vs $V_{out}$ (common source)"
third.labelX="$V_{in}$ [V]"
third.labelY="$V_{out}$(drain) [V]"
third.compute()
third.draw('data', line, resid=True)
		
# plt.savefig( os.path.join(folder, 'Figs-Tabs', 'fig_{}.png'.format(fignum)), bbox_inches='tight', dpi = 1000) # format='ps', papertype='a4', orientation='landscape'


XX  = Vin
YY  = Vs 
dXX = mme(XX, "volt", "oscil")
dYY = mme(YY, "volt", "oscil")

mask = XX < 5
mmx=np.mean(dXX[mask])
mmy=np.mean(dYY[mask])
linepars, linecov = fit_generic_xyerr(line, dline, XX[mask], YY[mask], dXX[mask], dYY[mask], p0=[1, 0])

print("Source follower")
print(linepars)
print(linecov[0, 0]**0.5)
print(linecov[1, 1]**0.5)
print(linecov[0, 1]/(linecov[0, 0]*linecov[1, 1])**0.5)

resd = (YY[mask] - line(XX[mask], *linepars)) /dYY[mask]
DoFs = len(XX[mask])-len(linepars)
tellChi2(resd, DoFs)

usualLines = {
	line: dict(
		pars = linepars,
		linetype = dict(color=greencol),
		mask = mask
	)
}


dVs = mme(Vs, 'volt','oscil')
fourth = Graph(Vin, Vs, dXX, dVs, funFacts=usualLines)
fourth.titolo="$V_{in}$ vs $V_{out}$ (source follower)"
fourth.labelX="$V_{in}$ [V]"
fourth.labelY="$V_{out}$(source) [V]"
fourth.compute()
fourth.draw('data', line, resid=True)


###### ending
plt.show()
Graph.countfigs.send(0)
