import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import scipy.stats.distributions as dists
from copy import deepcopy as clone


greencol = (0,1,0)

def compose2(f, g, unpack=False):
	if unpack:
		return (lambda *args: f(*g(*args)))
	else:
		return (lambda *args: f(g(*args)))
		


line = lambda x, m, q: m*x + q
logline = lambda x, m, q: m*np.log10(x) + q
dline = np.vectorize(lambda x, m, q: m)
dlogline = lambda x, m, q: m*np.log10(np.e)/x
const = np.vectorize(lambda x, q: q)

def startfignums(max=100, start=1):
	""" Trivial counting generator, used to produce numbering of Graph() figures
	"""
	while start < max:
		val = (yield start)
		start +=1
		if val is not None:
			start = val

def tellChi2(resd, DoFs):
	chi2msg = 'ChiSquare = {0} ({1} DoF, p = {2})\n'.format(np.sum(resd**2), DoFs, dists.chi2.sf(np.sum(resd**2), DoFs))
	print(chi2msg)
	return chi2msg

usualFun = {line, logline, const}



			
class Graph(object):
	countfigs = startfignums()
	def __init__(self, XX, YY, dXX=0, dYY=0, funFacts=dict(), additional_plots=dict(), what=set()):
		self.figNum = self.countfigs.__next__()
		self.dataX = XX
		self.dataY = YY
		self.errX = dXX
		self.errY = dYY
		self.titolo = 'SONO TITO'
		self.labelX = 'SONO LA X'
		self.labelY = 'SONO LA Y'
		self.typeX = 'linear'
		self.typeY = 'linear'
		self.reX = 1
		self.reY = 1
		self.edgePadding = 1/20
		self.additional = additional_plots.copy()
		self.funPars = funFacts.copy()
		self.whatplot = set(what)
		self.figObj = plt.figure(self.figNum)

	def findextr(self):
		self.minX = np.amin(self.dataX)
		self.maxX = np.amax(self.dataX)
		self.minY = np.amin(self.dataY)
		self.maxY = np.amax(self.dataY)
		return self.minX, self.maxX, self.minY, self.maxY

	def findpts(self, xtype=None, ytype=None):
		if xtype is None:
			xtype = self.typeX
		if ytype is None:
			ytype = self.typeY
		if xtype == 'linear':
			xwidth = self.maxX - self.minX
			xlows = self.minX - xwidth*self.edgePadding
			xhighs = self.maxX + xwidth*self.edgePadding
			self.pts = np.linspace(xlows, xhighs, num=max(len(self.dataX)*10, 200))
		elif xtype == 'log':
			xwidth = self.maxX/self.minX
			xlows = np.log10(self.minX / xwidth**self.edgePadding)
			xhighs = np.log10(self.maxX * xwidth**self.edgePadding)
			self.pts = np.logspace(xlows, xhighs, num=max(len(self.dataX)*10, 200))
		if ytype == 'linear':
			ylows = self.minY - (self.maxY - self.minY)*self.edgePadding
			yhighs = self.maxY + (self.maxY - self.minY)*self.edgePadding
		elif ytype == 'log':
			ylows = np.log10(self.minY / (self.maxY/self.minY)**self.edgePadding)
			yhighs = np.log10(self.maxY * (self.maxY/self.minY)**self.edgePadding)

		return xlows, xhighs, ylows, yhighs, self.pts

	def compute(self):
		self.findextr()
		self.findpts()

	def draw(self, *args, resid=False):
		mainsubs = 1
		self.figObj.clf()
		if args: 
			toplot = set(args)
		else:
			toplot = self.whatplot
		if not resid:
			mainAx = self.figObj.add_subplot(1,1,1)
			mainAx.set_xlabel(self.labelX)
			resid = ()
		else:
			if resid is True: resid = [f for f in toplot if callable(f)]
			if len(resid) < 3:
				mainsubs = 4
			else:
				mainsubs = len(resid)+2
			subGS = mpl.gridspec.GridSpec(5, 1) # TODO: make better
			mainAx = self.figObj.add_subplot(subGS[:4])
			resdAx = self.figObj.add_subplot(subGS[4:])
			resdAx.set_xlabel(self.labelX)
			resdAx.set_xscale(self.typeX)

		mainAx.set_ylabel(self.labelY)
		mainAx.set_title(self.titolo)
		mainAx.set_xscale(self.typeX)
		mainAx.set_yscale(self.typeY)
		
		residuals = dict()

		if 'data' in toplot:
			toplot -= {'data'}
			XX = self.dataX
			YY = self.dataY
			dXX = self.errX
			dYY = self.errY
			mainAx.errorbar(XX*self.reX, YY*self.reY, dYY*self.reY, dXX*self.reX, fmt='none', ecolor='black', label='data')

		for fun in toplot:
			if callable(fun):
				points = self.pts
				mask = self.funPars[fun].get('mask', np.ones(len(XX), bool))
				pars = self.funPars[fun]['pars']
				lineKw = self.funPars[fun].get('linetype', dict())
				mainAx.plot(points*self.reX, fun(points, *pars)*self.reY, **lineKw)
				if 'color' not in  lineKw:
					lineKw['color'] = mainAx[-1].get_color()
				if resid:
					residuals[fun] = (YY[mask] - fun(XX[mask], *pars)) / dYY[mask]
			else:
				self.additional[fun]()

		mainAx.legend()
		
		for fun in resid:
			mask = self.funPars[fun].get('mask', np.ones(len(XX), bool))
			resdAx.plot(XX[mask] * self.reX, residuals[fun], 'o', color=self.funPars[fun]['linetype']['color'])
		
		# plt.savefig( os.path.join(folder, 'Figs-Tabs', 'fig_{}.png'.format(fignum)), bbox_inches='tight', dpi = 1000)
		# format='ps', papertype='a4', orientation='landscape'


	
