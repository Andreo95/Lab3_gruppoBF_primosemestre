"""
pyan module for doing stuff.

Module contains:
	- objects used for fitting and graphing data
	- useful functions
	- random crap
	- the color green
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import inspect
import scipy.stats.distributions as dists
from copy import deepcopy as clone

__all__ = [		# things imported when you do "import *"
	'greencol',
	'dline',
	'dlogline',
	'nullfunc',
	'const',
	'createline',
	'startfignums',
	'tell_chi2',
	'maketab',
	'Fitter',
	'Graph',
]

greencol = (0,1,0)

dline = np.vectorize(lambda x, m, q: m)
dlogline = np.vectorize(lambda x, m, q: m*np.log10(np.e)/x)
nullfunc = np.vectorize(lambda *args: 0)
const = np.vectorize(lambda x, q: q)


def createline(type='linear', name=""):
	"""
	Factory of linear(-like) functions.

	Parameters
	----------
	type: str, optional (default 'linear')
		if 'log', returns a function linear in log10(x).
	name: str, optional
		if provided, the function returned will have this name,
		otherwise it will be named "line" (or "logline" if type is 'log').

	Returns
	-------
	func: callable
		func(x, m, q) = mx + q (or m*log10(x) + q if type is 'log').
	"""
	if type == 'log':
		def logline(x, m, q):
			"""f: (x, m, q) --> m*log10(x) + q ."""
			return m*np.log10(x) + q
		logline.deriv = dlogline
		func = logline
	elif type == 'linear':
		def line(x, m, q):
			"""f: (x, m, q) --> m*x + q ."""
			return m*x + q
		line.deriv = dline
		func = line
	elif type == 'const':
		def flatline(x, q):
			"""f: (x, q) --> q ."""
			return q*np.ones(len(x))
		flatline.deriv = nullfunc
		func = flatline
	if name:
		func.__name__ = name
	return func


def startfignums(max=100, start=1):
	"""Trivial counting generator, used to produce numbering of Graph() figures."""
	while start < max:
		val = (yield start)
		start += 1
		if val is not None:
			start = val


def compose2(f, g, unpack=False):
	"""Mostly useless."""
	if unpack:
		return (lambda *args: f(*g(*args)))
	else:
		return (lambda *args: f(g(*args)))


def tell_chi2(resd, dofs):
	"""
	Chi^2 prettyprinter.

	Calculate chi squared from normalized residuals
	and prints it along with p-value and degrees of freedom.

	Parameters
	----------
	resd: (iterable of) number
		the normalised residuals of a fitted dataset
	dofs: unsigned integer
		the number of degrees of freedom of the fit

	Returns
	-------
	chi2msg: string
		the message printed
	"""
	chi2 = np.sum(resd**2)
	pval = dists.chi2.sf(chi2, dofs)
	chi2msg = "ChiSquare = {0:.2f} ({1} DoF, p = {2:.4f})".format(chi2, dofs, pval)
	print(chi2msg)
	return chi2msg


def maketab(*columns, errors='all', precision=3):
	"""
	Print data in tabular form.

	Creates a LaTeX (table environment with S columns from siunitx) formatted
	table using inputs as columns; the input is assumed to be numerical.
	Errors can be given as columns following the relevant numbers column
	and their positions specified in a tuple given as the kwarg 'errors';
	they will be nicely formatted for siunitx use.
	The precision (number of signifitant digits) of the numbers representation
	will be inferred from errors, or when absent from the 'precision' kwarg.

	Parameters
	----------
	*columns: N iterables of lenght L with numerical items
		the lists of numbers that will make up the columns of the table,
		must be of uniform lenght.
	errors: tuple(-like) of ints, 'all' or 'none', optional (default 'all')
		the columns with positions corresponging to errors items will
		be considered errors corresponding to the previous column and formatted
		accordingly; if 'all' every other column will be considered an error,
		if 'none' no column will be considered an error.
	precision: (lenght-L tuple of) int, optional (default 3)
		number of significant digits to be used when error is not given;
		if a tuple(-like) is given, each column will use the correspondingly
		indexed precision.

	Returns
	-------
	tab: string
		the formatted text constituting the LaTeX table.
	"""
	vals = np.asarray(columns).T
	cols = np.alen(vals.T)
	precision = np.array(precision) * np.ones(cols)
	if errors == 'all':
		errors = range(1, cols, 2)
	if errors == 'none':
		errors = []
	beginning = (
		"\\begin{table}\n" +
		"\t\\begin{tabular}{*{" +
		str(cols if errors is None else cols-len(errors)) +
		"}{S}} \n" +
		"\t\t\n" +
		"\t\t\\hline \n"
	)
	inner = ""
	for i, row in enumerate(vals):
		rows = enumerate(row, start=1)
		tab += "\t"
		for pos, v in rows:
			num = v if np.isfinite(v) else "{-}"
			space = "&" if pos < cols else "\\\\ \n"
			err = ""
			prec = -precision[pos-1]
			if pos in errors:
				prec = np.floor(np.log10(row[pos]))
				if row[pos]/10**prec < 2.5:
					prec -= 1
				err = "({0:.0f})".format(round(row[pos]/10**prec))
				if next(rows, (-1, None))[0] >= cols:
					space = "\\\\ \n"
				num = round(num, int(-prec))
			inner += "\t{0:.{digits:.0f}f} {1}\t{2}".format(num, err, space, digits=max(0, -prec))
	ending = "\t\\end{tabular} \n" + "\t\\caption{some caption} \n" + "\t\\label{t:somelabel} \n" + "\\end{table}"
	return beginning + inner + ending


class Fitter(object):
	"""
	Object containing a dataset to be fitted and otherwise manipulated.

	Methods defined here:

	fit(self, *functions, verbose=True, **kwargs)
		Fit functions to dataset.

	Notes
	-----
	Class attribute _fitter_func needs to be initialized to the fitting function
	to be used.
	"""

	_fitter_func = None

	def __init__(self, *dataset, name="fit"):
		"""
		Initialize fitter object from a dataset (x, y, dx, dy).

		Parameters
		----------
		*dataset: four numpy array-like object
			the four arrays are interpreted as x, y, sigmax, sigmay; x and y
			must be of the same lenght L, sigmax and sigmay can each be either
			of lenght L or of lenght 1 (in which case they don't in fact need to
			be an array).
		name: string, optional (default "fit")
			name of the Fitter object, used when printing results and the like.
		"""
		self.data = np.array(dataset)
		self.fitted = dict()
		self.name = name

	def fit(self, *functions, verbose=True, **kwargs):
		"""
		Fit the functions to the dataset.

		Parameters
		----------
		*functions: any number of callables f
			if f has a (callable) .deriv attribute, it will be used as its
			derivative to propagate x errors; if f has a .mask attribute
			(numpy array-like of booleans with the same shape as the x data) it
			will be used to determine which datapoints should be used for the
			fit.
		verbose: boolean, optional (default True)
			if True, fit results are printed

		Additional keyword arguments are passed directly to the underlying
		fitting function.

		Returns
		-------
		This method doesn't return anything but stores fit results for each
		fitted callable in its attributes (.pars, .cov, .sigmas, .resd).
		"""
		xxs = self.data[0]
		if verbose:
			fitmsg = "Il {fitname} di {funname} ha dato i seguenti parametri:"
			print("Lavoro su {}\n".format(self.name))
		for f in functions:
			mask = getattr(f, 'mask', np.ones(len(xxs), dtype=bool))
			data = self.data[:, mask]
			df = getattr(f, 'deriv', nullfunc)
			p0 = getattr(f, 'pars', None)
			popt, pcov = self._fitter_func(f, df, *data, p0=p0, **kwargs)
			f.pars = popt
			f.cov = pcov
			f.sigmas = np.sqrt(np.diag(pcov))
			f.resd = (data[1] - f(data[0], *popt)) / np.sqrt(data[3]**2 + data[2]**2*df(data[0], *popt)**2)
			if verbose:
				print(fitmsg.format(fitname=self.name, funname=f.__name__))
				argnames = inspect.getargspec(f).args[1:]
				for name, par, err in zip(argnames, popt, f.sigmas):
					print("{0} = {1:.4f} \pm {2:.4f}".format(name, par, err))
				tell_chi2(f.resd, np.alen(data[0]) - len(popt))
				print("")
		if verbose:
			print("{} completo\n\n".format(self.name))


class Graph(object):
	"""
	Object containing data to be plotted and information about the plot style.

	It's an interface to matplotlib functions.

	Class methods defined here:

	from_fitter(cls, fit_obj, *args, **kw)
		Return Graph instance initialized with data from a Fitter object.

	Methods defined here:

	set_edges(self, xtype=None, ytype=None)
		Set the points to be used when plotting functions from data extremes.
		Return xlows, xhighs, ylows, yhighs, self.pts

	draw(self, *toplot, resid=False, docalcs=True)
		Draw toplot items (callables) on the figure of the object.

	Useful instance attributes:

	.title (string)
		Title written on the plot.
	.labelX/Y (string)
		Axes labels written on the plot.
	.typeX/Y ('linear' or 'log', defaults to 'linear')
		If 'log', plot will use logarithmic scale for the relevant axis.
	"""

	countfigs = startfignums()

	def __init__(self, x, y, dx=0, dy=0, additional_plots=dict()):
		"""Same old initializer."""
		self.fig_num = self.countfigs.__next__()
		self.dataX = x
		self.dataY = y
		self.errX = dx
		self.errY = dy
		self.title = "SONO TITO"
		self.labelX = "SONO LA X"
		self.labelY = "SONO LA Y"
		self.typeX = 'linear'
		self.typeY = 'linear'
		self.reX = 1
		self.reY = 1
		self.edge_padding = 1/20
		self.datakw = dict(fmt='none', ecolor='black', label='data')
		self.additional = additional_plots.copy()
		self.fig_obj = plt.figure(self.fig_num)

	@classmethod
	def from_fitter(cls, fit_obj, *args, **kw):
		"""Return Graph instance initialized with data from a Fitter object."""
		return cls(*fit_obj.data, *args, **kw)

	def set_edges(self, xtype=None, ytype=None):
		# TODO: document
		"""
		Implemented but undocumented, sry.

		Parameters
		----------


		Returns
		-------

		"""
		minX = np.amin(self.dataX)
		maxX = np.amax(self.dataX)
		minY = np.amin(self.dataY)
		maxY = np.amax(self.dataY)
		if xtype is None:
			xtype = self.typeX
		if ytype is None:
			ytype = self.typeY
		if xtype == 'linear':
			xwidth = maxX - minX
			xlows = minX - xwidth*self.edge_padding
			xhighs = maxX + xwidth*self.edge_padding
			self.pts = np.linspace(xlows, xhighs, num=max(len(self.dataX)*10, 200))
		elif xtype == 'log':
			xwidth = maxX/minX
			xlows = minX / xwidth**self.edge_padding
			xhighs = maxX * xwidth**self.edge_padding
			self.pts = np.logspace(np.log10(xlows), np.log10(xhighs), num=max(len(self.dataX)*10, 200))
		if ytype == 'linear':
			ylows = minY - (maxY - minY)*self.edge_padding*3
			yhighs = maxY + (maxY - minY)*self.edge_padding*3
		elif ytype == 'log':
			ylows = np.log10(minY / (maxY/minY)**self.edge_padding**3)
			yhighs = np.log10(maxY * (maxY/minY)**self.edge_padding**3)

		return xlows, xhighs, ylows, yhighs, self.pts

	def draw(self, *toplot, resid=False, docalcs=True):
		"""
		Plot items on the figure of the object.

		The data is always plotted as an errorbar plot.

		Parameters
		----------
		*toplot: any number of callables f
			functions to be plotted in addition to the data, must have a .pars
			attribute (tuple) containing function parameters; if f has a
			.bounds attribute (list-like of pairs-like) it will be used to
			determine where to draw the functions (in the intervals with
			extremes given by the pairs); if f has a .linekw attribute
			(dictionary), it will be passed as keyword arguments for the plot of
			the function.
		resid: boolean or list-like of callables f, optional (default False)
			if False, doesn't plot residuals, else plots residuals of the
			callables (masking data with f.mask if available, using f.resid
			if present, else calculating them); if True all the toplot callables
			are used.
		docalcs: boolean, optional (default True)
			if True, self.set_edges is called before drawing.


		Returns
		-------
		This method doesn't return anything but assigns f.linekw['color'] to all
		callables where it isn't already set and assigns plot and residuals axes
		to self attributes.
		"""
		if docalcs:
			minp, maxp, minh, maxh, _ = self.set_edges()
		if not resid:
			main_ax = self.fig_obj.add_subplot(1, 1, 1)
			main_ax.set_xlabel(self.labelX)
			resid = ()
		else:
			if resid is True:
				resid = toplot
			# if len(resid) < 3:
			# 	mainsubs = 4
			# else:
			# 	mainsubs = len(resid)+2
			subGS = mpl.gridspec.GridSpec(5, 1)		# TODO: make better
			main_ax = self.fig_obj.add_subplot(subGS[:4])
			resd_ax = self.fig_obj.add_subplot(subGS[4:])
			resd_ax.axhline(y=0, color='black')
			resd_ax.set_xlabel(self.labelX)
			resd_ax.set_xscale(self.typeX)
			resd_ax.set_xlim(minp*self.reX, maxp*self.reX)
			self.resd_ax = resd_ax

		self.main_ax = main_ax
		main_ax.set_xlim(minp*self.reX, maxp*self.reX)
		main_ax.set_ylim(minh*self.reY, maxh*self.reY)
		main_ax.set_ylabel(self.labelY)
		main_ax.set_title(self.title)
		main_ax.set_xscale(self.typeX)
		main_ax.set_yscale(self.typeY)

		x = self.dataX*self.reX
		y = self.dataY*self.reY
		dx = self.errX*self.reX
		dy = self.errY*self.reY
		main_ax.errorbar(x, y, dy, dx, **self.datakw)

		for fun in toplot:
			if callable(fun):
				mask = np.zeros(len(self.pts), dtype=bool)
				for lowest, highest in getattr(fun, 'bounds', [(-np.inf, np.inf)]):
					mask |= (self.pts > lowest) & (self.pts < highest)
				points = self.pts[mask]
				try:
					linekw = fun.linekw
				except AttributeError:
					linekw = fun.linekw = {}
				g, = main_ax.plot(points*self.reX, fun(points, *fun.pars)*self.reY, **linekw)
				if 'color' not in linekw:
					linekw['color'] = g.get_color()
			else:
				self.additional[fun]()

		main_ax.legend()
		x = self.dataX
		for fun in resid:
			mask = getattr(fun, 'mask', np.ones(len(x), dtype=bool))
			resdkw = dict(marker='o')
			if hasattr(fun, 'linekw'):
				resdkw.update(fun.linekw)
			resdkw.update(ls='none')
			try:
				res = fun.resd
			except AttributeError:
				res = (y - f(y, *fun.pars)) / dy
			resd_ax.plot(x[mask] * self.reX, res, **resdkw)
