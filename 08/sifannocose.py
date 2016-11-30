import sys, os
sys.path.append(os.path.join(os.path.realpath('..'), '00 - Risorse', 'Python'))
folder = os.path.realpath('.')
import numpy as np
from lab import *
from pyan import *
import Oscillografo
import matplotlib.pyplot as plt
import scipy.stats.distributions as dists
import re

Fitter._fitter_func = staticmethod(fit_generic_xyerr)

def newsine():
	def sinewave(t, w, phi, V0, C):
		return V0*np.sin(w*t+phi) + C
	sinewave.deriv = lambda t, w, phi, V0, C: V0*w*np.cos(w*t+phi)
	return sinewave

myreg = re.compile("parte1_(?P<freq>[0-9]+).csv")

csvdata = os.path.join(folder, "Dati", "parte1")

results1 = []
results2 = []

for filename in os.listdir(csvdata):
	w = myreg.match(filename)
	freq = float(w.group("freq"))
	o = Oscillografo.OscilloscopeData(os.path.join(csvdata, filename))

	f = Fitter(o.T2, o.CH2, np.ones(len(o.CH2))*mme(np.amax(o.T2)-np.amin(o.T2), 'time', 'oscil'), np.ones(len(o.CH2))*o.dCH2)
	fun = newsine()
	fun.pars = [2*np.pi*freq, 0, 1, 0]
	f.fit(fun, verbose=False)
	results1.append(fun.pars)
	f = Fitter(o.T1, o.CH1, np.ones(len(o.CH1))*mme(np.amax(o.T1)-np.amin(o.T1), 'time', 'oscil'), np.ones(len(o.CH2))*o.dCH1)
	f.fit(fun, verbose=False)
	results2.append(fun.pars)
	deltaphi = (results2[-1][1] - results1[-1][1])

	print(freq, 'Hz:\nCH1:', np.amax(o.CH1), np.amin(o.CH1), '\nCH2:', np.amax(o.CH2), np.amax(o.CH2), '\ndeltaphi:', deltaphi)

	#g = Graph.from_fitter(f)
	#g.draw(fun, resid = True)


	#plt.show()

#plt.show()