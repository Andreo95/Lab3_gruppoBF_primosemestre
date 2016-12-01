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

results = []

for filename in os.listdir(csvdata):
	w = myreg.match(filename)
	freq = float(w.group("freq"))
	o = Oscillografo.OscilloscopeData(os.path.join(csvdata, filename))

	f = Fitter(o.T2, o.CH2, np.ones(len(o.CH2))*mme(np.amax(o.T2)-np.amin(o.T2), 'time', 'oscil'), np.ones(len(o.CH2))*o.dCH2)
	vout = newsine()
	vout.pars = [2*np.pi*freq, 0, 1, 0]
	f.fit(vout, verbose=False)
	if vout.pars[2] < 0:
		vout.pars[2] *= -1
		if vout.pars[1] > 0:
			vout.pars[1] -= np.pi
		else:
			vout.pars[1] += np.pi

	f = Fitter(o.T1, o.CH1, np.ones(len(o.CH1))*mme(np.amax(o.T1)-np.amin(o.T1), 'time', 'oscil'), np.ones(len(o.CH2))*o.dCH1)
	vin = newsine()
	vin.pars = [2*np.pi*freq, 0, 1, 0]
	f.fit(vin, verbose=False)
	if vin.pars[2] < 0:
		vin.pars[2] *= -1
		if vin.pars[1] > 0:
			vin.pars[1] -= np.pi
		else:
			vin.pars[1] += np.pi

	phase = vout.pars[1] - vin.pars[1]
	dphase = np.sqrt(vout.sigmas[1]**2 + vin.sigmas[1]**2)
	Av = 20 * np.log10(vout.pars[2] / vin.pars[2])
	dAv = 20 * np.log10(np.e) * np.sqrt(vout.sigmas[2]**2 / vout.pars[2]**2 + vin.sigmas[2]**2 / vin.pars[2]**2)
	avgf = (vout.pars[0]/vout.sigmas[0]**2 + vin.pars[0]/vin.sigmas[0]**2)/(vout.sigmas[0]**-2 + vin.sigmas[0]**-2)/(2*np.pi)
	df = 1 / np.sqrt(vout.sigmas[0]**-2 + vin.sigmas[0]**-2) / (2*np.pi)

	results.append((avgf, df, phase, dphase, Av, dAv))

results = np.array(results).T

freqs, dfreqs = results[0:2]
phase, dphase = results[2:4]
gain, dgain = results[4:6]

sfasamento = Graph(freqs, phase, dfreqs, dphase)
aperbeta = Graph(freqs, gain, dfreqs, dgain)

sfasamento.draw()
aperbeta.draw()

plt.show()

# sorter = np.argsort(freqs)
# freqs, phase, gain = np.array([freqs, phase, gain])[:,sorter]
