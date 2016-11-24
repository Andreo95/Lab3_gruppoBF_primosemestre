import sys, os
sys.path.append( os.path.join(os.path.realpath('..'), '00 - Risorse', 'Python') )
folder = os.path.realpath('.')
import numpy as np, pylab, matplotlib.pyplot as plt, matplotlib as mpl, scipy.stats.distributions as dists
from lab import *
from pyan import *
#from solvers import *

plt.close('all')
Fitter._fitter_func = staticmethod(fit_generic_xyerr)

#### Parte 4

datafile = 'dati_4.txt'
rawdata = np.loadtxt(os.path.join(folder, 'Dati', datafile)).T

freq = rawdata[0]*1e3
dfreq = freq/100
Vout = rawdata[1] - rawdata[2]
Vin = 2.08      # Volt
Av = 20*np.log10(Vout/Vin)
dAv = 20*np.log10(np.e)*(mme(Vout, 'time', 'oscil')/Vout)

t = rawdata[3]*1e-6
dt = mme(t, 'volt', 'oscil')
phase = -2*1*t*freq
dphase = 2*1*np.sqrt(t**2 * dfreq**2)
dphase[dphase < phase/1e3] += 2*1*np.sqrt(dt**2 * freq**2)[dphase < phase/1e3]

def sfasinteg(f, ft):
    return np.arctan(-f/ft)/np.pi

# sfasinteg.mask = phase != 0
sfasinteg.deriv = lambda f, ft: +1/(1 + (f/ft)**2)/np.pi/ft

lowpass = createline('log')
lowpass.mask = (freq > 2e3)
lowpass.bounds = [(100, 1e10)]

flattop = createline('const')
flattop.mask = (freq < 50)
flattop.bounds = [(0, 5e3)]

four = Fitter(freq, Av, dfreq, dAv)
four.fit(lowpass, flattop)

four_phase = Fitter(freq, phase, dfreq, dphase)
four_phase.fit(sfasinteg, max_cycles = 30)

fourth = Graph.from_fitter(four)
fourth.typeX = 'log'
fourth.labelX = "Frequenza [Hz]"
fourth.labelY = "Guadagno [dB]"
fourth.title = ""
fourth.draw(lowpass, flattop, resid=True)

fourth_phase = Graph(freq, phase, dfreq, dphase)
fourth_phase.typeX = 'log'
fourth_phase.labelX = "Frequenza [Hz]"
fourth_phase.labelY = "$\Delta\phi$ [$\pi$ rad]"
fourth_phase.title = ""
fourth_phase.draw(sfasinteg, resid=True)
plt.draw_all()
plt.show()
Graph.countfigs.send(0)
