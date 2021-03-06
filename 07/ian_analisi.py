import sys, os
sys.path.append(os.path.join(os.path.realpath('..'), '00 - Risorse', 'Python'))
folder = os.path.realpath('.')
import numpy as np, pylab, matplotlib.pyplot as plt, matplotlib as mpl, scipy.stats.distributions as dists
from lab import *
from pyan import *

plt.close('all')

@staticmethod
def _fitfun(f, df, *args, **kwargs):
	return fit_generic_xyerr2(f, *args, **kwargs)

Fitter._fitter_func = _fitfun #staticmethod(fit_generic_xyerr)

#### Parte 000

datafile = 'dati_B_2.txt'
rawdata = np.loadtxt(os.path.join(folder, 'Dati', datafile)).T

vin = rawdata[0] - rawdata[1]
dvin = mme(vin, 'volt', 'oscil')
tot = 1e-6 * (rawdata[2]) # + rawdata[3])/2
dt = mme(tot, 'time', 'oscil')

c_t = 1.059e-9 #nF
R_2 = 99.68 #ohm
R_1 = 98.5e3 #kohm
c_f = 1.036e-9 #nF
tau = R_1 * c_f

charge = c_t * vin
dq = np.sqrt(c_t**2 * dvin**2 + mme(c_t, 'farad')**2 * vin)
vsh = charge / c_f

def q_rilevata(tot, vth):
	return c_f * vth * np.exp(tot / tau)

q_rilevata.pars = [.2]
q_rilevata.deriv = lambda tot, vth: 1/(tau) * q_rilevata(tot, vth)
# q_rilevata.mask = tot > 1e-4

fit_one = Fitter(tot, charge, dt, dq)
fit_one.fit(q_rilevata)

first = Graph.from_fitter(fit_one)
first.title = 'Rilevazione di carica'
first.labelX = 'Time over Threshold [ms]'
first.labelY = 'Carica immessa [nC]'
first.reX = 1e3
first.reY = 1e9
# first.typeY = 'log'
first.draw(q_rilevata, resid=True)
print(q_rilevata.cov)


###### ending
# plt.draw_all()
plt.show()
Graph.countfigs.send(0)
