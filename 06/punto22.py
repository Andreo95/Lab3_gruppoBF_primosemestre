#stima non veloce incertezze di calibrazione...ancora da finire

import sys, os
sys.path.append( os.path.join(os.path.realpath('..'), '00 - Risorse', 'Python') )
folder = os.path.realpath('.')
import numpy as np, pylab, matplotlib.pyplot as plt, matplotlib as mpl, scipy.stats.distributions as dists
from lab import *
from pyan import *
import numpy.random
#from solvers import *

plt.close('all')
Fitter._fitter_func = staticmethod(fit_generic_xyerr)

datafile = 'dati_2.txt'
rawdata = np.loadtxt(os.path.join(folder, 'Dati', datafile)).T

results=[]
for i in range(100):
    freq = rawdata[0]*1e3
    dfreq = freq/100
    dIn=(3/100)*np.random.standard_normal()
    dOut=(3/100)*np.random.standard_normal()
    Vout = (rawdata[1] - rawdata[2])*(1+dOut)
    dVout = mme(Vout, 'volt', 'oscil')
    Vin = (.528 - (-.512))*(1+dIn) # Volt
    Av = 20*np.log10(Vout/Vin)
    dAv = 20*np.log10(np.e)*(dVout/Vout)
    
    amplicut = createline('log')
    amplicut.mask = freq > 4e5
    amplicut.bounds = [(1e5, 1e10)]
    
    flattop = createline('const')
    flattop.mask = freq < 4e4
    
    two = Fitter(freq, Av, dfreq, dAv)
    two.fit(amplicut, flattop)
    d=two.data[:,amplicut.mask]
    
    results.append((amplicut.pars,tell_chi2((d[1] - amplicut(d[0], *amplicut.pars))**2 /( d[3]**2 + d[2]**2*amplicut.deriv(d[0], *amplicut.pars)**2), 4)))
    
    # second = Graph.from_fitter(two)
    # second.typeX = 'log'
    # second.labelX = "Frequenza [Hz]"
    # second.labelY = "Guadagno [dB]"
    # second.title = "Frequenza di taglio dell'opamp"
    # second.draw(amplicut, flattop, resid=True)

print(results)

plt.plot()