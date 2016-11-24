import sys, os
sys.path.append( os.path.join(os.path.realpath('..'), '00 - Risorse', 'Python') )
folder = os.path.realpath('.')
import numpy as np, pylab, matplotlib.pyplot as plt, matplotlib as mpl, scipy.stats.distributions as dists
from lab import *
#from pyan import *
import Oscillografo
import scipy
dati=Oscillografo.OscilloscopeData(os.path.join(folder, "Figs-Tabs", "slewrate.csv"))


results=[]
domain=np.linspace(1.5e-7, 3.0e-7, 100)


pylab.title("$V_{out}(t)$ con diversi cut_off")
pylab.xlabel("tempo [s]")
pylab.ylabel("$V_{out}$")
for cut_off in domain:
    espera=lambda t, Va, Vb, tau: Va-Vb*np.exp(-t/tau)
    a, b = scipy.optimize.curve_fit(espera, dati.T2[dati.T2>cut_off], dati.CH2[dati.T2>cut_off], p0=(5, 10,1e-6), sigma=mme(10, "volt", "oscil"))
    print(a, b)
    
    pylab.plot(dati.T2[dati.T2>cut_off],espera(dati.T2[dati.T2>cut_off], *a))
    pylab.errorbar(dati.T2[dati.T2>cut_off], dati.CH2[dati.T2>cut_off], mme(10, "volt"))
    resd=(espera(dati.T2[dati.T2>cut_off], *a)-dati.CH2[dati.T2>cut_off])/mme(10, "volt")
    chisq=sum(resd**2)
    dof=len(resd)
    pval = dists.chi2.sf(chisq, dof)
    results.append(pval)
    

print(results)
pylab.figure(2)
pylab.title("cut off vs p-value")
pylab.xlabel("tempo dopo il triggering [s]")
pylab.ylabel("p-value")
pylab.plot(domain, results)

cut_off=1.95e-7
retta=lambda x, m, q: m*x+q
newT=dati.T2[(dati.T2>1.5e-7) &(dati.T2<cut_off)]
newY=dati.CH2[(dati.T2>1.5e-7) &(dati.T2<cut_off)]
a, b=scipy.optimize.curve_fit(retta, newT, newY, p0=None, sigma=mme(10, "volt", 'oscil'))
pylab.figure(3)
pylab.title("retta")
pylab.xlabel("tempo [s]")
pylab.ylabel("$V_{out}$ [V]")
pylab.errorbar(newT, newY, mme(10, "volt", 'oscil'))
pylab.plot(np.linspace(1.5e-7, cut_off, 1000), retta(np.linspace(1.5e-7, cut_off, 1000), *a))

resd=(retta(newT, *a)-newY)/mme(10, "volt", 'oscil')
chisq=sum(resd**2)
dof=len(resd)
pval = dists.chi2.sf(chisq, dof)
print(a, b, chisq, dof, pval)
pylab.show()

