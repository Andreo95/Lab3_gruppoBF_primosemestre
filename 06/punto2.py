import sys, os
sys.path.append( os.path.join(os.path.realpath('..'), '00 - Risorse', 'Python') )
folder = os.path.realpath('.')
import numpy as np, pylab, matplotlib.pyplot as plt, matplotlib as mpl, scipy.stats.distributions as dists
from lab import *
from pyan import *
#from solvers import *

plt.close('all')
Fitter._fitter_func = staticmethod(fit_generic_xyerr)


#### Parte 1

datafile = 'dati_1.txt'
rawdata = np.loadtxt(os.path.join(folder, 'Dati', datafile)).T

R1 = 2.27e3	#kohm
R2 = 22.1e3	#kohm

vin = rawdata[0] - rawdata[1]
vout = rawdata[2] - rawdata[3]
dvin, dvout = mme([vin,vout], 'volt', 'oscil')

linamp = createline()
linamp.mask = vin < 1.1


one = Fitter(vin, vout, dvin, dvout)
one.fit(linamp)


first = Graph.from_fitter(one)
first.labelY = "$V_{out}$ [V]"
first.labelX = "$V_{in}$ [V]"
first.title = "Amplificatore invertente"
first.draw(linamp, resid=True)


#### Parte 2

datafile = 'dati_2.txt'
rawdata = np.loadtxt(os.path.join(folder, 'Dati', datafile)).T


freq = rawdata[0]*1e3
dfreq = freq/100
Vout = rawdata[1] - rawdata[2]
dVout = mme(Vout, 'volt', 'oscil')
Vin = .528 - (-.512) # Volt
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
#print(sum((d[1] - amplicut(d[0], *amplicut.pars))**2 /( d[3]**2 + d[2]**2*amplicut.deriv(d[0], *amplicut.pars)**2)))
tell_chi2((d[1] - amplicut(d[0], *amplicut.pars))**2 /( d[3]**2 + d[2]**2*amplicut.deriv(d[0], *amplicut.pars)**2), 4)

second = Graph.from_fitter(two)
second.typeX = 'log'
second.labelX = "Frequenza [Hz]"
second.labelY = "Guadagno [dB]"
second.title = "Frequenza di taglio dell'opamp"
second.draw(amplicut, flattop, resid=True)


#### Parte 3

datafile = 'dati_3.txt'
rawdata = np.loadtxt(os.path.join(folder, 'Dati', datafile)).T

res = rawdata[0]*1e3
dres = mme(res, 'ohm')
Vin = .472 - (-.568)    # Volt
Vout = rawdata[1]
dVout = mme(Vout, 'volt', 'oscil')
ft = rawdata[2]*1e3
dft = ft/20
Av = (Vout/Vin)
dAv = (dVout/Vin)


GBW  = Av*ft
dGBW = np.sqrt(dAv**2 * ft**2 + dft**2 * Av**2)

three = Fitter(res, GBW, dres, dGBW)

third = Graph.from_fitter(three)
third.reY = 1e-6
third.reX = 1e-3
third.labelX = "Resistenza [k$\Omega$]"
third.labelY = "GBW [MHz]"
third.title = ""
third.draw(resid=True)


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
phase = 2*1*t*freq
dphase = 2*1*np.sqrt(t**2 * dfreq**2)
# dphase[dphase < 1e-3] = 2*1*np.sqrt(dt**2 * freq**2)[dphase < 1e-3]

def sfasinteg(f, ft):
    return -np.arctan(f/ft)/np.pi

sfasinteg.mask = phase != 0
sfasinteg.deriv = lambda f, ft: +1/(1 + (f/ft)**2)/np.pi

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


#### Parte 5

datafile = 'dati_5.txt'
rawdata = np.loadtxt(os.path.join(folder, 'Dati', datafile)).T
sorting = np.argsort(rawdata[0])
rawdata = rawdata[:,sorting]

freq = rawdata[0]
dfreq = freq/100
Vout = rawdata[1] - rawdata[2]
Vin = 2.08      # Volt
Av = 20*np.log10(Vout/Vin)
dAv = 20*np.log10(np.e)*(mme(Vout, 'time', 'oscil')/Vout)


t = rawdata[3]*1e-6
dt = mme(t, 'volt', 'oscil')
phase = 2*1*(t*freq)
dphase = 2*1*np.sqrt(t**2 * dfreq**2)
# dphase[dphase < 1e-3 ] = 2*1*np.sqrt(dt**2 * freq**2)[dphase < 1e-3]

def sfasderiv(f, ft):
    return -np.arctan(ft/f)/np.pi
    
sfasderiv.mask = freq < 3e3
sfasderiv.deriv = lambda f, ft: -1/(1 + (ft/f)**2)/np.pi


highpass = createline('log')
highpass.mask = (freq < 860)
highpass.bounds = [(0, 1e4)]

flattop = createline('const')
flattop.mask = (freq < 5e4) & (freq > 13e3)
flattop.bounds = [(100, 1e10)]

five = Fitter(freq, Av, dfreq, dAv)
five.fit(highpass, flattop)

five_phase = Fitter(freq, phase, dfreq, dphase)
five_phase.fit(sfasderiv, max_cycles=30)

fifth = Graph.from_fitter(five)
fifth.typeX = 'log'
fifth.labelX = "Frequenza [Hz]"
fifth.labelY = "Guadagno [dB]"
fifth.title = ""
fifth.draw(highpass, flattop, resid=True)

fifth_phase = Graph(freq, phase, dfreq, dphase)
fifth_phase.typeX = 'log'
fifth_phase.labelX = "Frequenza [Hz]"
fifth_phase.labelY = "$\Delta\phi$ [$\pi$ rad]"
fifth_phase.title = ""
fifth_phase.draw(sfasderiv, resid=True)


###### ending
plt.draw_all()
plt.show()
Graph.countfigs.send(0)
