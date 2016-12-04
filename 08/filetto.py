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

## cose

def newsine():
    def sinewave(t, w, phi, V0, C):
        return V0*np.sin(w*t+phi) + C
    sinewave.deriv = lambda t, w, phi, V0, C: V0*w*np.cos(w*t+phi)
    return sinewave

myreg = re.compile("parte1_(?P<freq>[0-9]+).csv")
csvdata = os.path.join(folder, "Dati", "parte1")

R_1 = 9.95e3 #kohm
R_2 = 9.94e3#0.94e3
R_3 = 9.90e3
R_4 = 9.94e3
R_5 = 9.90e3
C_1 = 10.81e-9 #nF
C_2 = 10.11e-9 #nF

def ampligain(pot, x, rdiodes='ignore'):
    # guadagno dell'opamp come ampli non invertente
    if rdiodes == 'ignore':
        R_d = R_3
    else:
        R_d = rdiodes
    return 1 + (R_4 + R_d + (1-x)*pot) / (R_5 + x*pot)


def beth(f):
    # attenuazione della rete di feedback (di wien)
    # OCCHIO: è complessa, ci sarà da prenderne modulo e argomento...
    w = 2*np.pi*f
    ws = 1/(R_1 * C_1)
    wp = 1/(R_2 * C_2)
    invb = 1 + R_1/R_2 * (1 + ws/wp + 1j*(w/wp - ws/w))
    return 1/invb

def par(R1, R2):
    return (R1*R2)/(R1+R2)

def beth1(f):
    w=2*np.pi*f
    RC_1=1j/(C_1*w)
    RC_2=1j/(C_2*w)
    para=par(R_2, RC_2)
    return para/(para+RC_1+R_1)


def agg(x):
    return x+(2*(x<0)-1)*np.pi


## lettura

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
    
    
    Av=vout.pars[2]/vin.pars[2]
    dAv=((vout.sigmas[2]/vout.pars[2])**2+ vin.sigmas[2]**2 / vin.pars[2]**2)**0.5 
    avgf =(vout.pars[0]/vout.sigmas[0]**2 + vin.pars[0]/vin.sigmas[0]**2)/(vout.sigmas[0]**-2 + vin.sigmas[0]**-2)/(2*np.pi)
    df = 1 / np.sqrt(vout.sigmas[0]**-2 + vin.sigmas[0]**-2) / (2*np.pi)

    results.append((avgf, df, phase, dphase, Av, dAv))

results = np.array(results).T

freqs, dfreqs = results[0:2]
phase, dphase = results[2:4]
phase=agg(phase)
gain, dgain = results[4:6]
dgain=dgain#+gain*(1.44/100)
dgaincal=dgain+gain*(2**0.5/100) #1/100 errore calibrazione oscilloscopio? Penso che le due tracce siano acquisite in maniera indipendente.... 
## graphing

sfasamento = Graph(freqs, phase, dfreqs, dphase)
sfasamento.typeX = 'log'
sfasamento.title="Frequenza vs fase"
aperbeta = Graph(freqs, gain, dfreqs, dgain)
aperbeta.typeX = 'log'

aperbeta.title="Frequenza vs guadagno"

sfasamento.draw()
aperbeta.draw()

#plt.show()

## fits senza calibrazione

def amplificazione(f, A):
    x=beth(f)*A   #ampligain(pot, x, diodes)
    return np.absolute(x)

amplificazione.pars=[3.0]
fitt = Fitter(freqs, gain, dfreqs, dgain) 
fitt.fit(amplificazione)
terzo = Graph.from_fitter(fitt)
terzo.typeX = 'log'
terzo.title="Fit A ($Ab(f)$) senza errori di calibrazione"
terzo.labelX="frequenza [Hz]"
terzo.labelY="$A\beta(f)$"
terzo.draw(amplificazione)

##fit con calibrazione
amplificazione.pars=[3.0]
fitt = Fitter(freqs, gain, dfreqs, dgaincal) 
fitt.fit(amplificazione)
terzo = Graph.from_fitter(fitt)
terzo.typeX = 'log'
terzo.title="Fit A ($Ab(f)$)"
terzo.labelX="frequenza [Hz]"
terzo.labelY="$A\beta(f)$"
terzo.draw(amplificazione)

print("si vede come il chiq torni tantissimo!")

###overfitting?...
def amplificazione2(f, pot , x, diodes):
    x=beth(f)*ampligain(pot, x, diodes)
    return np.absolute(x)

amplificazione2.pars=[10000, 0.1, 4000]
fitt = Fitter(freqs, gain, dfreqs, dgain) 
fitt.fit(amplificazione2)
quarto = Graph.from_fitter(fitt)
quarto.typeX = 'log'
quarto.title="Fit"
quarto.draw(amplificazione2)



print("chiaramente non cambia nulla fra i due modelli, forse è meglio mettere direttamente la A")


def fase(f, g):
    print(g)
    return np.arctan(np.imag(beth(f)), np.real(beth(f)))

quinto=Graph(freqs, phase, dfreqs, dphase)
fase.pars=[1]
quinto.draw(fase)


plt.show()


# sorter = np.argsort(freqs)
# freqs, phase, gain = np.array([freqs, phase, gain])[:,sorter]
