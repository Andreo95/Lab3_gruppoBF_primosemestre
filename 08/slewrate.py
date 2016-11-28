import sys, os
sys.path.append(os.path.join(os.path.realpath('..'), '00 - Risorse', 'Python'))
folder = os.path.realpath('.')
from pyan import *
import Oscillografo
import numpy as np
from lab import *
import matplotlib.pyplot as plt
import scipy.stats.distributions as dists
from os import listdir
import re


Fitter._fitter_func = staticmethod(fit_generic_xyerr)
reg = re.compile("parte1_(?P<freq>[0-9]+).csv")

folder = os.path.realpath('.')
path=os.path.join(folder ,"Dati\\parte1")

results1=[]
results2=[]

for f in listdir(path):
    w = reg.match(f)
    print(w.group('freq'))
    freq=float(w.group("freq"))
    o=Oscillografo.OscilloscopeData(os.path.join(path, f))
    
    
    f=Fitter(o.T2 ,o.CH2, np.ones(len(o.CH2))*(np.amax(o.T2)-np.amin(o.T2))/5000 ,np.ones(len(o.CH2))*o.dCH2)
    fun=lambda t, w, phi, V0, C: V0*np.sin(w*t+phi) + C
    fun.pars = [2*np.pi*freq, 0, .2, 0]
    fun.deriv=lambda t, w, phi, V0, C: V0*w*np.cos(w*t+phi)
    f.fit(fun)
    results1.append(fun.pars)    
    f=Fitter(o.T1 ,o.CH1, np.ones(len(o.CH1))*(np.amax(o.T1)-np.amin(o.T1))/5000 ,np.ones(len(o.CH2))*o.dCH1)
    f.fit(fun)
    results2.append(fun.pars)
    
    #g=Graph.from_fitter(f)
    #g.draw(fun, resid=True)
    
    
    #plt.show()

#plt.show()