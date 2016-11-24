import Oscillografo
import os
from pyan import *
import numpy as np
from lab import *
import matplotlib.pyplot as plt
Fitter._fitter_func = staticmethod(fit_generic_xyerr)


folder = os.path.realpath('.')
datafile="Figs-Tabs\slew_1.csv"
path=os.path.join(folder, datafile)
o=Oscillografo.OscilloscopeData(path)

f=Fitter(o.T2 ,o.CH2, mme(o.T2 , "time", "oscil"), np.ones(len(o.CH2))*o.dCH2)
fun=createline()
fun.mask = (o.T2 < 0) & (o.T2 > -1.2e-6)
f.fit(fun)

Graph.from_fitter(f).draw(fun, resid=True)

plt.show()



