import sys, os
sys.path.append(os.path.join(os.path.realpath('..'), '00 - Risorse', 'Python'))
folder = os.path.realpath('.')
from pyan import *
import Oscillografo
import numpy as np
from lab import *
import matplotlib.pyplot as plt
Fitter._fitter_func = staticmethod(fit_generic_xyerr)


folder = os.path.realpath('.')
datafile="Figs-Tabs\slew_1.csv"
path=os.path.join(folder, datafile)
o=Oscillografo.OscilloscopeData(path)
o.plot()

f=Fitter(o.T2 ,o.CH2, np.ones(len(o.CH2)), np.ones(len(o.CH2)))
fun=createline("const")
fun.mask = (o.T2 < -1.6e-6)# o.T2>5e-7
f.fit(fun)
chiq=sum(fun.resd**2)
mult=(chiq/len(fun.resd))**0.5
print("mult={0}".format(mult))



f=Fitter(o.T2 ,o.CH2, np.ones(len(o.CH2))*(np.amax(o.T2)-np.amin(o.T2))/5000 , np.ones(len(o.CH2))*mult)
fun=createline()
fun.mask = (o.T2 < 0) & (o.T2 > -1.2e-6)
f.fit(fun)
#chiq=sum(fun.resd**2/(np.ones(len(fun.mask))*mult)**2)

g=Graph.from_fitter(f)
g.title="Slow ratio e discriminatore"
g.labelX="tempo [s]"
g.labelY="tensione [V]"
g.draw(fun, resid=True)


plt.show()
