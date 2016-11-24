import sys, os
sys.path.append( os.path.join(os.path.realpath('..'), '00 - Risorse', 'Python') )
folder = os.path.realpath('.')
import numpy as np, pylab, matplotlib.pyplot as plt, matplotlib as mpl, scipy.stats.distributions as dists
from lab import *
from pyan import *

plt.close('all')
rawdata = None

#### Parte 000

datafile = 'dati_?.txt'
rawdata = np.loadtxt(os.path.join(folder, 'Dati', datafile)).T



fit_one = Fitter(x, y, dx, dy)
fit_one.fit()

first = Graph.from_fitter(fit_one)
first.set_edges()
first.draw('data', fitfun, line)


###### ending
# plt.draw_all()
plt.show()
Graph.countfigs.send(0)
