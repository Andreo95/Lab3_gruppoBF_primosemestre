import numpy as np, pylab, matplotlib.pyplot as plt, matplotlib as mpl, scipy.stats.distributions as dists
from lab import *
from ian import *
from PyTeX import *
from solvers import *
from scipy.optimize import curve_fit, brentq
from pynverse import inversefunc as invfcn


Rs = 4.70e6
Vse = 2.48
Vc = 1.184
Rin = Rs/(Vse/Vc-1)
w = 1020*2*np.pi
print(Rin)
C = (1/Rin**2-1/Rs**2)**0.5/w
print(C)
Vse = 2.54 #v
Vc = .520 #v
Rin = Rs/(Vse/Vc-1)
w = 10008*2*np.pi
print(Rin)
C = (1/Rin**2-1/Rs**2)**0.5/w
print(C)



def F(V2, Zc):
    return abs(V2*(1+Rs*(1/Rs+1/(-1j*Zc+227))))

popt, pcov = curve_fit(F,np.array([.520]), np.array([2.54]), p0=(1e7))
print(popt)

def vratio(freq, C):
    Rs = 4.70e6
    Ri = 4.57e6
    w = freq*2*np.pi
    Rp = 227
    Rd = 551
    Cin = 100e-9
    Zin = parallel(Ri, Rp + 1/1j/w/C) + 1/1j/w/Cin
    return 1 + abs(Rs/Zin)
    
# def freq1(C):
#     return vratio(1020, C)
#     
# def freq2(C):
#     return vratio(10080, C)

freqs = np.array((1020, 10080))
ratios = np.array((2.48/1.184, 2.54/.520))
errs = np.array((.08, .19))
popt, pcov = curve_fit(vratio, freqs, ratios, p0=(1e-12), sigma=errs, absolute_sigma=True)

# def test1(C):
#     return ratios[0] - freq1(C)
# 
# def test2(C):
#     return ratios[1] - freq2(C)
#     

qq = Graph(freqs, ratios, dYY=errs, funFacts = {vratio : dict(pars = popt, linetype = dict(color='red'))})
qq.compute()
qq.draw('data', vratio)