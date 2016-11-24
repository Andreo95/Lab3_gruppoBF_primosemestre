from pylab import *
from scipy.optimize import curve_fit

v1, v2, dv1, dv2 = loadtxt('data-1.txt', unpack=True)

fun = lambda v1, tv: tv * v1

par, cov = curve_fit(fun, v1, v2, sigma=dv2, p0=(0))

sigmaeff = sqrt(dv2**2 + (dv1*par[0])**2)

par, cov = curve_fit(fun, v1, v2, sigma=dv2, p0=par)

sigmaeff = sqrt(dv2**2 + (dv1*par[0])**2)

par, cov = curve_fit(fun, v1, v2, sigma=dv2, p0=par)

tv = par[0]
dtv = sqrt(cov[0])

print("tv = %f +- %f" % (tv, dtv))

errorbar(v1, v2, yerr=dv2, xerr=dv1, fmt=',', capsize=0, label='Data')
l = linspace(min(v1), max(v1), 256)
plot(l, fun(l, tv), label='Fit')
legend(loc=0)
show()

