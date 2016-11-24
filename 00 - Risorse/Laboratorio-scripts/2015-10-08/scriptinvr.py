from pylab import *
from scipy.optimize import curve_fit

v, dv, i, di = loadtxt('data10.txt', unpack=True)

def ohm_eq(v, invr, i0):
	return i0 + v * invr

errorbar(v, i, yerr=di, xerr=dv, fmt='.')

savefig('plot-raw-invr.pdf')

par, cov = curve_fit(ohm_eq, v, i, p0=(0.003, 0.002), sigma=di)

invr, i0 = par
dinvr, di0 = sqrt(diag(cov))

c = cov[0, 1] / (dinvr * di0)

print("1/r = %.6f +/- %.6f" % (invr, dinvr))
print("i0 = %.6f +/- %.6f" % (i0, di0))
print("c = %f" % (c))

vv = linspace(min(v), max(v), 1000)

plot(vv, ohm_eq(vv, invr, i0))

savefig('plot-fit-invr.pdf')

show()
