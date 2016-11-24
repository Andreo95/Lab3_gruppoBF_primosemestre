from pylab import *
from scipy.optimize import curve_fit

def ohm_eq(v, r, i0):
	return i0 + v / r

v, dv, i, di = loadtxt('data10.txt', unpack=True)

sidi2 = sum(i / di**2)
sv2di2 = sum(v**2 / di**2)
svdi2 = sum(v / di**2)
sivdi2 = sum(i * v / di**2)
s1di2 = sum(1 / di**2)

denom = s1di2 * sv2di2 - svdi2**2
i0 = (sidi2 * sv2di2 - svdi2 * sivdi2) / denom
invr = (s1di2 * sivdi2 - sidi2 * svdi2) / denom
di0 = sqrt(sv2di2 / denom)
dinvr = sqrt(s1di2 / denom)
r = 1 / invr
dr = dinvr / invr**2

print("r = %.1f +/- %.1f (%.1f %%) [kohm]" % (r, dr, dr / r * 100))
print("i0 = %.5f +/- %.5f (%.1f %%) [mA]" % (i0, di0, di0 / i0 * 100))
print('c = 1')

vx = 2.0 # V
ip = ohm_eq(vx, r, i0)
C = dinvr * di0
dip = sqrt(di0**2 + (dinvr * vx)**2 + 2 * C * vx)

print("i' = %.5f +/- %.5f (%.1f %%) [mA]" % (ip, dip, dip / ip * 100))
