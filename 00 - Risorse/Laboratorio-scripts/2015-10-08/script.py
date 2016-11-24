from pylab import *
from scipy.optimize import curve_fit

v, dv, i, di = loadtxt('data10.txt', unpack=True)

def ohm_eq(v, r, i0):
	return i0 + v / r

title('fake xperiment')
xlabel('potential / V')
ylabel('current / mA')
xlim(min(v) - 0.2, max(v) + 0.2)

errorbar(v, i, yerr=di, xerr=dv, fmt='.', label='data')

savefig('plot-raw.pdf')

par, cov = curve_fit(ohm_eq, v, i, p0=(300.0, 0.002), sigma=di)

r, i0 = par
dr, di0 = sqrt(diag(cov))

c = cov[0, 1] / (dr * di0)

print("r = %.1f +/- %.1f (%.1f %%) [kohm]" % (r, dr, dr / r * 100))
print("i0 = %.5f +/- %.5f (%.1f %%) [mA]" % (i0, di0, di0 / i0 * 100))
print("c = %.2f" % (c))
print(cov)

vv = linspace(min(v), max(v), 1000)

plot(vv, ohm_eq(vv, r, i0), label='fit')

legend(loc=0)

savefig('plot-fit.pdf')

show()

