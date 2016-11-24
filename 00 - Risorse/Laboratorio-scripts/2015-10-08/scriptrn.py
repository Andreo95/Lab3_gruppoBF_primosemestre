from pylab import *
from scipy.optimize import curve_fit

v, dv, i, di = loadtxt('data10.txt', unpack=True)

def ohm_eq(v, r, i0):
	return i0 + v / r

par, cov = curve_fit(ohm_eq, v, i, p0=(300.0, 0.002), sigma=di)

r, i0 = par
dr, di0 = sqrt(diag(cov))

c = cov[0, 1] / (dr * di0)

rn = (i - ohm_eq(v, r, i0)) / di

figure(1)
title('fake xperiment: residuals')
xlabel('potential / V')
ylabel('normalized residuals')
xlim(min(v) - 0.2, max(v) + 0.2)
plot(v, rn, 'o')
grid()
savefig('plot-rn.pdf')

figure(2)
title('fake xperiment: residuals')
xlabel('normalized residuals')
ylabel('occurences')
hist(rn, bins=4)
savefig('plot-rnhist.pdf')

show()

