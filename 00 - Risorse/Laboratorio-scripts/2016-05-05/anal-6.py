from pylab import *
from scipy.optimize import curve_fit

v1 = 0.800
dv1 = 0.032
r2, v2, dv2, dr2 = loadtxt('data-6.txt', unpack=True)

p2 = v2**2 / (2*r2)
dp2 = p2 * sqrt(2)*(dv2/v2)

#par, cov = curve_fit(fun, v1, v2, sigma=dv2, p0=(0))

#sigmaeff = sqrt(dv2**2 + (dv1*par[0])**2)

#par, cov = curve_fit(fun, v1, v2, sigma=dv2, p0=par)

#sigmaeff = sqrt(dv2**2 + (dv1*par[0])**2)

#par, cov = curve_fit(fun, v1, v2, sigma=dv2, p0=par)

#tv = par[0]
#dtv = sqrt(cov[0])

#print("tv = %f +- %f" % (tv, dtv))

title('Power transfer')
xscale('log')
yscale('log')
xlabel('resistive load [ohm]')
ylabel('power used by load [watt]')
errorbar(r2, p2, yerr=dp2, xerr=dr2, fmt=',', capsize=0)
savefig('plot-6.pdf')
show()
