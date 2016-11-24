from pylab import *
from scipy.optimize import curve_fit
from scipy.odr import Model, Data, RealData, ODR

v, dv, i, di = loadtxt('data10.txt', unpack=True)

def ohm_eq(par, v):
	return par[1] + v / par[0]

title('fake xperiment')
xlabel('potential / V')
ylabel('current / mA')
xlim(min(v) - 0.2, max(v) + 0.2)

errorbar(v, i, yerr=di, xerr=dv, fmt='.', label='data')

linear = Model(ohm_eq)
mydata = RealData(v, i, sx=dv, sy=di)
myodr = ODR(mydata, linear, beta0=(300.0, 0.002))
myoutput = myodr.run()

r, i0 = myoutput.beta
dr, di0 = myoutput.sd_beta
cov = myoutput.cov_beta

c = cov[0, 1] / sqrt(cov[0, 0] * cov[1, 1])

print("r = %.1f +/- %.1f (%.1f %%) [kohm]" % (r, dr, dr / r * 100))
print("i0 = %.5f +/- %.5f (%.1f %%) [mA]" % (i0, di0, di0 / i0 * 100))
print("c = %.2f" % (c))
print(cov)

vv = linspace(min(v), max(v), 1000)

plot(vv, ohm_eq((r, i0), vv), label='fit')

legend(loc=0)

savefig('plot-fitodr.pdf')

show()

