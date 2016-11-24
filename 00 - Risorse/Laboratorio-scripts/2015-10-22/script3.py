from pylab import *
from scipy.optimize import curve_fit

def fit_const_yerr(y, sigmay):
	from math import fsum
	from numpy import array
	dy2 = sigmay ** 2
	sydy2 = fsum(y / dy2)
	s1dy2 = fsum(1 / dy2)
	a = sydy2 / s1dy2
	vara = 1 / s1dy2
	return array([a, vara])

def fit_affine_yerr(x, y, sigmay):
	"""
		fit y = a * x + b
		
		Parameters
		----------
		x : M-length array
			independent data
		y : M-length array
			dependent data
		sigmay : M-length array
			standard deviation of y
		
		Returns
		-------
		a : float
			optimal value for a
		b : float
			optimal value for b
		vara : float
			variance of a
		varb : float
			variance of b
	"""
	from math import fsum
	from numpy import array
	dy2 = sigmay ** 2
	sy = fsum(y / dy2)
	sx2 = fsum(x ** 2 / dy2)
	sx = fsum(x / dy2)
	sxy = fsum(x * y / dy2)
	s1 = fsum(1 / dy2)
	denom = s1 * sx2 - sx ** 2
	a = (s1 * sxy - sy * sx) / denom
	b = (sy * sx2 - sx * sxy) / denom
	vara = s1 / denom
	varb = sx2 / denom
	return array([a, b, vara, varb])

# 1: calcola Rg noto V0 e una serie di Rl, Vl

V0 = 5.03
dV0 = V0 * 0.005 + 0.01

Rl, Vl = loadtxt('rlvl.txt', unpack=True)
dRl, dVl = zeros((2, len(Rl)))
for i in range(len(Rl)):
	if Rl[i] < 200:
		dRl[i] = Rl[i] * 0.008 + 0.3
	elif Rl[i] < 2000:
		dRl[i] = Rl[i] * 0.008 + 1
	else:
		dRl[i] = Rl[i] * 0.008 + 10
	dVl[i] = Vl[i] * 0.005 + 0.01

Rgs = Rl * (V0 / Vl - 1)
dRgs = sqrt((V0/Vl-1)**2 * dRl**2 + (Rl/Vl*dV0)**2 + (Rl*V0/Vl**2*dVl)**2)

Rg, dRg = fit_const_yerr(Rgs, dRgs)
dRg = sqrt(dRg)

print('compute Rgs separately and average:')
print("Rg = %.3g +- %.3g" % (Rg, dRg))

subplot(311)
errorbar(Rl, Rgs, yerr=dRgs, xerr=dRl, label='Data', fmt='.')
spaceRl = linspace(min(Rl), max(Rl), 1000)
plot(spaceRl, [Rg] * 1000, label='Fit')
title('Internal resistance versus load')
xscale('log')
xlabel('Rl')
ylabel('Rg')
legend(loc=0)
savefig('script3fig1.pdf')

# 2: calcola Rg e V0 fittando numericamente Vl(Rl) con gli errori su Vl

def fun_Vl(Rl, Rg, V0):
	return V0 / (1 + Rg/Rl)

par, cov = curve_fit(fun_Vl, Rl, Vl, p0=(20.0, 5.0), sigma=dVl)

Rg, V0 = par
dRg, dV0 = sqrt(diag(cov))
c = cov[0, 1] / sqrt(cov[0, 0] * cov[1, 1])

print('numeric fit:')
print("Rg = %.3g +- %.3g" % (Rg, dRg))
print("V0 = %.3g +- %.3g" % (V0, dV0))
print("c = %.2f" % (c))

subplot(312)
errorbar(Rl, Vl, yerr=dVl, xerr=dRl, label='Data', fmt='.')
plot(spaceRl, fun_Vl(spaceRl, *par), label='Fit')
title('Tension versus load')
xscale('log')
xlabel('Rl')
ylabel('Vl')
legend(loc=0)
savefig('script3fig2.pdf')

# 3: calcola Rg e V0 linearizzando e fittando analiticamente 1/Vl(1/Rl)

iVl = 1 / Vl
diVl = dVl / Vl**2
iRl = 1 / Rl
diRl = dRl / Rl**2

a, b, da, db = fit_affine_yerr(iRl, iVl, diVl)
da = sqrt(da)
db = sqrt(db)

Rg = a / b
dRg = Rg * sqrt((da/a)**2 + (db/b)**2)
V0 = 1 / b
dV0 = db / b**2

print('linearized analytical fit:')
print("Rg = %.3g +- %.3g" % (Rg, dRg))
print("V0 = %.3g +- %.3g" % (V0, dV0))

subplot(313)
errorbar(iRl, iVl, yerr=diVl, xerr=diRl, label='Data', fmt='.')
spaceiRl = linspace(min(iRl), max(iRl), 1000)
plot(spaceiRl, a * spaceiRl + b, label='Fit')
title('Tension versus load (linearized)')
xscale('log')
#yscale('log')
xlabel('1 / Rl')
ylabel('1 / Vl')
legend(loc=0)
savefig('script3fig3.pdf')

show()

