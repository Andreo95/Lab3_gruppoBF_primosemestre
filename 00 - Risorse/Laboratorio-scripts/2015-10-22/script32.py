from pylab import *

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

chi2 = sum((iVl - (a * iRl + b)) ** 2 / diVl ** 2)

print('linearized analytical fit:')
print("Rg = %.3g +- %.3g" % (Rg, dRg))
print("V0 = %.3g +- %.3g" % (V0, dV0))
print("chi2/dof = %.3g/%d" % (chi2, len(Rl) - 2))

figure(3)
errorbar(iRl, iVl, yerr=diVl, xerr=diRl, label='Data', fmt='.')
spaceiRl = linspace(min(iRl), max(iRl), 1000)
plot(spaceiRl, a * spaceiRl + b, label='Fit')
title('Tension versus load resistance (linearized)')
xscale('log')
#yscale('log')
xlabel('1 / R [1/ohm]')
ylabel('1 / V [1/volt]')
legend(loc=0)
savefig('script3fig3.pdf')

show()

