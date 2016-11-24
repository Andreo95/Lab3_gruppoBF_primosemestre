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

Rl, Vl, I = loadtxt('rlvli.txt', unpack=True)
dRl, dVl, dI, Ra = zeros((4, len(Rl)))
for i in range(len(Rl)):
	if Rl[i] < 200:
		dRl[i] = Rl[i] * 0.008 + 0.3
	elif Rl[i] < 2000:
		dRl[i] = Rl[i] * 0.008 + 1
	else:
		dRl[i] = Rl[i] * 0.008 + 10
	if I[i] < 2:
		dI[i] = I[i] * 0.005 + 0.001
		Ra[i] = 100
	elif I[i] < 20:
		dI[i] = I[i] * 0.005 + 0.01
		Ra[i] = 10
	else:
		dI[i] = I[i] * 0.012 + 0.1
		Ra[i] = 1
	dVl[i] = Vl[i] * 0.005 + 0.01

def fun_I(Rl, Rg, V0):
	return V0 / (Rg + Ra + Rl)

par, cov = curve_fit(fun_I, Rl, I, p0=(16.0, 5000.0), sigma=dI)

Rg, V0 = par
dRg, dV0 = sqrt(diag(cov))
c = cov[0, 1] / sqrt(cov[0, 0] * cov[1, 1])

chi2 = sum((I - fun_I(Rl, *par)) ** 2 / dI ** 2)

errorbar(Rl, I, yerr=dI, xerr=dRl, label='Data', fmt='.')
spaceRl = linspace(min(Rl), max(Rl), 1000)
Rl = sort(Rl)
Ra = sort(Ra)
plot(Rl, fun_I(Rl, *par), label='Fit')
title('Current versus load')
xscale('log')
xlabel('R [ohm]')
ylabel('I [A]')
legend(loc=0)
savefig('script4fig1.pdf')

print('numeric fit:')
print("Rg = %.3g +- %.3g" % (Rg, dRg))
print("V0 = %.3g +- %.3g" % (V0/1000, dV0/1000))
print("c = %.2f" % (c))
print("chi2/dof = %.3g/%d" % (chi2, len(Rl) - 2))
print('errors on current:')
print(dI)

show()

