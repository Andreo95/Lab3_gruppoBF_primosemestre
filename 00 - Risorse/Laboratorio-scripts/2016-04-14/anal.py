from pylab import *
from lab import *

f, vo, dvo, vi, dvi = loadtxt('data.txt', unpack=True)
df = array([1.0] * len(f))
a = vo / vi
da = sqrt((dvo/vo)**2 + (dvi/vi)**2) * a

# R = 343
r = 40.3
A = lambda f, f0, R, C: R*C / sqrt((R+r)**2 * C**2 + (1/f - f/f0**2)**2 / (2*pi)**2)

p0 = (616.0, 343, 0.1e-6)
par, cov = fit_generic_xyerr2(A, f, a, df, da, p0)
sigma = sqrt(diag(cov))
chi2 = sum((A(f,*par)-a)**2/(da**2))
print("chi2 = %g, ndof = %d" % (chi2, len(f) - len(par)))
# names = ['f0', 'R', 'r', 'C']
names = ['f0', 'R', 'C']
for i in range(len(par)):
	print("%s = %s" % (names[i], util_format(par[i], sigma[i])))
ncov = cov.copy()
for i in range(len(par)):
	for j in range(len(par)):
		ncov[i,j] = cov[i,j]/(sigma[i]*sigma[j])
print(ncov)
ff = linspace(min(f), max(f), 1000)
fa = A(ff, *par)
for i in range(len(fa) - 1):
	if (fa[i+1]>max(fa)/2) and (fa[i]<=max(fa)/2):
		fmin = ff[i]
	if (fa[i+1]<max(fa)/2) and (fa[i]>=max(fa)/2):
		fmax = ff[i]
fwhm = fmax - fmin
print("fwhm = %d" % (fwhm))
title('RLC circuit gain vs frequency')
plot([fmin] * 2, [min(a), max(a)], 'r')
plot([fmax] * 2, [min(a), max(a)], 'r')
plot(ff, fa, label='Fit')
errorbar(f, a, yerr=da, xerr=df, fmt=',', capsize=0, label='Data')
xlabel('Frequency [Hz]')
ylabel('Gain')
legend(loc=0)
savefig("FiguraBella.pdf")
