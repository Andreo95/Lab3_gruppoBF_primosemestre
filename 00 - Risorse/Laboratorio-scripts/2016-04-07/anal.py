from pylab import *
from lab import *
from scipy.optimize import curve_fit

filelist = loadtxt('Data/cippa.txt', dtype='|S32')
allah = 1e+10
cuts = array([allah, allah, allah, 3, 7, 7, 6, 8, 8, 7, 5, 6, 3]) * 1000
C = 0.1e-6
dC = C / 10
ifilename = 0
for filename in filelist:
	t, v = loadtxt('Data/%s' % filename, unpack=True)
	nt = []
	nv = []
	for i in range(len(t)):
		if t[i] <= cuts[ifilename]:
			nt.append(t[i])
			nv.append(v[i])
	t = array(nt)
	v = array(nv)
	dv = array([1.0] * len(v))
	f_fit = lambda t, omega, tau, bias, ampl, phi: bias - ampl * exp(-t/tau) * sin(omega*t + phi)
	p0 = (6.3e-3, 2000, 500, 500, 0)
	try:
		par, cov = curve_fit_patched(f_fit, t, v, sigma=dv, p0=p0, absolute_sigma=False)
		omega, tau, _, _, _ = par
		domega, dtau, _, _, _ = sqrt(diag(cov))
		L = 1/C * 1 / ((omega*1e6)**2 + 1/(tau*1e-6)**2)
		dL = (dC/C + (2*omega*domega - 2*dtau/tau)) * L
		normres = (v - f_fit(t, *par)) / dv
		chi2 = sum(normres ** 2)
		ndof = len(v) - len(par)
		print("\n%s" % filename)
		print("omega = %s" % xe(omega, domega))
		print("tau = %s" % xe(tau, dtau))
		print("chi2/ndof = %g / %d" % (chi2, ndof))
		print("L = %s" % xe(L, dL))
		clf()
		suptitle(filename)
		subplot(211)
		errorbar(t, v, yerr=dv, fmt=',', capsize=0, label='Data')
		ft = linspace(min(t), max(t), 1024)
		plot(ft, f_fit(ft, *par), label='Fit')
		ylabel('voltage [arb. un.]')
		interv = (max(t) - min(t)) * 0.1
		xlim(min(t) - interv, max(t) + interv)
		subplot(212)
		plot(t, normres, '.')
		xlabel('time [us]')
		ylabel('norm. res.')
		xlim(min(t) - interv, max(t) + interv)
	except:
		clf()
		errorbar(t, v, yerr=dv, fmt=',', capsize=0, label='Data')
		ft = linspace(min(t), max(t), 1024)
		plot(ft, f_fit(ft, *p0), label='Initial params')
		print("\nFAILED %s" % filename)
	savefig("plot-%s.pdf" % filename)
	ifilename += 1

