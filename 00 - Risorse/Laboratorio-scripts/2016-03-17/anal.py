from lab import *
from pylab import *

for ifile in [1, 22, 47]:
	C = float("0.%d" % ifile)
	t, v = loadtxt("smorz%d.txt" % ifile, unpack=True)[:,1:]
	
	dt = array([diff(t).std(ddof=1)] * len(v))
	dv = array([1] * len(v))
	
	def f_v(t, v0, tau, omega, phi, bias):
		return v0 * exp(-t/tau) * cos(omega*t + phi) + bias
	
	def df_v(t, v0, tau, omega, phi, bias):
		return v0 * (-1.0/tau * exp(-t/tau) * cos(omega*t + phi) - omega * exp(-t/tau) * sin(omega*t + phi))
	
	p0 = (400, 17800, 4.5e-3, 1.5, 425)
	
	par, cov = fit_generic_xyerr(f_v, df_v, t, v, dt, dv, p0)
	sigma = sqrt(diag(cov))
	
	print('capacitor: %d' % ifile)
	
	res = (f_v(t, *par) - v) / sqrt(dv**2 + (df_v(t, *par) * dt)**2)
	chi2 = sum(res ** 2)
	ndof = len(t) - len(par)
	
	print("chi2/ndof = %.0f/%d" % (chi2, ndof))
	omega = par[2]
	domega = sigma[2]
	T = 2*pi / omega
	dT = T * domega / omega
	print("T = %s" % util_format(T, dT))
	tau = par[1]
	dtau = sigma[1]
	print("tau = %s" % util_format(tau, dtau))
	L = 1.0 / (omega**2 * C) * 1e-6
	dL = L * sqrt((2*domega/omega)**2 + 0.1**2)
	print("L = %s" % util_format(L, dL))
	r = 2 * L / tau * 1e6
	dr = r * sqrt((dL/L)**2 + (dtau/tau)**2)
	print("r = %s" % util_format(r, dr))
	bias = par[-1]
	dbias = sigma[-1]
	print("bias = %s" % util_format(bias, dbias))
	print('')
	
	clf()
	suptitle('Capacitor voltage in RLC circuit')
	
	subplot(211)
	
	ylabel('voltage [digit]')
	bd = (max(t) - min(t)) / 25
	xlim(min(t) - bd, max(t) + bd)
	
	errorbar(t, v, yerr=dv, xerr=dt, fmt='.', capsize=0, label='data')
	
	ft = linspace(min(t), max(t), 2000)
	
	plot(ft, f_v(ft, *par), label='fit')
	# plot(ft, f_v(ft, *p0), label='in. par.')
	
	legend(loc=0)
	
	subplot(212)
	
	plot(t, res, '.')
	
	ylabel('norm. res.')
	xlabel('time [us]')
	xlim(min(t) - bd, max(t) + bd)
	
	savefig("plot%d.pdf" % ifile)


