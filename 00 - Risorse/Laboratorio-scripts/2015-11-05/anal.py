from lab import *
from pylab import *

t = loadtxt('tensmult.txt', unpack=True)
dt = mme(t, 'volt')

v, dv = zeros((2, len(t)))

for i in range(len(t)):
	file = "data%d.txt" % (i + 1)
	data = loadtxt(file, unpack=True)
	avg = data.mean()
	sigma = sqrt(data.var(ddof=1) / len(data))
	v[i] = avg
	dv[i] = sigma
	figure(i + 1)
	clf()
	title("Tension measured with Arduino")
	subplot(211)
	xlabel('tension [AU]')
	ylabel('occurences')
	hist(data, bins=max(data)-min(data))
	subplot(212)
	xlabel('sample')
	ylabel('tension [AU]')
	errorbar(range(len(data)), data, yerr=1, fmt=',')
	savefig("graf%d.png" % (i + 1))
	xi = t[i] * 1000 / avg
	dxi = xi * sqrt((sigma / avg)**2 + (dt[i] / t[i])**2)
	vdata = data * xi / 1000.0
	dvdata = sqrt((vdata * dxi / xi)**2 + 1)
	vavg, vvar = fit_const_yerr(vdata, dvdata)
	vsigma = sqrt(vvar)
	print("data %d: avg = %g, sigma = %g, xi = %g, dxi = %g, V = %g, dV = %g" % (i + 1, avg, sigma, xi, dxi, vavg, vsigma))

m, q, vm, vq = fit_affine_xyerr(v, t, dv, dt, print_info=True)

dm = sqrt(vm)
dq = sqrt(vq)
chi2 = sum((t - m * v + q)**2 / (dv**2 + dt**2))
dof = len(v) - 2

print("m = %g +- %g\nq = %g +- %g\nchi2/dof = %g / %d = %g" % (m, dm, q, dq, chi2, dof, chi2/dof))

figure(len(t) + 1)
clf()
xlabel('tension [AU]')
ylabel('tension [V]')
errorbar(v, t, yerr=dt, xerr=dv, fmt='.')
plot(v, m * v + q)
savefig('fit.pdf')


