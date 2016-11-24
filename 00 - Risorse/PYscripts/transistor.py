from pylab import *

ibs = array([22.3, 20.4, 17.3, 13.7, 8.0])
dibs = ibs * 0.005 + 0.1 + 0.1

title('NPN transistor: collector current vs. collector-emitter tension')
ylabel('collector current [mA]')
xlabel('collector-emitter tension [V]')

ics = zeros(5)
dics = zeros(5)

for i in range(5):
	ib = ibs[i]
	dib = dibs[i]
	xi = 4.76e-3
	dxi = 0.033e-3
	data = loadtxt('data%02d.txt' % (5 - i), unpack=True) * xi
	v1, vce = data
	dv1 = v1 * sqrt((xi/v1)**2 + (dxi/xi)**2)
	dvce = vce * sqrt((xi/vce)**2 + (dxi/xi)**2)
	rc = 989.0
	drc = rc * 0.008 + 1.0
	ic = (v1 - vce) / rc
	dic = ic * sqrt(((dv1 + dvce)/(v1 - vce))**2 + (drc/rc)**2)
	dicp = sqrt(2.0) * xi / rc
	errorbar(vce, ic*1000.0, yerr=dicp*1000.0, xerr=xi, fmt=',', capsize=0.0, label='Ib = %.1f $\mu$A' % ib)
	legend(loc=0)
	for j in range(len(vce)):
		if vce[j] <= .7 and vce[j + 1] >= .7:
			ics[i] = ic[j] * 1000
			dics[i] = dic[j] * 1000
			break

for i in range(5):
	ic = ics[i]
	dic = dics[i]
	ib = ibs[i]
	dib = dibs[i]
	bf = ic / ib * 1000
	dbf = bf * sqrt((dic/ic)**2 + (dib/ib)**2)
	print("%.1f +- %.1f\t%.2f +- %.2f\t%.0f +- %.0f" % (ib, dib, ic, dic, bf, dbf))

savefig('plot.pdf')

############ part 2

clear()

from scipy.optimize import curve_fit

title('NPN transistor: collector current vs. collector-emitter tension')
ylabel('collector current [mA]')
xlabel('collector-emitter tension [V]')

ib = 8.0
dib = ib * 0.005 + 0.1 + 0.1

xi = 4.76e-3
dxi = 0.033e-3

data = loadtxt('data01.txt', unpack=True)[:,75:] * xi

v1, vce = data
dv1 = v1 * sqrt((xi/v1)**2 + (dxi/xi)**2)

dvce = vce * sqrt((xi/vce)**2 + (dxi/xi)**2)
dvcep = xi

rc = 989.0
drc = rc * 0.008 + 1.0

ic = (v1 - vce) / rc
dic = ic * sqrt(((dv1 + dvce)/(v1 - vce))**2 + (drc/rc)**2)
dicp = sqrt(2.0) * xi / rc

errorbar(vce, ic*1000.0, yerr=dic*1000.0, xerr=dvce, fmt=',', capsize=0.0, label='Ib = %.1f $\mu$A' % ib)

def f(vce, a, b):
	return a * vce + b
par, cov = curve_fit(f, vce, ic, sigma=dic)
a, b = par
da, db = sqrt(diag(cov))
c = cov[0,1]/(da*db)
print('corr = %.2f' % c)

ib *= 1e-6
dib *= 1e-6

bf = b / ib
dbf = bf * sqrt((db/b)**2 + (dib/ib)**2)

vea = b / a
dvea = vea * sqrt((db/b)**2 + (da/a)**2 + 2*(da/a)*(db/b)*c)

print('bf = %.0f +- %.0f' % (bf, dbf))
print('vea = %.0f +- %.0f' % (vea, dvea))

fvce = linspace(min(vce), max(vce), 1000)
plot(fvce, f(fvce, *par) * 1000, label='Fit')

legend(loc=0)
savefig('plot2.pdf')

