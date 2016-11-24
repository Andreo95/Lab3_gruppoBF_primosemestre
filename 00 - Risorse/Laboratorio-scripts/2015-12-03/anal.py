from lab import *
from pylab import *

f, Vin, Vout, phi = loadtxt('data.txt', unpack=True)

dpVout = array([(4.0 if fr < 100000 else 6.0) for fr in f])
dVout = dpVout / 100 * Vout
A = Vout / Vin
dA = A * sqrt((dVout/Vout)**2 + 0.04**2)
dB = 20 * log10(A)
ddB = 20 / log(10) * dA / A

fun_A = lambda f, fT: sqrt(1 / (1 + (f/fT)**2))
fun_phi = lambda f, fT: arctan(f/fT) / pi

par, cov = curve_fit_patched(fun_A, f, A, p0=(365.0), sigma=dA, absolute_sigma=True)

fT = par[0]
dfT = sqrt(cov[0,0])

par, cov = curve_fit_patched(fun_phi, f, phi, p0=(365.0))

fTp = par[0]
dfTp = sqrt(cov[0,0])

ff = logspace(floor(log10(min(f))), ceil(log10(max(f))), 256)

suptitle('Bode diagram of low-pass RC filter')

subplot(211)
ylabel('gain [dB]')
xscale('log')
ylim(-70, 10)
errorbar(f, dB, yerr=ddB, fmt=',', label='Data')
plot(ff, 20 * log10(fun_A(ff, fT)), label='Fit')
legend(loc=0)

subplot(212)
xlabel('frequency [Hz]')
ylabel('phase [pi rad]')
xscale('log')
plot(f, -phi, 'o', label='Data')
plot(ff, -fun_phi(ff, fT), label='Fit')
legend(loc=0)

savefig('plot.pdf')

print("\ncut frequency: %s Hz" % util_format(fT, dfT))
print("\ncut frequency (from phase): %s Hz\n" % util_format(fTp, dfTp))

for i in range(len(Vout)):
	print(util_format(Vout[i], dVout[i]))

chi2 = sum(((A - fun_A(f, fT)) / dA)**2)

print("chi square / dof: %f/%d" % (chi2, len(f) - 1))

show()

