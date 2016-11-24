import sys, os
sys.path.append( os.path.join(os.path.realpath('..'), '00 - Risorse', 'Python') )

import numpy as np, pylab, lab, matplotlib.pyplot as plt, matplotlib as mpl
from normalTools import chiAvg

folder = os.path.realpath('.')

gs=mpl.gridspec.GridSpec(5,1)


####################### Parte 2.b

datafile = '{0}{1}_dati.txt'.format(2,'b')

rawdata = np.loadtxt( os.path.join(folder, 'Dati', datafile) ).T

R1 = 810 	# ohm
R2 = 1159	# ohm

Vin, Vout = rawdata

A = Vin
B = Vout

errA = lab.mme(A, unit="volt")
errB = lab.mme(B, unit="volt")

XX = A/B
errXX = XX * np.sqrt( (errA/A)**2 + (errB/B)**2 )

X, sigmaX = chiAvg(XX, errXX)

print("misurato: ", X , "+-" , sigmaX)
print("atteso: ", (R1 + R2)/R2, 'pm', (R1/R2)*(lab.mme(R1, unit='ohm')/R1 + lab.mme(R2, unit='ohm')/R2), '= misurato + ', (R1 + R2)/R2 - X )

### test 3


pars = lab.fit_affine_xyerr(B, A, errB, errA )
ratio = pars[0]
erratio = np.sqrt(pars[2])
print(ratio, 'pm', erratio)


### test 2

f = lambda x, m, q: m*x + q
df = lambda x, m, q: m

par, pcov = lab.fit_generic_xyerr(f, df, B, A, errB, errA, p0 = np.array([ratio, pars[1]]) )
ratio = par[0]
erratio = np.sqrt(np.diag(pcov)[0])

print(ratio, 'pm', erratio, '\n')


res = (A - f(B, ratio, par[1]))/errA
print('chi2 = ', sum(res**2))

plt.figure(1)
plt.clf()
plt.subplot(gs[:4,:])
plt.title('Partitore di tensione, bassa R')
plt.ylabel('Vin [$V$]')
plt.errorbar(B, A, errA, errB, fmt='none', ecolor='black')
plt.plot(B, f(B, ratio, par[1]), color='green')
plt.plot(B, f(B, (R1 + R2)/R2, 0), color='red')
plt.subplot(gs[4:,:])
plt.ylabel('Norm. res.')
plt.xlabel('Vout [$V$]')
plt.plot(B, (A - f(B, ratio, par[1]))/errA , 'ob')
plt.plot(B, np.zeros(len(B)), color='k')

#plt.savefig( os.path.join(folder, 'Figs-Tabs', 'fig_{}.png'.format(1)), bbox_inches='tight', dpi = 400) #format='ps', papertype='a4', orientation='landscape'

####################### parte 2.c/d

datafile = '{0}{1}_dati.txt'.format(2,'c')

rawdata = np.loadtxt( os.path.join(folder, 'Dati', datafile) ).T

Vin, Vout = rawdata

R2 = 1.03	*1e6	# Mohm
R1 = 1.518	*1e6	# Mohm

errR1, errR2 = lab.mme(np.array([R1,R2]), unit='ohm')

A = Vin
B = Vout

errA = lab.mme(A, unit="volt")
errB = lab.mme(B, unit="volt")

XX = A/B
errXX = XX * np.sqrt( (errA/A)**2 + (errB/B)**2 )

X, sigmaX = chiAvg(XX, errXX)

pars = lab.fit_affine_xyerr(B, A, errB, errA )
ratio = pars[0]
erratio = np.sqrt(pars[2])

print("misurato: ", X , "+-" , sigmaX)
print("atteso: ", (R1 + R2)/R2, '= misurato + ', (R1 + R2)/R2 - X )

Rmm = 1/((ratio - 1)/R1 - 1/R2)
errRmm = Rmm*np.sqrt( ( (ratio - 1)**2/R1**2 * (erratio**2/(ratio-1)**2 + (errR1/R1)**2 ) + errR2**2/R2**4) / ((ratio - 1)/R1 - 1/R2)**2)

print(Rmm, 'pm', errRmm)
print(lab.util_mm_esr2(B, metertype='analog', what='res') , '\n')


plt.figure(2)
plt.clf()
plt.title('Partitore di tensione, alta R')
plt.ylabel('Vin [$V$]')
plt.xlabel('Vout [$V$]')
plt.errorbar(B, A, errA, errB, fmt='none', ecolor='black')
plt.plot(B, f(B, ratio, par[1]), color='green')
plt.plot(B, f(B, (R1 + R2)/R2, 0), color='red')
plt.plot(B, f(B, (R1 + 1/(1/R2 + 1/Rmm))*(1/R2 + 1/Rmm), 0), ls='--', color='blue')

####################### parte 2e


datafile = '{0}{1}_dati.txt'.format(2,'e')

rawdata = np.loadtxt( os.path.join(folder, 'Dati', datafile) ).T

V = rawdata[0]
I2, I1 = rawdata[1:] * 1e-6

R1 = 560
R2 = 220
R3 = 98.9 * 1e3

Rz = np.array([R1,R2,R3])

errRz = lab.mme(Rz, unit='ohm')

A = I1
B = I2

errA = lab.mme(I1, unit='ampere', metertype='analog')
errB = lab.mme(I2, unit='ampere', metertype='analog')

f = lambda x, m : m*x
df = lambda x, m : m

pars, pcov= lab.fit_generic_xyerr(f, df, B, A, errB, errA )
ratio = pars[0]
erratio = np.sqrt(np.diag(pcov)[0])

print("misurato: ", ratio , "+-" , erratio)
print("atteso: ", R2/R1, '+-', (R2/R1 * np.sqrt((errRz[0]/R1)**2 + (errRz[1]/R1)**2) ))

linpars = lab.fit_affine_xyerr(B, A, errB, errA )


res = (I1 - f(I2, pars[0]))/errA
print('chi2 = ', sum(res**2))

res =  (I1 - (f(I2, linpars[0]) + linpars[1]) )/errA
print('chi2 (affine) = ', sum(res**2))

plt.figure(3)
plt.clf()
plt.subplot( mpl.gridspec.GridSpec(6,1)[:4] )
plt.title('Partitore di corrente')
plt.ylabel('I1 [$\mu A$]')
plt.errorbar(I2*1e6, I1*1e6, errA*1e6, errB*1e6, fmt='none', ecolor='black')
plt.plot(I2*1e6, f(I2, pars[0])*1e6, color='green')
plt.plot(I2*1e6, (f(I2, linpars[0]) + linpars[1])*1e6, color='blue')
plt.plot(I2*1e6, f(I2, R2/R1)*1e6, color='red')
plt.subplot( mpl.gridspec.GridSpec(6,1)[4:5] )
plt.ylabel('Norm. res. (linear)')
plt.plot(I2*1e6, (I1 - f(I2, pars[0]))/errA , 'og')
plt.plot(I2*1e6, np.zeros(len(I2)), color='k')
plt.subplot( mpl.gridspec.GridSpec(6,1)[5:] )
plt.ylabel('Norm. res. (affine)')
plt.xlabel('I2 [$\mu A$]')
plt.plot(I2*1e6, (I1 - (f(I2, linpars[0]) + linpars[1]) )/errA , 'ob') 
plt.plot(I2*1e6, np.zeros(len(I2)), color='k')

print( lab.util_mm_esr2(I1, metertype='analog', unit='ampere', what='res'), '\n' )


####################### parte 3b


datafile = '{0}{1}_dati.txt'.format(3,'')

rawdata = np.loadtxt( os.path.join(folder, 'Dati', datafile) ).T

R1 = 9910 	# ohm
R2 = 9900	# ohm

Vin, Vout = rawdata[:2]

A = Vin
B = Vout

errA , errB = rawdata[2:] * np.sqrt(2) + 0.03*np.array([A,B]) + 0.001 

XX = A/B
errXX = XX * np.sqrt( (errA/A)**2 + (errB/B)**2 )

X, sigmaX = chiAvg(XX, errXX)

print("misurato: ", X , "+-" , sigmaX)
print("atteso: ", (R1 + R2)/R2, 'pm', (R1/R2)*(lab.mme(R1, unit='ohm')/R1 + lab.mme(R2, unit='ohm')/R2), '= misurato + ', (R1 + R2)/R2 - X )

### test 3

f = lambda x, m, q: m*x + q
df = lambda x, m, q: m

pars = lab.fit_affine_xyerr(B, A, errB, errA )
ratio = pars[0]
erratio = np.sqrt(pars[2])
print(ratio, 'pm', erratio)


### test 2

f = lambda x, m, q: m*x + q
df = lambda x, m, q: m

par, pcov = lab.fit_generic_xyerr(f, df, B, A, errB, errA, p0 = np.array([ratio, pars[1]]) )
ratio = par[0]
erratio = np.sqrt(np.diag(pcov)[0])

print(ratio, 'pm', erratio)

res = (A - f(B, ratio, par[1]))/errA
print('chi2 = ', sum(res**2))

plt.figure(4)
plt.clf()
plt.subplot(gs[:4,:])
plt.title('Partitore di tensione, misura da oscilloscopio')
plt.ylabel('Vin [$V$]')
plt.errorbar(B, A, errA, errB, fmt='none', ecolor='black')
plt.plot(B, f(B, ratio, par[1]), color='green')
plt.plot(B, f(B, (R1 + R2)/R2, 0), color='red')
plt.subplot(gs[4:,:])
plt.ylabel('Norm. res.')
plt.xlabel('Vout [$V$]')
plt.plot(B,  res, 'ob')
plt.plot(B, np.zeros(len(B)), color='k')


plt.show()

