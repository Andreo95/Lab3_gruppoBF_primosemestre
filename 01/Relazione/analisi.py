import sys, os
sys.path.append( os.path.join(os.path.realpath('..'), '00 - Risorse', 'Python') )

import numpy as np, pylab, lab, matplotlib.pyplot as plt, matplotlib as mpl
from normalTools import chiAvg

folder = os.path.realpath('.')


## Luzio, usa i dannati spazi

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

for i in range(len(A)):
    print(Vin[i], "&", 1000*errA[i], "&", Vout[i], "&", 1000*errB[i], "\\"+"\\")



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


####################### parte 2.c/d


print("parte 2.c/d")

datafile = '{0}{1}_dati.txt'.format(2,'c')

rawdata = np.loadtxt( os.path.join(folder, 'Dati', datafile) ).T

Vin, Vout = rawdata

R2 = 1.03	*1e6	# Mohm
R1 = 1.518	*1e6	# Mohm
#R2=40*1e3
errR1, errR2 = lab.mme(np.array([R1,R2]), unit='ohm')

A = Vin
B = Vout
print(A)
print(B)
print(B/A)

errA = lab.mme(A, unit="volt")
errB = lab.mme(B, unit="volt")

XX = B/A
errXX = XX * np.sqrt( (errA/A)**2 + (errB/B)**2 )

X, sigmaX = chiAvg(XX, errXX)


XX=A/B
pars = lab.fit_affine_xyerr(B, A, errB, errA )
ratio = pars[0]
erratio = np.sqrt(pars[2])

print("misurato: ", X , "+-" , sigmaX)
print("atteso: ", R2/(R1+R2), '= misurato + ', R2/(R1+R2) - X )



####fin quì
Rmm = 1/((ratio - 1)/R1 - 1/R2)
errRmm = Rmm*np.sqrt( ( (ratio - 1)**2/R1**2 * (erratio**2/(ratio-1)**2 + (errR1/R1)**2 ) + errR2**2/R2**4) / ((ratio - 1)/R1 - 1/R2)**2)

print(Rmm, 'pm', errRmm)
print(lab.util_mm_esr2(B, metertype='analog', what='res') )

for i in range(len(A)):
    print(Vin[i], "&", 1000*errA[i], "&", Vout[i], "&", 1000*errB[i], "\\"+"\\")
    

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

print(errRz)

A = I1
B = I2

errVi=lab.mme(V, unit="volt", metertype="digital")
errA = lab.mme(I1, unit='ampere', metertype='analog')
errB = lab.mme(I2, unit='ampere', metertype='analog')

for i in range(len(A)):
    print(V[i], "&",errVi[i] , "&", I2[i]*1e6, "&", errB[i]*1e6, "&", I1[i]*1e6 , "&",errA[i]*1e6 ,"\\"+"\\")



f = lambda x, m : m*x
df = lambda x, m : m

pars, pcov= lab.fit_generic_xyerr(f, df, B, A, errB, errA )
ratio = pars[0]
erratio = np.sqrt(np.diag(pcov)[0])

print("misurato: ", ratio , "+-" , erratio)
print("atteso: ", R2/R1, '+-', (R2/R1 * np.sqrt((errRz[0]/R1)**2 + (errRz[1]/R1)**2) ))

print("cosa:", (I1+I2)/(V/R3))
plt.figure(1)
plt.title("I_1 vs I_2")
plt.xlabel("I_2 [uA]")
plt.ylabel("I_1 [uA]")
plt.errorbar(I2*1e6, I1*1e6, errA*1e6, errB*1e6, fmt='.', color='black')
plt.plot(I2*1e6, f(I2, pars[0])*1e6, color='green')
plt.plot(I2*1e6, f(I2, R2/R1)*1e6, color='red')

plt.figure(2)
plt.title("residui normalizzati I_2 vs I_1")
plt.xlabel("I_2 [uA]")
plt.ylabel("(I_2-I_2)at/I_2")
errRes=np.sqrt(((pars[0]/I1)*errB)**2+((I2*pars[0]/I1**2)*errA)**2)
plt.errorbar(I2*1e6, (I1- f(I2, pars[0]))/I1,errRes , errA*1e6)
 
print((R1+R2)*((V/R3)/(I1+I2)-np.ones(len(V))), "!!!!!!!!!!")


print( lab.util_mm_esr2(I1, metertype='analog', unit='ampere', what='res') )


####################### parte 3b


datafile = '{0}{1}_dati.txt'.format(3,'d')

rawdata = np.loadtxt( os.path.join(folder, 'Dati', datafile) ).T

R1 =0.985 # ohm
R2 =0.560	# ohm

Vin, Vout = rawdata[:2]





A = Vin
B = Vout

errA , errB = rawdata[2:] * np.sqrt(2) + 0.03*np.array([A,B]) + 0.001 

# for i in range(len(Vin)):
#     print(Vin[i], "&", np.sqrt((0.03*Vin[i])**2+errA[i]**2), "&", Vout[i],"&" ,  np.sqrt((0.03*Vout[i])**2+errB[i]**2), "\\"+"\\")

# print("...",1/((A/B-1)/R1-1/R2))
# 
# XX = A/B
# errXX = XX * np.sqrt( (errA/A)**2 + (errB/B)**2 )
# 
# X, sigmaX = chiAvg(XX, errXX)
# 
# print("misurato: ", X , "+-" , sigmaX)
# print("atteso: ", (R1 + R2)/R2, 'pm', (R1/R2)*(lab.mme(R1, unit='ohm')/R1 + lab.mme(R2, unit='ohm')/R2), '= misurato + ', (R1 + R2)/R2 - X )

### test 3

# f = lambda x, m, q: m*x + q
# df = lambda x, m, q: m
# 
# pars = lab.fit_affine_xyerr(B, A, errB, errA )
# ratio = pars[0]
# erratio = np.sqrt(pars[2])
# print(ratio, 'pm', erratio)
# 
# 
# ### test 2
# 
# f = lambda x, m, q: m*x + q
# df = lambda x, m, q: m
# 
# par, pcov = lab.fit_generic_xyerr(f, df, B, A, errB, errA, p0 = np.array([ratio, pars[1]]) )
# ratio = par[0]
# erratio = np.sqrt(np.diag(pcov)[0])
# 
# print(ratio, 'pm', erratio)


#############4


datafile = '{0}{1}_dati.txt'.format(4,'b')

rawdata = np.loadtxt( os.path.join(folder, 'Dati', datafile) ).T

print(rawdata)
freq, sp, dsp=rawdata
sp=sp*1e-3
dsp=dsp*1e-3
for i in range(len(freq)):
    print(freq[i], "&",1/(2*sp[i]), "&", dsp[i]/(2*sp[i]**2), "\\"+"\\")