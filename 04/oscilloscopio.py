import sys, os

sys.path.append( os.path.join(os.path.realpath('..'), '00 - Risorse', 'Python') )
folder = os.path.realpath('.')

import numpy as np, pylab, matplotlib.pyplot as plt, matplotlib as mpl, scipy.stats.distributions as dists
from lab import *
from ian import *
from PyTeX import *
import re
plt.close('all')
rawdata = None

####### setup

R_2=17.54e3
R_1=178.0e3
R_C=9.98e3
R_e=989 
R_list = [R_1, R_2, R_C, R_e]
R_errs = mme([R_1, R_2, R_C, R_e],'ohm')

cin=238e-9
cout=103.0e-9
ce= 100e-6
c_list = [cin,cout,ce]
c_errs = mme([cin,cout,ce], 'farad')

things = R_list + c_list
dthg = np.append(R_errs, c_errs)
qq = [x for pair in zip(things,dthg) for x in pair]
tt = np.array(qq)[:,np.newaxis]
# print(maketab(*tt,errors='all'))

## quiescenti 
V_e=1.10 
V_c=9.12 
V_cc=20.0 
V_b=1.68 
V_list = [V_e, V_c, V_cc, V_b]
V_errs = mme([V_e, V_c, V_cc, V_b], 'volt','oscil')

#### Parte 2

datafile = "datarelli1.txt"
rawdata=[]
t1=[]
ch1=[]
t2=[]
ch2=[]

try:
	f=os.path.join(folder, 'Dati', datafile)
	print(f, type(f))
	ff=open(f)
	for l in ff.readlines():
		try:
			raw=re.findall(r"[\+.\-.\w]+", l)
			t1.append(float(raw[0]))
			ch1.append(float(raw[1]))
			t2.append(float(raw[2]))
			ch2.append(float(raw[3]))
		except:
			print("la stringa "+l+"non è interpretabile")
except FileNotFoundError:
	print('Il file <'+f+">non sembra esistere")
T1=np.array(t1)
T2=np.array(t2)
CH1=np.array(ch1)
CH2=np.array(ch2)
dCH1=mme(CH1, unit="volt", metertype="oscil")
dCH2=mme(CH2, unit="volt", metertype="oscil")

pylab.figure(1)
pylab.title("Clipping estremo")
pylab.xlabel("Tempo [s]")
pylab.ylabel("Tensione [V]")
pylab.errorbar(T1, CH1, dCH1, label="CH1")
pylab.errorbar(T2, CH2, dCH2, label="CH2")

d=[np.amin(T1), np.amax(T1)]
AAnd=vectorize(lambda a, b: a and b)
print(np.amax(CH2[AAnd(T2<10e-6,T2>-10e-6)]))
x=np.mean(CH2[AAnd(T2<10e-6,T2>-10e-6)])
print(x)

pylab.plot(d, [x, x])

print(np.amax(CH2[AAnd(T2<-0.00008,T2>-0.00012)]))
x=np.mean(CH2[AAnd(T2<-0.00008,T2>-0.00012)])
print(x)

pylab.plot(d, [x, x])
pylab.show()

# Vin = (rawdata[0] - rawdata[1])
# dVin = rawdata[2]
# Vout = rawdata[3] - rawdata[4]
# dVout = rawdata[5]
# ammax=1.6
# used = Vin<ammax
# XX  = Vin
# YY  = Vout
# dXX = dVin/2*.7
# dYY = dVout/2*.7
# 
# linepars, linecov = fit_generic_xyerr(line, dline, XX[used], YY[used], dXX[used], dYY[used], p0=[10, 0])
# 
# resd = (YY[used] - line(XX[used], *linepars)) /dYY[used]
# DoFs = len(XX[used])-len(linepars)
# tellChi2(resd, DoFs)
# 
# plotLines = {
# 	line: dict(
# 		pars = linepars,
# 		linetype = dict(color=greencol),
# 		mask = used
# 	)
# }
# 
# 
# first = Graph(XX, YY, dXX, dYY, funFacts = plotLines)
# first.compute()
# first.whatplot |= {'data', line}
# #first.titolo="$V_{{in}}$ vs $V_{{out}} [V_{{cutoff}}={l})$".format("{in}","{out}","{cutoff}",ammax)
# first.titolo="$V_{{in}}$ vs $V_{out}$ [$V_{cutoff}$="+str(ammax)+"]"
# first.labelX="$V_{in}$ [V]"
# first.labelY="$V_{out}$ [V]"
# first.draw(resid=True)
# 
# print(linepars, linecov)
# print("cor={0}".format((linecov[0,0]*linecov[1, 1]/linecov[1, 0]**2)**-.5) )
# 
# # plt.savefig( os.path.join(folder, 'Figs-Tabs', 'fig_{}.png'.format(fignum)), bbox_inches='tight', dpi = 1000) # format='ps', papertype='a4', orientation='landscape'
# 
# 
# #### Parte 3
# 
# 
# 
# ###### ending
# #plt.draw_all()
# plt.show()
# Graph.countfigs.send(0)