import sys, os

sys.path.append( os.path.join(os.path.realpath('..'), '00 - Risorse', 'Python') )
folder = os.path.realpath('.')

import numpy as np, pylab, matplotlib.pyplot as plt, matplotlib as mpl, scipy.stats.distributions as dists
from lab import *
from ian import *
from PyTeX import *

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

## calcolo punto lavoro 

def workpt(V_BE_wish):
	beth = 178
	qq = 1/R_1 + 1/R_2
	lll = (1 + (1 + beth)*R_e*qq)
	mmm = V_cc/R_1 - V_BE_wish*qq
	I_b_att = mmm/lll
	I_c_att = beth*I_b_att
	I_e_att = (beth + 1)*I_b_att
	V_c_att = V_cc - R_C*I_c_att
	V_e_att = I_e_att * R_e
	V_b_att = (V_cc/R_1 - I_b_att)/qq
	return [I_b_att, I_c_att, I_e_att, V_c_att, V_e_att, V_b_att]

#### Parte 2

datafile = 'dati_{0}{1}.txt'.format(2,'')

try:
	f=os.path.join(folder, 'Dati', datafile)
	rawdata = np.loadtxt( os.path.join(folder, 'Dati', datafile) ).T
except FileNotFoundError:
	print('Il file <'+f+">non sembra esistere")

Vin = (rawdata[0] - rawdata[1])
dVin = rawdata[2]
Vout = rawdata[3] - rawdata[4]
dVout = rawdata[5]
ammax=1.6
used = Vin<ammax
XX  = Vin
YY  = Vout
dXX = dVin/2*.7
dYY = dVout/2*.7

linepars, linecov = fit_generic_xyerr(line, dline, XX[used], YY[used], dXX[used], dYY[used], p0=[10, 0])

resd = (YY[used] - line(XX[used], *linepars)) /dYY[used]
DoFs = len(XX[used])-len(linepars)
tellChi2(resd, DoFs)

plotLines = {
	line: dict(
		pars = linepars,
		linetype = dict(color=greencol),
		mask = used
	)
}


first = Graph(XX, YY, dXX, dYY, funFacts = plotLines)
first.compute()
first.whatplot |= {'data', line}
#first.titolo="$V_{{in}}$ vs $V_{{out}} [V_{{cutoff}}={0})$".format(ammax)
first.titolo="$V_{{in}}$ vs $V_{out}$ [$V_{cutoff}$="+str(ammax)+"]"
first.labelX="$V_{in}$ [V]"
first.labelY="$V_{out}$ [V]"
first.draw(resid=True)

print(linepars, linecov)
print("cor={0}".format((linecov[0,0]*linecov[1, 1]/linecov[1, 0]**2)**-.5) )

# plt.savefig( os.path.join(folder, 'Figs-Tabs', 'fig_{}.png'.format(fignum)), bbox_inches='tight', dpi = 1000) # format='ps', papertype='a4', orientation='landscape'


#### Parte 3

datafile = 'dati_{0}{1}.txt'.format(3,'')

try:
	rawdata = np.loadtxt( os.path.join(folder, 'Dati', datafile) ).T
except FileNotFoundError:
	print('Hai sbagliato file')


freq = rawdata[0]
dfreq = freq/1000
Vout = rawdata[1] - rawdata[2]
dVout = rawdata[3]

#guadagno atteso										################################à
f=np.linspace(10,1e11,100000)
Z_in= 1/(2*np.pi*1j*f*cin)
Z_out=1/(2*np.pi*1j*f*cout)
hfe = 1000
hie=10e3
#s=1j*f/(2*np.pi)
R_B=(R_1*R_2)/(R_1+R_2)
G=abs( ( (hfe-1)*(R_C+Z_out) +R_C )/( (((R_B + Z_in)/R_B)*( (R_e*(1+hfe)) +hie)) +Z_in))
#primo= (R_C*hfe)/(hie +R_e*(1+hfe))
#secondo= (1+hie +R_e*(1+hfe))/(cin*(hie + R_e*(1+hfe))*R_B )
#terzo= ((R_C*hfe)+hie+R_B)/(cout*(hie+R_B))
#G= abs(primo*s*(1/(s+secondo))*(1/(s+terzo)))


plt.figure(8)
plt.plot(np.log10(f),G)

#for i in range (len(rawdata[0])):
#	print(freq[i],'&',dfreq[i],'&',Vout[i],'&',dVout[i],'\ \ ')

Vin = .352 + .440	#Volt
dVin = mme(Vin, 'volt', 'oscil')

Av = 20*np.log10(Vout/Vin)
dAv = 20*np.log10(np.e)*np.sqrt(dVin**2/Vin**2, dVout**2/Vout**2)

used = (freq>500) & (freq<11000)
XX  = freq
YY  = Av
dXX = dfreq
dYY = dAv

flatpar, flatvar = fit_const_yerr(YY[used], dYY[used])
print('centro',flatpar,flatvar)
resd = (YY[used] - const(XX[used], flatpar)) /dYY[used]
DoFs = len(XX[used])-1
tellChi2(resd, DoFs)


plotLines.update({
	logline: dict(
		pars = linepars,
		linetype = dict(color=greencol)
	),
	const: dict(
		pars = (flatpar,),
		linetype = dict(color='green'),
		mask = used
	)
})

third = Graph(XX, YY, dXX, dYY)
third.titolo=" Guadagno in tensione al variare della frequenza "
third.labelX="log(f) [Hz]"
third.labelY="$A_v$ [dB]"

###retta a destra
f = lambda x,a,b: a*x +b
df = lambda x,a,b: a

XX  = np.log10(freq)
YY  = Av
dXX = np.ones(len(freq))*(1/1000)
dYY = dAv

used =(freq> 350000)
pars, pcov = fit_generic_xyerr(f, df, XX[used], YY[used], dXX[used], dYY[used], p0=None)
print('dx',pars,pcov)

resd = (YY[used] - f(XX[used], *pars)) /dYY[used]
DoFs = len(XX[used])-len(pars)
print('ChiSquare = {0} ({1} DoF, p = {2})\n'.format( np.sum(resd**2), DoFs, dists.chi2.sf(np.sum(resd**2), DoFs) ))

def rightline(*args):
	return logline(*args)

plotLines.update({
	rightline: dict(
		pars = pars,
		linetype = dict(color='blue')
	)
})


rightcross = FindCross(pars, pcov, flatpar, flatvar)

print( 'incrocio a destra: ', FindCross(pars, pcov, flatpar, flatvar))

## retta a sinistra
f = lambda x,a,b: a*x +b
df = lambda x,a,b: a

XX  = np.log10(freq)
YY  = Av
dXX = np.ones(len(freq))*(1/1000)
dYY = dAv

used =(freq< 34)
pars, pcov = fit_generic_xyerr(f, df, XX[used], YY[used], dXX[used], dYY[used], p0=None)
print('sx',pars,pcov)

resd = (YY[used] - f(XX[used], *pars)) /dYY[used]
DoFs = len(XX[used])-len(pars)
print('ChiSquare = {0} ({1} DoF, p = {2})\n'.format( np.sum(resd**2), DoFs, dists.chi2.sf(np.sum(resd**2), DoFs) ))

def leftline(*args):
	return logline(*args)

plotLines.update({
	leftline: dict(
		pars = pars,
		linetype = dict(color=greencol)
	)
})


third.funPars = plotLines
third.typeX = 'log'
third.compute()
third.whatplot |= {'data', const, rightline, leftline}
third.draw(resid=False)

leftcross = FindCross(pars, pcov, flatpar, flatvar)

print( 'incrocio a sinistra: ', FindCross(pars, pcov, flatpar, flatvar))


# plt.savefig( os.path.join(folder, 'Figs-Tabs', 'fig_{}.png'.format(fignum)), bbox_inches='tight', dpi = 1000) # format='ps', papertype='a4', orientation='landscape'


#### Parte 4

datafile = 'dati_{0}{1}.txt'.format(4,'')

try:
	rawdata = np.loadtxt( os.path.join(folder, 'Dati', datafile) ).T
except FileNotFoundError:
	print('Hai sbagliato file')

Vin = (rawdata[0] - rawdata[1])
dVin = rawdata[2]
Vout = rawdata[3] - rawdata[4]
dVout = rawdata[5]

used = Vin<1000
XX  = Vin
YY  = Vout
dXX = dVin/2*.7
dYY = dVout/2*.7

linepars, linecov = fit_generic_xyerr(line, dline, XX[used], YY[used], dXX[used], dYY[used], p0=[10, 0])

resd = (YY[used] - line(XX[used], *linepars)) /dYY[used]
DoFs = len(XX[used])-len(linepars)
tellChi2(resd, DoFs)

plotLines = {
	line: dict(
		pars = linepars,
		linetype = dict(color=greencol),
		mask = used
	)
}

second = Graph(XX, YY, dXX, dYY, funFacts = plotLines)
second.compute()
second.whatplot |= {'data', line}
second.draw(resid=True)

###### ending
#plt.draw_all()
plt.show()
Graph.countfigs.send(0)
