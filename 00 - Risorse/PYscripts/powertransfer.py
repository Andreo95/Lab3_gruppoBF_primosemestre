#Trasferimento di Potenza (trasformatori)
import sys, os
import pylab as pl
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.gridspec as gridspec
sys.path.append(os.path.abspath('D:\\Valerio\\L2\\PYscripts\\'))
from lab import *

Directory='D:\\Valerio\\L2\\18\\'
file='data.txt'
FileName=(Directory+file)
R2,V2,dV2 = pl.loadtxt(FileName).T
dR2 = mme( R2, 'ohm')
#first manipulation
P2j=0.5*V2**2/R2
dP2j=P2j*(2*dV2/V2+dR2/R2)
# P2j*=1e6
# dP2j*=1e6
#definiamo parametri
R1=33
dR1=0.3 #?#
V1=0.4
a=20
#inizio a plottare prima che mi faccia errori
r2=pl.logspace(0,7,1000)
gs=gridspec.GridSpec(5,1)
pl.subplot(gs[0:4,:])
pl.xscale('log')
pl.title('Trasferimento di potenza al circuito secondario')
pl.xlabel('R$_2$ [$\\Omega$]')
pl.ylabel('P$_2$ [$\mu$W]')
pl.plot(R2,P2j*1e6,'.',color='blue',label='dati')
pl.errorbar(R2,P2j*1e6, dP2j*1e6, dR2, linestyle='',ecolor='blue')

#via col fit
def P2(r2,a,V1):
    return ((V1**2)/2)*(r2/a**2)/(R1+r2/a**2)**2
    
#definiamo guess iniziali
start=np.array([a,V1])#da vedere per bene!!!

#control print:
# pl.plot(r2,P2(r2,start[0],start[1]),color='red')

#famoce ste azzo dde fitte
popt,pcov= curve_fit(P2,R2,P2j,start,absolute_sigma=False)
a,V1=popt
da,dV1=pl.sqrt(pcov.diagonal())

#Ghettizzazione OUTLIERS
i=0
P2in=np.array([])
R2in=np.array([])
dP2in=np.array([])
dR2in=np.array([])
scarto=np.array([])
for i in range(len(P2j)):
    if (np.abs(P2j[i]-P2(R2[i],a,V1))<=5*dP2j[i]):
        P2in=np.insert(P2in,len(P2in),P2j[i])
        R2in=np.insert(R2in,len(R2in),R2[i])
        dR2in=np.insert(dR2in,len(dR2in),dR2[i])
        dP2in=np.insert(dP2in,len(dP2in),dP2j[i])
    else:
        pl.plot(R2[i],P2j[i],'.',color='red')
        scarto=np.insert(scarto,len(scarto),R2[i])
        scarto=np.insert(scarto,len(scarto),P2j[i])
if(len(scarto)>=1):
    pl.plot(scarto[0],scarto[1],'.',color='red',label='outliers')
#remaking the fit
popt,pcov= curve_fit(P2,R2in,P2in,start,absolute_sigma=False)
a,V1=popt
da,dV1=pl.sqrt(pcov.diagonal())
#estrapolo info statistiche
corrV1a=pcov[1,0]/np.sqrt(pcov[0,0]*pcov[1,1])
#
chisq=0.
h=0.
for h in range(len(R2in)):
    chisq+=((P2in[h]-P2(R2in[h],a,V1))/dP2in[h])**2

#stampo risultati    
print(file+'\n'+'a=',a,'+-',da,'\n','V1=',V1,'+-',dV1,'\n''corrV1a=',corrV1a,'\n','X=',chisq,'dof=',len(R2in)-len(start))

#graficamus
pl.plot(r2,1e6*P2(r2,a,V1),color='orange',label='fit')
pl.legend(loc=1,fontsize=10)

#print dei residui
pl.subplot(gs[4:,:])
pl.xscale('log')
pl.xlabel('R$_2$ [$\\Omega$]')
pl.ylabel('norm.res.')
pl.rc('ytick',labelsize=8)
pl.plot(R2in,(P2in-P2(R2in,a,V1))/dP2in,color='red')
pl.plot(R2in,(P2in-P2(R2in,a,V1))/dP2in,'.',color='blue')
gs.update(hspace=1)