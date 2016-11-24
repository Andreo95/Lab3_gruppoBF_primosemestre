#polaroid
import sys, os
import pylab as pl
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.gridspec as gridspec
sys.path.append(os.path.abspath('C:\\Users\\user\\Desktop\\laboratorio\\L2\\PYscripts\\'))
from lab import *
### un polaroid
Directory='C:\\Users\\user\\Desktop\\laboratorio\\L2\\19\\'
file='data1.txt'
FileName=(Directory+file)
I,the,dthe = pl.loadtxt(FileName).T
#converto I da microampere in ampere:
I=I*1e-6
dI = mme( I, 'ampere','digital')
#converto in radianti le misure angolari
the=the*(2*np.pi)/360
dthe=dthe*(2*np.pi)/360
#definisco funzione di fit
def Ip(Th,A,B,Th0):
    return A*(np.cos(Th-Th0))**2+B
#ottengo un errore complessivo per il fit, non considero il seno nella derivata di Ip rispetto a Th, tanto è minore di 1 quindi al massimo sovrastimo
dIn=np.sqrt(dI**2+1e-8*dthe**2)

#definiamo guess iniziali
start=np.array([1e-4,2*1e-6,1])#da vedere per bene!!!

#inizio a plottare prima che mi faccia errori
th=pl.linspace(min(the)-0.5,max(the)+0.5,1000)
gs=gridspec.GridSpec(5,1)
pl.subplot(gs[0:4,:])
pl.title('singolo filtro polaroid')
pl.xlabel('$\\theta$ [deg]')
pl.ylabel('I$_p$ [$\mu$A]')
pl.plot(the*360/(2*np.pi),I*1e6,'.',color='blue',label='dati')
pl.errorbar(the*360/(2*np.pi),I*1e6, dI*1e6, dthe*360/(2*np.pi), linestyle='',ecolor='blue')
# #control print:
#pl.plot(th*360/(2*np.pi),Ip(th,start[0],start[1],start[2])*1e6,color='red',label='checkprint')

#famoce ste azzo dde fitte
popt,pcov= curve_fit(Ip,the,I,start,absolute_sigma=False)
A,B,Th0=popt
dA,dB,dTh0=pl.sqrt(pcov.diagonal())

#estrapolo info statistiche
corrAB=pcov[1,0]/np.sqrt(pcov[0,0]*pcov[1,1])
corrATh0=pcov[2,0]/np.sqrt(pcov[0,0]*pcov[2,2])
corrBTh0=pcov[2,1]/np.sqrt(pcov[1,1]*pcov[2,2])
#
chisq=0.
h=0.
for h in range(len(the)):
    chisq+=((I[h]-Ip(the[h],A,B,Th0))/dIn[h])**2

#stampo risultati    
print(file+'\n'+'A=',A,'+-',dA,'\n','B=',B,'+-',dB,'\n','Th0=',Th0,'+-',dTh0,'\n','corrAB=',corrAB,'\n','corrATh0=',corrATh0,'\n','corrBTh0=',corrBTh0,'\n','X=',chisq,'dof=',len(the)-len(start),'\n',dI)

#graficamus
pl.plot(th*360/(2*np.pi),Ip(th,A,B,Th0)*1e6,color='orange',label='fit')
pl.legend(loc=1,fontsize=10)

#print dei residui
pl.subplot(gs[4:,:])
pl.xlabel('$\\theta$ [rad]')
pl.ylabel('norm.res.')
pl.rc('ytick',labelsize=8)
pl.plot(the*360/(2*np.pi),(I-Ip(the,A,B,Th0))/dIn,color='red')
pl.plot(the*360/(2*np.pi),(I-Ip(the,A,B,Th0))/dIn,'.',color='blue')
gs.update(hspace=1)

###due polaroid
file='data2.txt'
FileName=(Directory+file)
I,the,dthe = pl.loadtxt(FileName).T
#converto I da microampere in ampere:
I=I*1e-6
dI = mme( I, 'ampere','digital')
#converto in radianti le misure angolari
the=the*(2*np.pi)/360
dthe=dthe*(2*np.pi)/360
#definisco funzione di fit
def Ip2(Th,A,B,Th0,w):
    return A*(np.sin(w*(Th-Th0)))**2+B
#ottengo un errore complessivo per il fit, non considero il seno nella derivata di Ip rispetto a Th, tanto è minore di 1 quindi al massimo sovrastimo
dIn=np.sqrt(dI**2+2*1e-8*dthe**2)

#definiamo guess iniziali
start2=np.array([8*1e-7,1*1e-7,0,2])#da vedere per bene!!!

#inizio a plottare prima che mi faccia errori
th=pl.linspace(min(the)-0.5,max(the)+0.5,1000)
pl.figure()
gs=gridspec.GridSpec(5,1)
pl.subplot(gs[0:4,:])
pl.title('doppio filtro polaroid')
pl.xlabel('$\\theta$ [deg]')
pl.ylabel('I$_p$ [$\mu$A]')
pl.plot(the*360/(2*np.pi),I*1e6,'.',color='blue',label='dati')
pl.errorbar(the*360/(2*np.pi),I*1e6, dI*1e6, dthe*360/(2*np.pi), linestyle='',ecolor='blue')
# #control print:
#pl.plot(th*360/(2*np.pi),Ip2(th,start2[0],start2[1],start2[2])*1e6,color='red',label='checkprint')

#famoce ste azzo dde fitte
popt,pcov= curve_fit(Ip2,the,I,start2,absolute_sigma=False)
A,B,Th0,w=popt
dA,dB,dTh0,dw=pl.sqrt(pcov.diagonal())

#estrapolo info statistiche
corrAB=pcov[1,0]/np.sqrt(pcov[0,0]*pcov[1,1])
corrATh0=pcov[2,0]/np.sqrt(pcov[0,0]*pcov[2,2])
corrBTh0=pcov[2,1]/np.sqrt(pcov[1,1]*pcov[2,2])
#
chisq=0.
h=0.
for h in range(len(the)):
    chisq+=((I[h]-Ip2(the[h],A,B,Th0,w))/dIn[h])**2

#stampo risultati    
print(file+'\n'+'A=',A,'+-',dA,'\n'+'w=',w,'+-',dw,'\n','B=',B,'+-',dB,'\n','Th0=',Th0,'+-',dTh0,'\n','corrAB=',corrAB,'\n','corrATh0=',corrATh0,'\n','corrBTh0=',corrBTh0,'\n','X=',chisq,'dof=',len(the)-len(start),'\n',dI)

#graficamus
pl.plot(th*360/(2*np.pi),Ip2(th,A,B,Th0,w)*1e6,color='orange',label='fit')
pl.legend(loc=1,fontsize=10)

#print dei residui
pl.subplot(gs[4:,:])
pl.xlabel('$\\theta$ [rad]')
pl.ylabel('norm.res.')
pl.rc('ytick',labelsize=8)
pl.plot(the*360/(2*np.pi),(I-Ip2(the,A,B,Th0,w))/dIn,color='red')
pl.plot(the*360/(2*np.pi),(I-Ip2(the,A,B,Th0,w))/dIn,'.',color='blue')
gs.update(hspace=1)
