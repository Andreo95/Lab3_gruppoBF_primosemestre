import numpy as np
import pylab as pl
from scipy.optimize import curve_fit
import matplotlib.gridspec as gridspec

Directory='D:\\Valerio\\L2\\14\\'
FileName=(Directory+'data1.txt')
t,Vc = pl.loadtxt(FileName).T
#definiamo parametri
dVc=1
r=42.6
C=2.2e-7
B=270
#inizio a plottare prima che mi faccia errori
tt=pl.linspace(0,18000,300)
gs=gridspec.GridSpec(5,1)
pl.subplot(gs[0:4,:])
pl.title('Vc(t) ddp per circuito RLC')
pl.xlabel('t [us]')
pl.ylabel('Vc [u.a.]')
pl.plot(t,Vc,'.',color='blue',label='dati')
pl.errorbar(t, Vc, dVc, linestyle='',ecolor='blue')
#definiamo funzione
def V(t,tau,w,phi,C):
    return (B)*np.exp(-t/tau)*np.cos(w*t+phi)+C
def V1(t,tau,C):
    return (B)*np.exp(-t/tau)+C
#definiamo guess iniziali
startingval=np.array([14000,(2*np.pi/1300),1.7,500])#da vedere per bene!!!
#famoce ste azzo dde fitte
popt,pcov= curve_fit(V,t,Vc,startingval,dVc)
tau,w,phi,C=popt
dtau,dw,dphi,dC=pl.sqrt(pcov.diagonal())
#estrapolo info statistiche
#corr=pcov[1,0]/numpy.sqrt(pcov[0,0]*pcov[1,1])
chisq=0.
h=0.
for h in range(len(t)):
    chisq+=((Vc[h]-V(t[h],tau,w,phi,C))/dVc)**2
#stampo risultati    
print('tau=',tau,'+-',dtau,'\n','T=',2*np.pi/w,'+-',dw,'\n','phi=',phi,'+-',dphi,'\n','X=',chisq,'dof=',len(t)-4)
#estrapole stime cristiche di R e L da tau e w:
L=tau*r/2*1e-6
R=2*L/(np.sqrt(1/(L*C-w**2)))
# L=tau*R/2
print('R=',R,'\n','L=',L)
Q0=(B/C)*np.cos(phi)#stima rozza
#grafichiamo
pl.plot(tt,V(tt,tau,w,phi,C),color='orange',label='fit')
pl.plot(tt,V1(tt,tau,C),'--',color='red')
pl.legend(loc=3,fontsize=10)
pl.subplot(gs[4:,:])
pl.xlabel('t [us]')
pl.ylabel('norm.res.')
pl.rc('ytick',labelsize=8)
pl.plot(t,(Vc-V(t,tau,w,phi,C))/dVc,'.',color='blue')
pl.plot(t,(Vc-V(t,tau,w,phi,C))/dVc,color='red')
gs.update(hspace=1)