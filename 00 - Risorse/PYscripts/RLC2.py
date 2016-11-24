import numpy as np
import pylab as pl
from scipy.optimize import curve_fit
import matplotlib.gridspec as gridspec

Directory='D:\\Valerio\\L2\\15\\data\\'
file='ferrlamn.txt'
FileName=(Directory+file)
t,dt,Vc,dVc = pl.loadtxt(FileName).T
#definiamo parametri
r=42.6
C=2.2e-7
B=270
#inizio a plottare prima che mi faccia errori
tt=pl.linspace(0,max(t)+500,300)
gs=gridspec.GridSpec(5,1)
pl.subplot(gs[0:4,:])
pl.title('Vc(t) ddp per circuito RLC; nucleo:'+file)     ###aggiungere nucleo al titolo
pl.xlabel('t [us]')
pl.ylabel('Vc [u.a.]')
pl.plot(t,Vc,'.',color='blue',label='dati')
pl.errorbar(t, Vc, dVc, dt, linestyle='',ecolor='blue')
#definiamo funzione
def V(t,tau,w,phi,C,A):
    return A*np.exp(-t/tau)*np.cos(w*t+phi)+C
def V1(t,tau,C):
    return (B)*np.exp(-t/tau)+C
#definiamo guess iniziali
start=np.array([14000,(2*np.pi/3000),1.7,500,350])#da vedere per bene!!!
# #control print:
# pl.plot(t,Vc,'.',color='blue')
# pl.plot(tt,V(tt,start[0],start[1],start[2],start[3],start[4]),color='red')
#famoce ste azzo dde fitte
popt,pcov= curve_fit(V,t,Vc,start,absolute_sigma=True)
tau,w,phi,C,A=popt
dtau,dw,dphi,dC,dA=pl.sqrt(pcov.diagonal())
#estrapolo info statistiche
#corr=pcov[1,0]/numpy.sqrt(pcov[0,0]*pcov[1,1])
#Ghettizzazione OUTLIERS
i=0
inn=np.array([])
tinn=np.array([])
dinn=np.array([])
for i in range(len(Vc)):
    if (np.abs(Vc[i]-V(t[i],tau,w,phi,C,A))<=5*dVc[i]):
        inn=np.insert(inn,len(inn),Vc[i])
        tinn=np.insert(tinn,len(tinn),t[i])
        dinn=np.insert(dinn,len(dinn),dVc[i])
    else:
        pl.plot(t[i],Vc[i],'.',color='red')
pl.plot(max(tt)+500,max(Vc)+100,'.',color='red',label='outliers')
#remaking the fit
popt,pcov= curve_fit(V,tinn,inn,start,absolute_sigma=True)
tau,w,phi,C,A=popt
dtau,dw,dphi,dC,dA=pl.sqrt(pcov.diagonal())
#
chisq=0.
h=0.
for h in range(len(tinn)):
    chisq+=((inn[h]-V(tinn[h],tau,w,phi,C,A))/dVc[h])**2
#stampo risultati    
print(file+'\n'+'A=',A,'+-',dA,'\n','tau=',tau,'+-',dtau,'\n','w=',w,'+-',dw,'\n','phi=',phi,'+-',dphi,'\n','X=',chisq,'dof=',len(tinn)-5)
#graficamus
pl.plot(tt,V(tt,tau,w,phi,C,A),color='orange',label='fit')
pl.legend(loc=3,fontsize=10)
#stampo in text le info del fit
param=np.array([tau,w*1000,chisq])
err=np.array([dtau,dw*1000])
err=np.around(err,3)#careful!!
param=np.around(param,2)#arrotondare a seconda degli errori
pl.text(max(tt)*4/5,min(Vc)*3/5, '$\\tau=$'+str(param[0])+'+-'+str(err[0])+',\n$\omega=$'+str(param[1])+'+-'+str(err[1])+',\n$\chi^2_{rid}=$'+str(param[2])+'/'+str(len(tinn)-5),verticalalignment='bottom', horizontalalignment='left',color='black', fontsize=10)
#print dei residui
pl.subplot(gs[4:,:])
pl.xlabel('t [us]')
pl.ylabel('norm.res.')
pl.rc('ytick',labelsize=8)
pl.plot(tinn,(inn-V(tinn,tau,w,phi,C,A))/dinn,'.',color='blue')
pl.plot(tinn,(inn-V(tinn,tau,w,phi,C,A))/dinn,color='red')
gs.update(hspace=1)