#Diodo p-n
import numpy as np
import pylab as pl
import matplotlib.gridspec as gsp
from scipy.optimize import curve_fit

Directory='C:\\Users\\user\\Desktop\\laboratorio\\L2\\9\\'
FileName=(Directory+'dio1.txt')

#NOTA: se avrò errori è sicuramente perchè ho messo i nomi interi per le librerie al posto che np,pl,gsp.. ti auguro di non perderci troppi secoli

#Dati e manipolazione
V1d,V2d = pl.loadtxt(FileName)[3:].T
Vref=4.97
dVref=0.03485      #comprende incertezza di lettura +-1 digit e incertezza di calibrazione del tester
Rd=3510
 #valori probabili! SOSTITUIRE!
dRd=30
Ecal=Vref/1023
dEcal=dVref/1023
V1=V1d*Ecal
V2=V2d*Ecal
dV1=dEcal*V1d
dV2=dEcal*V2d
nVt=52*1e-3     #guess sui parametri
Io=5*1e-9       #
I=(V1-V2)/Rd
dI=np.sqrt((dV1+dV2)**2/(Rd**2)+(dRd**2)*(V1-V2)**2/(Rd**4)+(dV2**2)*((I+Io)/nVt)**2)#ci propago anche l'errore su V2; 5% di probabilità che ho sbagliato
#largo a Shockley!
def ISh(V2,nVt,Io):             #ISh sta per I di Shockley
    return Io*(np.exp(V2/nVt)-1)
#fittamus papam:
startingval=np.array([nVt,Io])
popt,pcov= curve_fit(ISh,V2,I,startingval,dI)
nVt,Io    =popt
dnVt,dIo  =pl.sqrt(pcov.diagonal())
#estrapolo info statistiche
corr=pcov[1,0]/np.sqrt(pcov[0,0]*pcov[1,1])
chisq=0.
for h in range(len(V2)):
    chisq+=((I[h]-ISh(V2[h],nVt,Io))/dI[h])**2
#stampo risultati    
print('nVt=',nVt,'+-',dnVt,'\n','Io=',Io,'+-',dIo,'\n','corr=',corr,'\n','X=',chisq,'dof=',len(V2)-2)
print(Ecal,dEcal)
#realizzo grafico
gs=gsp.GridSpec(5,1)
v=np.linspace(0,0.7,200)
pl.suptitle('Diodo p-n: curva caratteristica e residui normalizzati')
pl.subplot(gs[0:4,:])
pl.ylabel('I [mA]')
pl.plot(V2,I,'.',color='blue',label='data')
pl.errorbar(V2,I,dI, linestyle='',ecolor='blue')
pl.plot(v,ISh(v,nVt,Io),color='red',label='fit')
pl.legend(loc=2)
pl.subplot(gs[4:,:])
pl.xlabel('$\Delta$V [V]')
pl.ylabel('norm.res.')
pl.plot(V2,(I-ISh(V2,nVt,Io))/dI,'-',color='blue')
pl.plot(V2,(I-ISh(V2,nVt,Io))/dI,',',color='red')
gs.update(hspace=1)