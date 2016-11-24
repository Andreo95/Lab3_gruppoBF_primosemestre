import pylab as py
import numpy as np
from scipy.optimize import curve_fit
import random

Directory='D:\\Valerio\\L2\\5\\'
FileName=(Directory+'data_C.txt')

t,V = py.loadtxt(FileName,unpack='True')
#simulazione
# t=np.array([])
# err=np.array([])
# V=np.array([])
# i=0.
# for i in range(250):
#     t=np.insert(t,len(t),i*0.2)
#     err=np.insert(err,len(err),random.gauss(0.,4))
# for u in range(250):
#     V=np.insert(V,len(V),1023.2*(1-np.exp(-t[u]/43.2))-5.01+err[u])
dV = np.ones(len(V))

def exp(t,tau,Vo,off):
    return Vo*(1-np.exp(-t/tau))+off
    
#nel caso di scarica:    
#return Vo*np.exp(-t/TS))+Off

init=py.array([6000,1000,-4])
popt,pcov= curve_fit(exp,t,V,init,dV)
tau,Vo,off    =popt
dtau,dVo,doff  =py.sqrt(pcov.diagonal())

npcov = py.zeros((3,3))
for i in range(3):
    for j in range(3):
        npcov[i,j]=pcov[i,j]/py.sqrt(pcov[i,i]*pcov[j,j])
print(npcov)

print('tau=',tau,'+-',dtau,'; Vo=',Vo,'+-',dVo,'; off=',off,'+-',doff,'; MCov=',pcov)

tt=np.linspace(0,30000,800)    
py.subplot(2,1,1)
py.errorbar(t,V,dV,fmt='.',color='blue')
py.plot(tt,exp(tt,tau,Vo,off),color='red')
py.rc('font',size=14)
py.xlabel('Tempo [ms]',size='16')
py.ylabel('valore  [digit]',size='16')
py.title('carica del condensatore',size='18')
py.minorticks_on()
py.legend(['data','fit'],prop={'size':14})

py.subplot(2,1,2)
py.plot(t,(V-exp(t,tau,Vo,off)),color='red')
py.rc('font',size=14)
py.ylabel('residui norm. [digit]',size='16')
py.xlabel('Tempo [ms]',size='16')
py.minorticks_on()
py.show()

chi_quadro=0.
for h in range(len(t)):
    chi_quadro+=((V[h]-exp(t[h],tau,Vo,off))/dV[h])**2

print('Chi=',chi_quadro, 'dof=',len(t))
