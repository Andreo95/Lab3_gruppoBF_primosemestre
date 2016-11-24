import numpy as np
import pylab as pl
import matplotlib.gridspec as gridspec
#definisco i parametri caratteristici
mo=1e5
s=mo
m=1
u=1e7
ro=np.sqrt(m/s)
#definisco le funzioni cruciali per il calcolo numerico della posizione del centro di massa
def psi(t):
    return (2*np.pi/(ro*np.sqrt(np.pi*s))+(8/5)*(u/m)*np.sqrt(np.pi*s))*(np.sqrt(mo)-np.sqrt(mo-m*t))
def ay(t):
    return ((m*m/(ro*s)+m*u/5)/(mo-m*t))*np.cos(psi(t))
def ax(t):
    return -((m*m/(ro*s)+m*u/5)/(mo-m*t))*np.sin(psi(t))
    
#pronti a ciclare; nota: t<mo/m
dt=1e-15
#n=int(round(mo/(m*dt)))
n=15000
i=0
Vx=np.zeros(n)
Vy=np.zeros(n)
X=np.zeros(n)
Y=np.zeros(n)
Vx[0]=0
Vy[0]=0
for i in range(n-1):
    Vx[i+1]=Vx[i]+dt*ax(dt*i)
    Vy[i+1]=Vy[i]+dt*ay(dt*i)
    X[i+1]=X[i]+dt*Vx[i]
    Y[i+1]=Y[i]+dt*Vy[i]
#plot
#tt=np.linspace(0,5,1000000)
#gs=gridspec.GridSpec(5,1)
#pl.subplot(gs[0:4,:])
pl.plot(X,Y)
#pl.subplot(gs[4:,:])
#pl.plot(tt,psi(tt),'--',color='black')
#pl.plot(tt,ax(tt),color='blue')
#pl.plot(tt,ay(tt),color='orange')