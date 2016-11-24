import numpy as np
import pylab as pl


a=10
b=28
c=8/3
def Vct(Vi,Vc,a):
    return a*(Vi-Vc)
def Vit(Vc,Vi,Il,b):
    return b*Vc-Vi-Vc*Il
def Ilt(Vc,Vi,Il,c):
    return Vc*Vi-c*Il
dt=1/1000
n=108000
Vc=np.zeros(n)
Vi=np.zeros(n)
Il=np.zeros(n)
i=0
Vi[0]=2
Vc[0]=4.9
Il[0]=10.3
for i in range((n-1)):
    Vc[i+1]=Vc[i]+(1/2)*dt*Vct(Vi[i],Vc[i],a)
    Vi[i+1]=Vi[i]+(1/2)*dt*Vit(Vc[i],Vi[i],Il[i],b)
    Il[i+1]=Il[i]+(1/2)*dt*Ilt(Vc[i],Vi[i],Il[i],c)
    Vc[i+1]=Vc[i+1]+(1/2)*dt*Vct(Vi[i+1],Vc[i+1],a)
    Vi[i+1]=Vi[i+1]+(1/2)*dt*Vit(Vc[i+1],Vi[i+1],Il[i+1],b)
    Il[i+1]=Il[i+1]+(1/2)*dt*Ilt(Vc[i+1],Vi[i+1],Il[i+1],c)
    
pl.plot(Vc,Vi,'-',)
# pl.plot(0,0,label=a)
# pl.plot(0,0,label=b)
# pl.plot(0,0,label=c)
# pl.plot(0,0,label=R)
# pl.legend()
