import pylab
import numpy as np
import matplotlib.gridspec as gridspec

#simulazione filtro passa-basso con onde quadre

#definisco le caratteristiche del filtro
R=330
C=1e-7
ff=50
#definisco il termine k-esimo della serie trasformata in base alla funzione di trasferimento
def f(x,k,R,C,ff):
    return (1/np.sqrt(1+(R*C*ff*(2+4*k)*np.pi)**2))*4*np.sin(ff*(2+4*k)*np.pi*x+np.arctan(-R*C*ff*(2+4*k)*np.pi))/((2*k+1)*np.pi)
#mi preparo a ciclare
n=10
i=0
u=0
b=0
m=1000
#plotto l'output per diversi valori della frequenza, usando i primi 200 termini della serie
pylab.figure()
gs = gridspec.GridSpec(5, 2)
for i in range(n):
    #mi creo l' array 'x' calibrando l'intervallo alla semiampiezza di 2T:
    a=2/(ff*2**i)
    x=np.array([])
    for b in range(m):
        x=np.insert(x,len(x),-a+2*a*b/m)
    somma=np.zeros(len(x))
    #costruisco la somma:
    for u in range(201):
        somma+=f(x,u,R,C,ff*2**i)
    #plotto in modo che ogni frequenza vada in un subplot diverso: alla fine sono riuscito a non metterci if! 
    pylab.subplot(gs[min(i,9-i),int(max(0,i-4)/(i-4-0.0001))])
    pylab.rc('ytick',labelsize=6)
    pylab.rc('xtick',labelsize=10)
    pylab.plot(x,somma)
    pylab.text((7/12)*a,max(somma)*5/8,'f='+str((ff*2**i)/1000)+str('kHz'), size=10)
    