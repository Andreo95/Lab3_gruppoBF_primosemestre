import pylab
import numpy as np
import matplotlib.gridspec as gridspec

#simulazione filtro passa-basso con onde quadre

#definisco le caratteristiche del filtro
fa=50
fb=25e3
#definisco il termine k-esimo della serie trasformata in base alle funzioni di trasferimento di integratore e derivatore
def f(x,k,fa,fb,ff):
    return (1/np.sqrt(1+(fb/(ff*(2+4*k)))**2))*(1/np.sqrt(1+(ff*(2+4*k)/fa)**2))*4*np.sin(ff*(2+4*k)*np.pi*x+np.arctan(-ff*(2+4*k)*2/fa) +np.arctan(fb/(ff*(2+4*k))))/((2*k+1)*np.pi)

#mi preparo a ciclare
n=10
i=0
u=0
b=0
m=1000  #numero di punti usati
ff=50   #prima frequenza usata; si andrà avanti raddoppiando tale valore 9 volte
#plotto l'output per diversi valori della frequenza, usando i primi 200 termini della serie
pylab.figure()
gs = gridspec.GridSpec(5, 2)
for i in range(n):
    #mi creo l' array 'x' calibrando l'intervallo alla semiampiezza di 2T:
    a=2/(ff*2**i)
    x=np.array([])
    for b in range(m-1):
        x=np.insert(x,len(x),-a+2*a*b/m)
    somma=np.zeros(len(x))
    #costruisco la somma:
    for u in range(201):
        somma+=f(x,u,fa,fb,ff*2**i)
    #plotto in modo che ogni frequenza vada in un subplot diverso: alla fine sono riuscito a non metterci if! (anche se segue un giro "a ferro di cavallo")
    pylab.subplot(gs[min(i,9-i),int(max(0,i-4)/(i-4-0.0001))])
    pylab.rc('ytick',labelsize=6)
    pylab.rc('xtick',labelsize=10)
    pylab.plot(x,somma)
    pylab.text((7/12)*a,max(somma)*5/8,'f='+str((ff*2**i)/1000)+str('kHz'), size=10)

#Traffico di attenuazioni
#estrapolazione delle attenuazioni:
F=np.logspace(0,6,100) #frequenze campione che userò
#calcolo i relativi output:
M=100
p=0
V=np.array([])
for i in range(len(F)):
    #mi creo l' array 'x' calibrando l'intervallo alla semiampiezza di 2T:
    a=1/(F[i])
    xx=np.array([])
    for p in range(M):
        xx=np.insert(xx,len(xx),-a+2*a*p/M)
    somm=np.zeros(len(xx))
    #costruisco la somma:
    for u in range(20):#nota: la scelta del numero di termini di ogni sviluppo e del numero di punti M è determinante!!
        somm+=f(xx,u,fa,fb,F[i])
    V=np.insert(V,len(V),np.amax(somm))
pylab.title('attenuazione integratore+derivatore: attenuazione vs frequenza')
pylab.xscale('log')
pylab.yscale('log')
pylab.xlabel('f [Hz]')
pylab.ylabel('A [u.a.]')
pylab.plot(F,V)
pylab.plot(F,1/np.sqrt(1+(F/fa)**2),'--',color='grey')
pylab.plot(F,(fa/fb)*F**0,'--',color='grey')
pylab.plot(F,(1/np.sqrt(1+(F/fa)**2))/np.sqrt(1+(fb/F)**2),':k',color='grey')