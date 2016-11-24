import numpy
import pylab
import matplotlib.gridspec as gridspec

#mi creo un array al posto di usare linspace che mi dava problemi incomprensibili
m=5000
x=numpy.array([])
a=2
u=0
for u in range(m):
    x=numpy.insert(x,len(x),-a+2*a*u/m)
#definisco il termine k-esimo della serie
def f(x,k):
    return 4*numpy.sin((2+4*k)*numpy.pi*x)/((2*k+1)*numpy.pi)
somma=numpy.ndarray([])
n=5
i=0
u=0
b=0
#plotto le somme parziali e i residui in subplot separati
gs = gridspec.GridSpec(15, 1)
for i in range(n):
    for u in range(b,10*i+1):
        somma= somma+f(x,u)
        b=10*i+1
    pylab.subplot(gs[(3*i):(3*i+2),0])
    pylab.rc('ytick',labelsize=6)
    pylab.plot(x,somma, label='n='+str(10*i+1))
    pylab.legend()
    pylab.subplot(gs[(3*i+2):(3*i+3),0])
    pylab.rc('ytick',labelsize=6)
    pylab.plot(x,1-numpy.abs(somma),color='lightblue')
    gs.update(hspace=0.07)

##triangolare, stesso procedimento
    
# def f2(x,h):
#     return (8*numpy.cos((2+4*h)*numpy.pi*x))/(((2*h+1)*numpy.pi)**2)
# summ=numpy.ndarray([])
# summst=numpy.ndarray([])
# r=5
# k=0
# j=0
# w=0
# c=0
# for k in range(10**2):
#     summst= summst+f2(x,k)
# pylab.figure()
# gs = gridspec.GridSpec(15, 1)
# for j in range(r):
#     for w in range(c,10*j+1):
#         summ= summ+f2(x,w)
#         c=10*j+1
#     pylab.subplot(gs[(3*j):(3*j+2),0])
#     pylab.rc('ytick',labelsize=6)
#     pylab.plot(x,summ, label='n='+str(10*j+1))
#     pylab.legend()
#     pylab.subplot(gs[(3*j+2):(3*j+3),0])
#     pylab.rc('ytick',labelsize=6)
#     pylab.plot(x,numpy.abs(summst-summ),color='lightblue')
#     gs.update(hspace=0.07)