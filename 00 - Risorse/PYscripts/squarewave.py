import numpy
import pylab

#mi creo un array al posto di usare linspace così posso trafficare senza problemi
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
#plotto le somme parziali e i residui ad ogni potenza di 10
for i in range(n):
    for u in range(b,10**i):
        b=10**i
        somma= somma+f(x,u)
    pylab.subplot(2,1,1)
    pylab.title('onda quadra')
    pylab.xlabel('t')
    pylab.ylabel('g(t)')
    pylab.plot(x,somma, label='n=10^'+str(i))
    pylab.legend()
    pylab.subplot(2,1,2)
    pylab.title('residui')
    pylab.xlabel('t')
    pylab.ylabel('1-|g(t)|')
    pylab.plot(x,1-abs(somma),label='n=10^'+str(i))
    pylab.legend()
     
## onda triangolare; nota: i conti sono molti, sconsiglio di eseguire questa parte se non strettamente necessario
    
# def f2(x,h):
#     return (8*numpy.cos((2+4*h)*numpy.pi*x))/(((2*h+1)*numpy.pi)**2)
# summ=numpy.ndarray([])
# summst=numpy.ndarray([])
# r=5
# k=0
# j=0
# w=0
# c=0
# for k in range(10**5):
#     summst= summst+f2(x,k)
# pylab.figure()
# for j in range(r):
#     for w in range(c,10**j):
#         c=10**j
#         summ= summ+f2(x,w)
#     pylab.subplot(2,1,1)
#     pylab.title('onda triangolare')
#     pylab.xlabel('t')
#     pylab.ylabel('g(t)')
#     pylab.plot(x,summ, label='n=10^'+str(j))
#     pylab.legend()
#     pylab.subplot(2,1,2)
#     pylab.title('residui')
#     pylab.xlabel('t')
#     pylab.ylabel('|standard-g(t)|')
#     pylab.plot(x,abs(summst-summ),label='n=10^'+str(j))
#     pylab.legend()