import numpy as np
import pylab as pl
import matplotlib.gridspec as gridspec

wc=1e-3
n=10
X=1e-5
x=np.linspace(5,10000,10000)
z=np.linspace(0,1000,10000)
col=['red','orange','yellow','green','cyan','blue','violet','brown','black','gray']
def A(x):
    return wc*np.sqrt(x/2)
def B(x):
    return wc*np.sqrt((np.sqrt(1+x**2)-1)/2)
def C(x):
    return wc*np.sqrt((np.sqrt(1+x**2)+1)/2)

def a(z):
    return np.exp(-A(X)*z)
def b(z):
    return np.exp(-B(X)*z)
def c1(z):
    return np.cos(A(X)*z)
def c2(z):
    return np.cos(C(X)*z)

pl.figure(1)
pl.title('lunghezza di pelle in funzione di wp/wy')
gs=gridspec.GridSpec(5,1)
pl.subplot(gs[0:4,:])
pl.plot(x,1/A(x),color='orange',label='w<<y & yw<<wp^1/2')
pl.plot(x,1/B(x),color='blue',label='w<<y')
pl.legend()
pl.subplot(gs[4:,:])
pl.plot(x,1/A(x)-1/B(x),color='red')

pl.figure(2)
pl.title('campo nel conduttore in funzione di wp/wy')
# gss=gridspec.GridSpec(5,1)

i=0
for i in range(n):
    X=X*10
#     pl.subplot(gss[0:4,:])
    pl.plot(z,a(z)*c1(z),color=col[i],label='wp/yw='+str(X))
    pl.plot(z,b(z)*c2(z),color=col[i])
#     pl.subplot(gss[4:,:])
#     pl.plot(z,a(z)*c1(z)-b(z)*c2(z),color=col[i])
pl.legend()