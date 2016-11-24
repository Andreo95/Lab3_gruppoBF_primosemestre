import numpy
import pylab
from scipy.optimize import curve_fit
import matplotlib.gridspec as gridspec

#caratteristiche dell'esperimento
R=3280
C=1e-7
#estrazione misure da file .txt
Directory='C:\\Users\\user\\Desktop\\laboratorio\\L2\\8\\'
FileName=(Directory+'atten.txt')
f, df, Vi, dVi, Vo, dVo, lag, dlag = pylab.loadtxt(FileName).T #nota: trascuro df: verifico che una stima tramite le cifre variabili dà errore inferiore all'1%
#manipolazione dati: attenuazione
Att=Vo/Vi
dAtt=(dVo*Vi+dVi*Vo)/Vi**2
#preparazione fit per passa basso
def A(f,a,ft):
    return a/numpy.sqrt(1+(f/ft)**2)
    
startingval=numpy.array([1,1/(2*numpy.pi*R*C)])
#fit numerico
popt,pcov= curve_fit(A,f,Att,startingval,dAtt)
a,ft    =popt
da,dft  =pylab.sqrt(pcov.diagonal())
#estrapolo info statistiche
corr=pcov[1,0]/numpy.sqrt(pcov[0,0]*pcov[1,1])
chisq=0.
for h in range(len(f)):
    chisq+=((Att[h]-A(f[h],a,ft))/dAtt[h])**2
#stampo risultati    
print('a=',a,'+-',da,'\n','ft=',ft,'+-',dft,'\n','corr=',corr,'\n','X=',chisq,'dof=',len(f)-2)
#realizzo grafico
gs=gridspec.GridSpec(5,1)
ff=numpy.logspace(1,6,500)
pylab.subplot(gs[0:4,:])
pylab.title('attenuazione passa basso: dati e fit numerico')
pylab.xlabel('f [Hz]')
pylab.ylabel('A')
pylab.xscale('log')
pylab.yscale('log')
pylab.plot(ff,A(ff,a,ft),color='red',label='fit')
pylab.plot(f,Att,'.',color='black',label='data')
pylab.errorbar(f, Att, dAtt, linestyle='',ecolor='black')
pylab.legend()
pylab.subplot(gs[4:,:])
pylab.xlabel('f [Hz]')
pylab.ylabel('norm.res.')
pylab.xscale('log')
pylab.rc('ytick',labelsize=8)
pylab.plot(f,(Att-A(f,a,ft))/dAtt,'.',color='black')
pylab.plot(f,(Att-A(f,a,ft))/dAtt,color='red')
gs.update(hspace=1)

#diagramma di bode
#mi preparo a fittare come se non ci fosse un domani
def dB(f,k):
    return k-20*numpy.log10(f)
start=numpy.array([-20*numpy.log10(2*numpy.pi*R*C)])
AttdB=20*numpy.log10(Att)
dAttdB=20*dAtt/(Att*numpy.log(10))
#seleziono i soli punti corrispondenti a frequenze f>1000
f1=numpy.array([])
AttdB1=numpy.array([])
dAttdB1=numpy.array([])
for i in range(len(f)):
    if (f[i] >=1000):
        f1=numpy.insert(f1,len(f1),f[i])
        AttdB1=numpy.insert(AttdB1,len(AttdB1),AttdB[i])
        dAttdB1=numpy.insert(dAttdB1,len(dAttdB1),dAttdB[i])
#fit numerico
popt,pcov= curve_fit(dB,f1,AttdB1,start,dAttdB1)
k    =popt
dk  =pylab.sqrt(pcov.diagonal())
#estrapolo info varie
chisq1=0.
for h in range(len(f1)):
    chisq1+=((AttdB1[h]-dB(f1[h],k))/dAttdB1[h])**2
fT=10**(k/20)
dfT=(dk*k/20)*(10**((k/20)-1))
print('fT=',fT,'+-',dfT,'\n','X=',chisq1,'dof=',len(f1)-1)
#grafico
pylab.figure()
pylab.title('diagramma di Bode con filtro passa basso')
pylab.xlabel('f [Hz]')
pylab.grid()
pylab.ylabel('A [dB]')
pylab.xscale('log')
pylab.plot(ff,dB(ff,k),color='red',label='fit')
pylab.plot(f,AttdB,'.',color='black',label='data')
pylab.errorbar(f, AttdB, dAttdB, linestyle='',ecolor='black')
pylab.plot(ff,0*ff,'--',color='green')
pylab.legend()