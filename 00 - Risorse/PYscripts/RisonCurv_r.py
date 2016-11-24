#Curva di risonanza per RLC
import pylab as pl
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.gridspec as gridspec

Directory='D:\\Valerio\\L2\\16\\'
file='data.txt'
FileName=(Directory+file)
f, df, Vout, dVout, Vin, dVin = pl.loadtxt(FileName).T
#first manipulation
A=Vout/Vin
dA=((dVout/Vout)+(dVin/Vin))*A

#definiamo parametri
R=330
dR=8 #?#
rr=42.6
CC=2.2e-7
LL=0.5
#inizio a plottare prima che mi faccia errori
ff=pl.linspace(0,max(f)+min(f),1000)
gs=gridspec.GridSpec(5,1)
pl.subplot(gs[0:4,:])
pl.title('Curva di risonanza')
pl.xlabel('f [Hz]')
pl.ylabel('A(f)')
pl.plot(f,A,'.',color='blue',label='dati')
pl.errorbar(f, A, dA, df, linestyle='',ecolor='blue')

#via col fit
def Att(f,C,r,L):
    return 2*np.pi*R*C*f/np.sqrt((1-L*C*(2*np.pi*f)**2)**2+(2*np.pi*(R+r)*C*f)**2)
    
#definiamo guess iniziali
start=np.array([CC,rr,LL])#da vedere per bene!!!

# #control print:
# pl.plot(ff,Att(ff,start[0],start[1],start[2]),color='red')

#famoce ste azzo dde fitte
popt,pcov= curve_fit(Att,f,A,start,absolute_sigma=True)
C,r,L=popt
dC,dr,dL=pl.sqrt(pcov.diagonal())

#Ghettizzazione OUTLIERS
i=0
Ain=np.array([])
fin=np.array([])
dAin=np.array([])
for i in range(len(A)):
    if (np.abs(A[i]-Att(f[i],C,r,L))<=5*dA[i]):
        Ain=np.insert(Ain,len(Ain),A[i])
        fin=np.insert(fin,len(fin),f[i])
        dAin=np.insert(dAin,len(dAin),dA[i])
    else:
        pl.plot(f[i],A[i],'.',color='red')
pl.plot(3*max(ff)/2,3*max(A)/2,'.',color='red',label='outliers')
#remaking the fit
popt,pcov= curve_fit(Att,f,A,start,absolute_sigma=True)
C,r,L=popt
dC,dr,dL=pl.sqrt(pcov.diagonal())
#estrapolo info statistiche
corrCr=pcov[1,0]/np.sqrt(pcov[0,0]*pcov[1,1])
corrCL=pcov[2,0]/np.sqrt(pcov[0,0]*pcov[2,2])
corrrL=pcov[2,1]/np.sqrt(pcov[1,1]*pcov[2,2])
#
chisq=0.
h=0.
for h in range(len(fin)):
    chisq+=((Ain[h]-Att(fin[h],C,r,L))/dAin[h])**2
#fwhm
D=np.sqrt(3)*(R+r)/(L*2*np.pi)
dD=np.sqrt(3)*(dR+dr)/(L*2*np.pi)+dL*D/L
#fo
fo=1/(2*np.pi*np.sqrt(L*C))
dfo=(dL+dC)/(2*np.pi*(np.sqrt(L*C))**3)
#stampo risultati    
print(file+'\n'+'C=',C,'+-',dC,'\n','r=',r,'+-',dr,'\n','L=',L,'+-',dL,'\n','Dfwhm=',D,'+-',dD,'\n','fo=',fo,'+-',dfo,'\n','corrCr=',corrCr,'\n','corrCL=',corrCL,'\n','corrrL=',corrrL,'\n','X=',chisq,'dof=',len(fin)-len(start))

#graficamus
pl.plot(ff,Att(ff,C,r,L),color='orange',label='fit')
pl.legend(loc=1,fontsize=10)

#print dei residui
pl.subplot(gs[4:,:])
pl.xlabel('f [Hz]')
pl.ylabel('norm.res.')
pl.rc('ytick',labelsize=8)
#pl.plot(fin,(Ain-Att(fin,C,r,L))/dAin,color='red')
pl.plot(fin,(Ain-Att(fin,C,r,L))/dAin,'.',color='blue')
gs.update(hspace=1)