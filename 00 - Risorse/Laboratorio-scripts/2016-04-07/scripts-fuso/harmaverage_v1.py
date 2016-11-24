# Questo script serve per interfacciarsi con Arduino nell'esperienza
# di carica e scarica del condensatore e gestire l'acquisizione multipla
# e la produzione di un file contenente quattro colonne (con medie e deviazioni standard):
# t, Deltat, V, DeltaV

# L'interfacciamento avviene attraverso:
# 1. scrittura di un carattere (byte) che esprime l'intervallo di campionamento
# 2. lettura dei dati disponibili su porta seriale

import serial # libreria per gestione porta seriale (USB)
import time # libreria per temporizzazione
import numpy

nmis=5 # numero acquisizioni da effettuare << DA MODIFICARE COME SI VUOLE

Directory='../dati_arduino/'   # nome directory dati
                                                    # << DA CAMBIARE SECONDO NECESSITA'
FileName=Directory+'smorzaverage.txt'  # parte comune nome file << DA CAMBIARE SECONDO NECESSITA'


# crea gli arrays necessari
nomt=["" for x in range (0,nmis)]
nomv=["" for x in range (0,nmis)]

for i in range (0,nmis):
    nomt[i]='tim'+str(i)
    nomv[i]='vol'+str(i)
    nomt[i]=numpy.zeros(256)
    nomv[i]=numpy.zeros(256)

print ("start")

# avvia il loop per le acqusizioni
for imis in range (0,nmis):
    ard=serial.Serial('/dev/ttyACM0',19200)  # apre porta seriale (occhio alla sintassi, dipende
                                # dal sistema operativo!)
    time.sleep(2)   # aspetta due secondi per evitare casini
    
    ard.write(b'5') # scrive il carattere per l'intervallo di campionamento
                    # in unita' di 10 us << DA CAMBIARE A SECONDA DEI GUSTI
                    # l'istruzione b indica che e' un byte (carattere ASCII)
    
    time.sleep(2) # aspetta due secondi per evitare casini

    for i in range (0,256):
        data = ard.readline().decode() # legge il dato e lo decodifica
        if data:
            nomt[imis][i]=int(data[0:data.find(' ')])
            nomv[imis][i]=int(data[data.find(' '):len(data)])
	    
    ard.close() # chiude la comunicazione seriale con Arduino
    print('Part ',imis+1,' of ',nmis, ' done')

# processa i dati per ricavare media e standard deviation (sia sui tempi che sulle tensioni)
runningtime=numpy.zeros(nmis)
runningvolt=numpy.zeros(nmis)
t=numpy.zeros(256)
dt=numpy.zeros(256)
V=numpy.zeros(256)
dV=numpy.zeros(256)

for i in range (0,256):
    for imis in range (0,nmis):
        runningtime[imis]=nomt[imis][i]
        runningvolt[imis]=nomv[imis][i]
    t[i]=numpy.average(runningtime)
    dt[i]=numpy.std(runningtime)
    V[i]=numpy.average(runningvolt)
    dV[i]=numpy.std(runningvolt)

for i in range (0,256):
        if (dV[i]==0):
                dV[i]=numpy.average(dV)
        if (dt[i]==0): 
                dt[i]=numpy.average(dt)

outputFile = open(FileName, "w" ) # apre file dati per scrittura

for i in range (0,256):
    outputFile.write(str(t[i]))
    outputFile.write('  ')
    outputFile.write(str(dt[i]))
    outputFile.write('  ')
    outputFile.write(str(V[i]))
    outputFile.write('  ')
    outputFile.write(str(dV[i]))
    outputFile.write('\r')
        
outputFile.close() # chiude il file dei dati



print('end') # scrive sulla console che ha finito
