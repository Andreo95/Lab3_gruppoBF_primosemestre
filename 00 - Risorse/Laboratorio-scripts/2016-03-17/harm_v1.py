# Questo script serve per interfacciarsi con Arduino nell'esperienza
# di acquisizione delle oscillazioni smorzate di un circuito rLC

# L'interfacciamento avviene attraverso:
# 1. scrittura di un carattere (byte) che esprime l'intervallo di campionamento
# 2. lettura dei dati disponibili su porta seriale

import serial # libreria per gestione porta seriale (USB)
import time # libreria per temporizzazione
import numpy

Directory='./'   # nome directory dati
                                                    # << DA CAMBIARE SECONDO NECESSITA'
FileName=Directory+'smorz47.txt'  # parte comune nome file << DA CAMBIARE SECONDO NECESSITA'


ard=serial.Serial('/dev/ttyACM0',9600)  # apre porta seriale (occhio alla sintassi, dipende
                                # dal sistema operativo!)
time.sleep(2)   # aspetta due secondi per evitare casini

ard.write(b'5') # scrive il carattere per l'intervallo di campionamento
                # in unita' di 10 us << DA CAMBIARE A SECONDA DEI GUSTI
                # l'istruzione b indica che e' un byte (carattere ASCII)

time.sleep(2) # aspetta due secondi per evitare casini

outputFile = open(FileName, "w" ) # apre file dati carica per scrittura

print("")
print ("start")
# loop lettura dati da seriale (256 righe e due colonne)
runningtime=numpy.zeros(256)
for i in range (0,256):
    data = ard.readline().decode() # legge il dato e lo decodifica
    if data:
        outputFile.write(data) # scrive i 256 dati in file 
        runningtime[i]=int(data[0:data.find(' ')]) # legge i tempi (in unita' di us)
            
outputFile.close() # chiude il file dei dati di carica

ard.close() # chiude la comunicazione seriale con Arduino

deltat=numpy.zeros(255) # crea un array per deltat
for i in range (0,255):
    deltat[i]=runningtime[i+1]-runningtime[i]
deltatavg=numpy.average(deltat) # e lo analizza per trovare la media
deltatstd=numpy.std(deltat) # e la deviazione standard

print("Delta t average = ",deltatavg," us") # che scrive sulla console
print("Delta t stdev = ",deltatstd," us")


print('end \n') # scrive sulla console che ha finito
