# Esperienza sull'acquisizione di d.d.p. tramite Arduino.
# Lo script scrive sulla porta seriale a cui e' collegato arduino un numero che viene interpretato da Arduino come ritardo in unita' di 10 ms e fa partire l'acquisizione.
# Poi attende l'arrivo dei dati elaborati da Arduino sulla seriale e li salva in un file.

import serial # libreria per gestione porta seriale (USB)
import time   # libreria per temporizzazione

nacqs = 5 # numero di acquisizioni da registrare (ognuna da 600 punti)

Directory='' # nome directory dove salvare i file dati
FileName=(Directory+'data17.txt') # nomina il file dati <<<< DA CAMBIARE SECONDO GUSTO 

outputFile = open(FileName, "w+" ) # apre file dati in scrittura

for j in range (1,nacqs+1):
        print('Apertura della porta seriale\n') # scrive sulla console (terminale)
        ard=serial.Serial('/dev/ttyACM0',9600)  # apre la porta seriale (da controllare come viene denominata, in genere /dev/ttyACM0)
        
        time.sleep(2)   # aspetta due secondi per evitare casini
        
        ard.write(b'1') # intervallo (ritardo) tra le acqusizioni in unita' di 10 ms <<<< questo si puo' cambiare (default messo a 10 ms)
        
        print('Start Acquisition ',j, ' of ',nacqs) # scrive sulla console (terminale)
        
        # loop lettura dati da seriale (sono 600 righe, eventualmente da aggiustare)
        for i in range (0,600):
                data = ard.readline().decode() # legge il dato e lo decodifica
                if data:
                        outputFile.write(data) # scrive i dati nel file
                                
        ard.close() # chiude la comunicazione seriale con Arduino
        
        print('Acquisition ',j,' completed\n') # scrive sulla console (terminale)
        
outputFile.close() # chiude il file dei dati

print('End')
