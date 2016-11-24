/* 
Questo sketch serve per acquisire forme d'onda dipendenti dal tempo in modalita'
average. La porta analogica letta è A0 e il trigger digitale e' sulla porta 5.
In ogni acquisizione vengono acquisiti 256 punti. Lo script di Python tratta i dati 
in modo da determinare il valore medio e la deviazione standard sul campione di misure
*/

// Blocco definizioni
const unsigned int analogPin=0; // Definisice la porta A0 per la lettura
const unsigned int sincPin = 5;  //pin 5 ingresso digitale per la sincronizzazione con il generatore
int i; // Definisice la variabile intera i (contatore)
int delays; // Definisce la variabile intera delays
int V[256]; // Definisce l'array intero V
long t[256]; // Definisice l'array t
unsigned long StartTime; // Definisce il valore StartTime
int start=0; // Definisce il valore start (usato come flag)
int sinc; //variabile di sincronizzazione

// Istruzioni di inizializzazione
void setup()
  {
   Serial.begin(19200); // Inizializza la porta seriale a 19200 baud
   Serial.flush(); // Pulisce il buffer della porta seriale 
   pinMode(sincPin, INPUT);  //pin sincPin configurato come ingresso digitale
   analogReference(INTERNAL); // Sceglie il riferimento V_ref = 1.1 V (nominali)
   bitClear(ADCSRA,ADPS0); // Istruzioni necessarie per velocizzare
   bitClear(ADCSRA,ADPS2); // il rate di acquisizione analogica
  }

// Istruzioni del programma
void loop()
  {
    if (Serial.available() >0) // Controlla se il buffer seriale ha qualcosa
      {
        delays = (Serial.read()-'0')*10; // Legge il byte e lo interpreta come ritardo
        Serial.flush(); // Svuota la seriale
	start=1; // Pone il flag start a uno
      }
 
  if(!start) return // Se il flag e' start=0 non esegue le operazioni qui di seguito
                    // altrimenti le fa partire (quindi aspetta di ricevere l'istruzione
                    // di partenza
    delay(1000); // Aspetta 1000 ms per evitare casini
    sinc = digitalRead(sincPin);//legge sincPin
    while (sinc==HIGH) // ciclo di attesa iniziale per sincronizzazione, attende che sincPin vada basso
      {sinc = digitalRead(sincPin);} //legge sincPin
    while (sinc==LOW) // ciclo di attesa iniziale per sincronizzazione, attende che sincPin vada altoo
      {sinc = digitalRead(sincPin);}//legge sincPin    
     StartTime=micros(); // Misura il tempo iniziale con l'orologio interno
     for(i=0;i<2;i++) // Fa un ciclo di due letture a vuoto per "scaricare" l'analogPin
       {
        V[i]=analogRead(analogPin);
       }
    for(i=0;i<256;i++) // Loop di misura
      {
          t[i]=micros()-StartTime; // Legge il timestamp e lo mette in array t
          V[i]=analogRead(analogPin); // Legge analogPin e lo mette in array V
          delayMicroseconds(delays); // Aspetta tot us
      }  
    for(i=0;i<256;i++) // Loop per la scrittura su porta seriale
      {
        Serial.print(t[i]); // Scrive t[i]
        Serial.print(" "); // Mette uno spazio
        Serial.println(V[i]); // Scrive V[i] e va a capo
      }
 
    start=0; // Annulla il flag
    Serial.flush(); // Pulsice il buffer della porta seriale (si sa mai)
}

