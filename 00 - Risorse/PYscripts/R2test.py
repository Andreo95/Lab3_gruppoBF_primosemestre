#test di bontà del fit: "R^2"#
import numpy as np

def f(x):
    return x

x=np.array([])
y=np.array([])
dy=np.array([])
u=0
i=0
D=0
N=0
for u in range(len(y)):#WARNING: in generale ym è il valore di y 'ripulito' in qualche modo dagli errori; se y non è costante sarà diverso dalla media qui fatta!
    ym+=y[i]/len(y)
for i in range(len(y)):
    D+=(y[i]-ym)**2/(dy[i]**2)
    N+=(y[i]-f(x[i])**2/(dy[i]**2)
R=1-N/D
#dà un'indicazione della capacità della funzione f di riprodurre la variazione dei dati osservati#
#più è vicino a 1, meglio il fit riproduce la variazione#