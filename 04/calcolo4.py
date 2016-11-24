from uncertainties import *
from lab import *

alpha=ufloat(77.20605946, 2.85542836**0.5)
hfe=ufloat(178, 3)
R_C=ufloat(9.98e3, mme(9.98e3, unit="ohm"))
R_e=ufloat(989, mme(989, unit="ohm"))
R_res=ufloat(100, mme(100, unit="ohm"))

def par(x, y):
    return x*y/(x+y)

hie=R_C*hfe/alpha - (hfe+1)*par(R_e, R_res)
print(hie)