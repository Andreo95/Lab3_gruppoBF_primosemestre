#stima veloce delle incertezze di calibrazione....

import math
import uncertainties
from uncertainties.umath import *


V_1=uncertainties.ufloat(3.14, 3/100*3.14)
V_2=uncertainties.ufloat(2.04, 3/100*2.04)
f1=1.028e6
f2=0.647e6
A=20*(log10(V_1)-log10(V_2))/math.log10(f1/f2)
print(A)