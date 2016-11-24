import uncertainties
from uncertainties import ufloat

g=ufloat(4.4e-3, 0.08e-3)
R2=ufloat(227, 3)
R1=ufloat(551, 5)
A1=-g*R1/(1+g*R2)
print(A1)
A2=g*R2/(1+g*R2)
print(A2)

