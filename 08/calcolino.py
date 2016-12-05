#da eseguire dopo sifannocose
import uncertainties
import uncertainties.umath

C1, C2=uncertainties.correlated_values(ampliphase.pars, ampliphase.cov)
C1=C_1*C1
C2=C_2*C2
R1=uncertainties.ufloat(R_1, mme(R_1, unit="ohm"))
R2=uncertainties.ufloat(R_2, mme(R_2, unit="ohm"))
f=(2*np.pi*uncertainties.umath.sqrt(C1*C2*R_1*R2))**-1