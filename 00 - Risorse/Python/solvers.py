import sys
import os
import uncertainties as un
import numpy as np
import warnings
sys.path.append(os.path.join(os.path.realpath('.'), '00 - Risorse', 'Python'))
folder = os.path.realpath('.')
from lab import *

"""
R_2 = 17.54e3
R_1 = 178.0e3
R_C = 9.98e3
R_e = 0.989e3
R_2 = un.ufloat(R_2, mme(R_2, unit = "ohm"))
R_1 = un.ufloat(R_1, mme(R_1, unit = "ohm"))
R_C = un.ufloat(R_C, mme(R_C, unit = "ohm"))
R_e = un.ufloat(R_e, mme(R_e, unit = "ohm"))

def Solve(Vcc, Vb, Ve, Vc):
	'''return (Vbe, Vce, Ib, Ic, Ie, I_r1, I_r2)'''
	if(type(Vcc) ==  un.core.Variable):
		Vbe = Vb-Ve
		Vce = Vc-Ve
		I_r1 = (Vcc-Vb)/R_1
		I_r2 = Vb/R_2
		Ic = (Vcc-Vc)/R_C
		Ib = I_r1-I_r2
		Ie = Ve/R_e
		return (Vbe, Vce, Ib, Ic, Ie, I_r1, I_r2)
	else:
		return Solve(un.ufloat(Vcc, mme(Vcc, unit = "volt", metertype = 'oscil')), un.ufloat(Vb, mme(Vb, unit = "volt", metertype = 'oscil')), un.ufloat(Ve, mme(Ve, unit = "volt", metertype = 'oscil')), un.ufloat(Vc, mme(Vc, unit = "volt", metertype = 'oscil')))
"""



def normalize_cinf(x):
	x = np.array(x)
	x[np.abs(x) ==  np.inf] = np.inf
	return x


def parallel(*Zn):
	# Zn = normalize_cinf(Zn)
	invZ = 0
	for Zi in Zn:
		invZ +=  1/Zi
	return 1/invZ

def cross(m, q, c):
	warnings.warn(
		"cross has been superseded by lineheight, use that instead",
		DeprecationWarning
	)
	return  (c-q)/m

def dcross(linepars, linecov, c, dc):
	'''DEPRECATED, use lineheight instead.'''
	warnings.warn(
		"dcross has been superseded by lineheight, use that instead",
		DeprecationWarning
	)
	d = [-(c-linepars[1])/linepars[0]**2, -1/linepars[0]]
	return d @ linecov @ d + dc*d[0]**2

def FindCross(linepars, linecov, c, dc):
	'''DEPRECATED, use lineheight instead.'''
	warnings.warn(
		"FindCross has been superseded by lineheight, use that instead",
		DeprecationWarning
	)
	return cross(linepars[0], linepars[1], c), dcross(linepars, linecov, c, dc)

def lineheight(linepars, linecov, height, hvar=0):
	"""Find intersection (with variance) of straight line and constant height."""
	m, q = linepars
	h = height
	x = (h - q) / m
	derivs = [(q - h)/m**2, -1/m, 1/m]
	covmat = np.zeros((3, 3))
	covmat[:2,:2] = linecov
	covmat[2,2] = hvar
	varx = derivs @ covmat @ derivs
	return x, varx
