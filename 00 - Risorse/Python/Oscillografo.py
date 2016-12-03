import sys, os

sys.path.append( os.path.join(os.path.realpath('..'), '00 - Risorse', 'Python') )
folder = os.path.realpath('.')

import numpy as np, pylab, matplotlib.pyplot as plt, matplotlib as mpl, scipy.stats.distributions as dists
from lab import *
import re

debug=False
reg=re.compile(r"(?P<nome>[\w]+\s?[\w]*),(?P<dato>[\w\+\-\.]+)")
reg2=re.compile(r"(?P<nome>[\w]+\s?[\w]*),(?P<dato>[\w\+\-\.]+[\.*,\.*])+[,\n]")

def tryparse(s):
    try:
        return float(s)
    except:
        print("string \"", s, "\" is not a float, is it?")
        return s

class OscilloscopeData():

    def __init__(self,filename, verbose=True, graphicose=True, getall=False):
        self.source=filename
        file=open(filename)
        t1=[]
        t2=[]
        ch1=[]
        ch2=[]
        self.CH1params={}
        self.CH2params={}
        for l in file.readlines():
            try:
                w=re.findall(r"[\+\-\w\.]+",l)
                t1.append(float(w[0]))
                ch1.append(float(w[1]))
                t2.append(float(w[2]))
                ch2.append(float(w[3]))
                if debug:
                    print("Trovata roba in ", l, " --->", (w[0], float(w[1]), float(w[2]), float(w[3])))
            except:
              #  print("datino")
                if(getall):
                    k=reg.findall(l)
                    print(k)
                    for i in k:
                        if(i[0] in self.CH1params):
                            self.CH2params[i[0]]=i[1]#tryparse(i[1])
                        else:
                            self.CH1params[i[0]]=i[1]#tryparse(i[1])

        self.T1=np.array(t1)
        self.T2=np.array(t2)
        self.CH1=np.array(ch1)
        self.CH2=np.array(ch2)
        self.dCH1=mme(np.amax(np.abs(self.CH1)), "volt", "oscil")
        self.dCH2=mme(np.amax(np.abs(self.CH2)), "volt", "oscil")
    
    def plot(self):
        pylab.figure(0)
        pylab.errorbar(self.T1, self.CH1, self.dCH1)
        pylab.figure(1)
        pylab.errorbar(self.T2, self.CH2, self.dCH2)
        pylab.show()
        
        

