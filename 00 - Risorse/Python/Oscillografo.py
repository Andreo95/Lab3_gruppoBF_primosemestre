import sys, os

sys.path.append( os.path.join(os.path.realpath('..'), '00 - Risorse', 'Python') )
folder = os.path.realpath('.')

import numpy as np, pylab, matplotlib.pyplot as plt, matplotlib as mpl, scipy.stats.distributions as dists
from lab import *
import re

debug=False


class OscilloscopeData():

    def __init__(self,filename, verbose=True, graphicose=True):
        self.source=filename
        file=open(filename)
        t1=[]
        t2=[]
        ch1=[]
        ch2=[]
        self.params={}
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
                pass #TODO: prendere il resto delle informazioni...
                #print("datino")
                # k=re.findall(r"[\w.\s.\+.\-]+", l)
                # try:
                    # self.params[k[0]+"_ch1"]=[]
                    # for i in range(1, len(k)):# pass
                # print(k)
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
        
        

