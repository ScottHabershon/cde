#!/usr/bin/python
# -*- coding: latin-1 -*-
          
from pylab import *
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import sys

ndir = 15 
nimage = 10
data = np.zeros((ndir,nimage,2))
maxb = -1000.0
imax = -1
for i in range(ndir):
    str1 = "./path_fix_pt_" + str(i+1) + "/input_neb.energy-neb-end"
    data1 = genfromtxt(str1,comments="#")
    data[i,:,0] = data1[:,0]/data1[nimage-1,0]
    data[i,:,1] = data1[:,1]
    maxbar = max(data1[:,1]-data1[0,1])
    if maxbar > maxb:
       maxb = maxbar
       imax = i
print 'MAX BARRIER: ',imax,maxb * 2625.5

fig,ax = plt.subplots(figsize=(12,6))
ax.set_ylabel("Energy / kJ mol$^{-1}$",fontsize=16)
ax.tick_params(labelsize=16)
ax.axes.get_xaxis().set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)


#plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.purples(np.linspace(0, 1, ndir))))
#plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.brg(np.linspace(0, 1, ndir+1))))

for i in range(ndir):
    ax.plot(i+data[i,:,0],2625.5*(data[i,:,1]-data[0,0,1]),'-',alpha=1.0,markersize=8,linewidth=2)
    ax.plot(i+data[i,:,0],2625.5*(data[i,:,1]-data[0,0,1]),'o',alpha=0.6,markersize=8,linewidth=2)

#ax.set_ylim([-300,300])
#plt.xlabel('Reaction progress',fontsize=18)
plt.ylabel('Energy [kJ/mol]',fontsize=18)
plt.yticks(fontsize=16)
plt.xticks(fontsize=16)
#plt.legend(frameon=False)
plt.savefig("plotall.pdf")
plt.show()

