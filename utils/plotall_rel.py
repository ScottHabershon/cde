#!/usr/bin/python
# -*- coding: latin-1 -*-
          
from pylab import *
#import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import sys
import os

ndir = 15 
dirs = []
filn = []
sh = 0 
for i in range(ndir):
  if os.path.exists("final_path_rx_" + str(i+1) + "_opt_profile.dat"):
   dirs.append(i)
   filn.append("final_path_rx_" + str(i+1) + "_opt_profile.dat")
   sh = 2
   continue
  if os.path.exists("final_path_rx_" + str(i+1) + "_opt_profile.dat"):
   dirs.append(i)
   Filn.append("final_path_rx_" + str(i+1) + "_opt_profile.dat")
   if sh != 3:
    sh = 2

maxim = 0
for nam in filn:
  if maxim < len(file(nam,'r').readlines())-sh:
    maxim = len(file(nam,'r').readlines())-sh
nimage = maxim
data = np.zeros((ndir,nimage,3))
maxbarrier = -1000.0

for i, nam in enumerate(filn):
    #str1 = "./final_path_rx_" + str(j+1) + "/input.energy-neb-end"
#    print(nam)
    data1 = genfromtxt(nam,comments="#")
    n = len(data1[:,0])
    data[i,:n,0] = data1[:,0]/data1[n-1,0]
    #data[i,:n,0] = data1[:,0]
    data[i,:n,1] = data1[:,1]
    data[i,:n,2] = data1[:,2]
    #data[i,n:,0] = [ data1[n-1,0] for j in range(nimage-n) ]
    data[i,n:,0] = [ data1[n-1,0]/data1[n-1,0] for j in range(nimage-n) ]
    data[i,n:,1] = [ data1[n-1,1] for j in range(nimage-n) ]
    data[i,n:,2] = [ data1[n-1,2] for j in range(nimage-n) ]
    maxb = max(data1[:,1]-data1[0,1])
    maxbarrier = max(maxb,maxbarrier)
 
#    print 'DATA = ', data[i,:,0]

print 'MAX BARRIER: ',maxbarrier * 2625.5

#Create Gaussian intensity profile.
#
#plt.figure(figsize=(12,6))

fig,ax = plt.subplots(figsize=(12,6))
# = plt.figure(figsize=(12,6))
#ax = fig.add_subplot(111)
ax.set_ylabel("Energy / kJ mol$^{-1}$",fontsize=16)
ax.tick_params(labelsize=16)
ax.axes.get_xaxis().set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)

dataX = np.zeros((ndir,nimage,3))

for i,k in enumerate(dirs):
    if i == 0:
       dataX[i,:,0] = data[i,:,0]
       dataX[i,:,1] = data[i,:,1] - data[0,0,1] 
       dataX[i,:,2] = data[i,:,2] - data[0,0,2]
    else:
       dataX[i,:,0] = data[i,:,0]  
       dataX[i,:,1] = data[i,:,1] - data[0,0,1]  
       dataX[i,:,2] = data[i,:,2] - data[0,0,2]  
       #dataX[i,:,1] = data[i,:,1] + dataX[i-1,nimage-1,1] 
       #dataX[i,:,2] = data[i,:,2] + dataX[i-1,nimage-1,2] 


#plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.purples(np.linspace(0, 1, ndir))))
#plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.brg(np.linspace(0, 1, ndir+1))))

#for i in range(ndir):
for i, k in enumerate(dirs):
    ax.plot(i+data[i,:,0],2625.5*(dataX[i,:,1]),'-',alpha=1.0,markersize=8,linewidth=2)
    ax.plot(i+data[i,:,0],2625.5*(dataX[i,:,1]),'o',alpha=0.6,markersize=8,linewidth=2)

#ax.set_ylim([-300,300])
plt.xlabel('Reaction progress',fontsize=18)
plt.ylabel('Energy [kJ/mol]',fontsize=18)
plt.yticks(fontsize=16)
plt.xticks(fontsize=16)
#plt.legend(frameon=False)
plt.savefig("plotall.pdf")
#plt.show()


