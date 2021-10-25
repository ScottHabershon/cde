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
filn1 = []
filn2 = []
sh = 0
runid = 101 
#str1 = "./Benzene_Set_3/run_" + str(runid) +"/NEB/"
str1 = "./MCP_Set_3/run_" + str(runid) +"/NEB/"
#str2 = "./Benzene_Set_3_idpp/run_" + str(runid) + "/NEB/"
str2 = "./MCP_Set_3_idpp/run_" + str(runid) + "/NEB/"
str3 = "path_fix_pt_" 
str4 = "/input_neb.energy-neb-end"

# Get filename lists.
for i in range(ndir):
  t1 = False
  t2 = False
  if os.path.exists(str1 + str3 + str(i+1) + str4):
   t1 = True
  if os.path.exists(str2 + str3 + str(i+1) + str4):
   t2 = True
  print t1,t2

  if t1 is True and t2 is True:
   dirs.append(i)
   filn1.append(str1 + str3 + str(i+1) + str4)
   filn2.append(str2 + str3 + str(i+1) + str4)

print filn1

# Populate with lowest barriers.
filn = []
ilen = len(filn1)
for i in range(ilen):
    data1 = genfromtxt(filn1[i],comments="#")
    data2 = genfromtxt(filn2[i],comments="#")
    max1 = max(data1[:,1]-data1[0,1])
    max2 = max(data2[:,1]-data2[0,1])
    if max1 < max2:
       filn.append(filn1[i])
    else:
       filn.append(filn2[i])

maxim = 0
for nam in filn:
  if maxim < len(file(nam,'r').readlines())-sh:
    maxim = len(file(nam,'r').readlines())-sh
nimage = maxim
data = np.zeros((ndir,nimage,3))

maxb = -1000.0
avebar = 0.00
for i, nam in enumerate(filn):
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
    maxbar = max(data1[:,1]-data1[0,1])
    maxb = max( maxbar, maxb )
    avebar = avebar + maxbar

print 'MAXIMUM BARRIER: ',maxb * 2625.5
print 'AVERAGE BARRIER: ',(avebar/len(filn)) * 2625.5


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
plt.show()



