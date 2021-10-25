#!/usr/bin/python
# -*- coding: latin-1 -*-

##############################################################################################
# Compares two energy profiles from NEB calculations performed in CDE.
#
# Usage: CompareProfiles.py file1 file2 
#
# Assumes that the relative energies in the third column are the things to compare!
# Output is printed to file 'compare.pdf'
###############################################################################################
          
from pylab import *
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import sys

file1 = str(sys.argv[1])
file2 = str(sys.argv[2])

data1 = genfromtxt(file1,comments="#")
x1 = data1[:,0]
y1 = data1[:,2]
y1 = y1 * 2625.49963

data2 = genfromtxt(file2,comments="#")
x2 = data2[:,0]
y2 = data2[:,2]
y2 = y2 * 2625.49963

#Create Gaussian intensity profile.
#
#fig, ax = plt.subplots(2, sharex=True)

sum = 0.0
nn = len(y1)
for i in range(len(y1)):
    sum = sum + (y1[i]-y2[i])**2

sum = sqrt(sum /nn)
print ('\nRMSE energy difference = ',sum,' kJ/mol\n')

plt.plot(x1,y1,'-o')
plt.plot(x2,y2,'--o')
plt.xlabel('Reaction coordinate',fontsize=16)
plt.ylabel('Energy [kJ/mol]',fontsize=16)
plt.yticks(fontsize=16)
plt.xticks(fontsize=16)
#plt.legend(frameon=False)
plt.savefig("compare.pdf")
plt.show()

