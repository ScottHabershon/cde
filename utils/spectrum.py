#!/usr/bin/python
# -*- coding: latin-1 -*-
          
from pylab import *
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

data = genfromtxt('spectrum.dat',comments="#")
print(data)
xs = data[:,0]
ys = data[:,1]


#Create Gaussian intensity profile.
#
sigma = 10.0
alpha = 1.0 / (2.0 * sigma * sigma)
npts = 250
xmin = 200
xmax = 800
dx = (xmax-xmin)/float(npts-1)
ysc = np.empty(npts)
ytc = np.empty(npts)

xx = np.empty(npts)
for i in range(npts):
    xx[i] = xmin + (i-1)*dx
    
    ysc[i] = 0.0
    for j in range(len(ys)):
        ysc[i] += ys[j] * exp(-alpha*(xx[i]-xs[j])**2)
        

# Plotting.
#
#fig, ax = plt.subplots(2, sharex=True)

plt.plot(xx,ysc)
plt.xlabel('Wavelength [nm]',fontsize=18)
plt.ylabel('Oscillator strength',fontsize=18)
plt.yticks(fontsize=18)
plt.xticks(fontsize=18)
#plt.legend(frameon=False)
plt.savefig("spectrum.pdf")
plt.show()

