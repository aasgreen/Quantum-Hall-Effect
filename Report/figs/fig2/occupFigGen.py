from numpy import *
import math
import sys
import matplotlib.pyplot as plt
#This python script will generate a graph using the enthought
#ipython distribution's numpy, ipython, and matplotlib.
#For those familiar with matlab, the syntax is very similar.

#Figure 1 for the Report:

fig = plt.figure()
ax = {}
ax[1] = fig.add_subplot(2,1,1)
ax[2] = fig.add_subplot(2,1,2)
#Import data
data = {} #keys will be flux values

data[1] = loadtxt('S0_-3.291subband.dat',unpack=True)
data[2] = loadtxt('S1_-2.0subband.dat', unpack=True)
data[3] = loadtxt('S2_-1.082subband.dat', unpack=True)
#Create figure and tweak so that it looks nice

###FIRST SUBPLOT#####
line, = ax[1].plot(data[1][0], data[1][1], label='First Subband')
line2, = ax[1].plot(data[2][0], data[2][1], label='Second Subband')
line3, = ax[1].plot(data[3][0], data[3][1], label='Third Subband')

xtext = ax[1].set_xlabel('Temperature')
ytext = ax[1].set_ylabel('Occupation Probability')
l1 = ax[1].legend()
#SECOND SUBPLOT###########
line, = ax[2].plot(data[1][0], data[1][1], label='Subband 1')
line2, = ax[2].plot(data[2][0], data[2][1], label='Subband 2')
line3, = ax[2].plot(data[3][0], data[3][1], label='Subband 3')
line4 = ax[2].axvline(x=.24,color='black',ls='--')
text1 = ax[2].text(.25,.7, 'First subbband filled up to here')
l2 = ax[2].legend(('First Subband','Second Subband','Third Subband'),loc=3)
ax[2].set_xlim([0,0.45])
xtext = ax[2].set_xlabel('Temperature')
ytext2 = ax[2].set_ylabel('Occupation Probability')
plt.show()
plt.savefig('occupProb.eps')



