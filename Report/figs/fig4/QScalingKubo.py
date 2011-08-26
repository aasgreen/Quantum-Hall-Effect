from numpy import *
import math
import sys
import matplotlib.pyplot as plt
#This python script will generate a graph using the enthought
#ipython distribution's numpy, ipython, and matplotlib.
#For those familiar with matlab, the syntax is very similar.

#Figure 3 for the Report:
#This figure shows what happens to the kubo formula as a function of the size of the flux
fig = plt.figure()
ax = {}
ax[1] = fig.add_subplot(1,1,1)
#ax[2] = fig.add_subplot(2,1,2)
#Import data
data = {} #keys will be flux values

data[1] = loadtxt('Q8_Sz8.dat',unpack=True)
data[2] = loadtxt('Q16_Sz8.dat', unpack=True)
data[3] = loadtxt('Q25_Sz8.dat', unpack=True)
data[4] = loadtxt('Q64_Sz8.dat', unpack=True)
#data[5] = loadtxt('Q128_Sz2.dat', unpack=True)
#data[6] = loadtxt('Q256_Sz2.dat', unpack=True) #
#Create figure and tweak so that it looks nice

###FIRST FIGURE#####
line, = ax[1].plot(data[1][0], data[1][1], label='Q=8')
line2, = ax[1].plot(data[2][0], data[2][1], label='Q=16')
line3, = ax[1].plot(data[3][0], data[3][1], label='Q=25')
line4, = ax[1].plot(data[4][0], data[4][1], label='Q=64')
#line5, = ax[1].plot(data[5][0], data[5][1], label='Q=128')
#line6, = ax[1].plot(data[6][0], data[6][1], label='Q=256')
xtext = ax[1].set_xlabel('Fermi Level')
ytext = ax[1].set_ylabel('Kubo Conductance')
l1 = ax[1].legend(loc=4)

#SECOND Figure###########
#This figure will show the scaling in Q that occurs.
#line, = ax[2].plot(data[1][0], data[1][1], label='Subband 1')
#line2, = ax[2].plot(data[2][0], data[2][1], label='Subband 2')
#line3, = ax[2].plot(data[3][0], data[3][1], label='Subband 3')
#line4 = ax[2].axvline(x=.24,color='black',ls='--')
#text1 = ax[2].text(.25,.7, 'First subbband filled up to here')
#l2 = ax[2].legend(('First Subband','Second Subband','Third Subband'),loc=3)
#ax[2].set_xlim([0,0.45])
#xtext = ax[2].set_xlabel('Temperature')
#ytext2 = ax[2].set_ylabel('Occupation Probability')
#plt.show()
plt.savefig('occupProb.eps')



