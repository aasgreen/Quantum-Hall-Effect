from numpy import *
import math
import sys
import matplotlib.pyplot as plt
#This python script will generate a graph using the enthought
#ipython distribution's numpy, ipython, and matplotlib.
#For those familiar with matlab, the syntax is very similar.

#Figure 3 for the Report:
#This figure shows what happens to the kubo formula as a function of the size of the lattice
fig = plt.figure()
ax = {}
ax[1] = fig.add_subplot(1,1,1)
#ax[2] = fig.add_subplot(2,1,2)
#Import data
data = {} #keys will be flux values

data[1] = loadtxt('Q8_Sz2.dat',unpack=True)
data[2] = loadtxt('Q8_Sz8.dat', unpack=True)
data[3] = loadtxt('Q8_Sz10.dat', unpack=True)
data[4] = loadtxt('Q8_Sz20.dat', unpack=True)
data[5] = loadtxt('Q8_Sz30.dat', unpack=True)
spectrum = loadtxt('Q8_Sz10spectrum.dat', unpack=True) #
#Create figure and tweak so that it looks nice

###FIRST FIGURE#####
line, = ax[1].plot(data[1][0], data[1][1], label='Lattice Length=2')
line2, = ax[1].plot(data[2][0], data[2][1], label='Lattice Length=8')
line3, = ax[1].plot(data[3][0], data[3][1], label='Lattice Length=10')
line3, = ax[1].plot(data[5][0], data[4][1], label='Lattice Length=20')
line3, = ax[1].plot(data[5][0], data[5][1], label='Lattice Length=30')

xtext = ax[1].set_xlabel('Fermi Level')
ytext = ax[1].set_ylabel('Kubo Conductance')
l1 = ax[1].legend(loc=4)








###OPTIONAL INSET to show energy spectrum
ax[3] = fig.add_axes([.23,.7,.3,.1])
line = ax[3].scatter(spectrum[3], spectrum[0], s=.2, c ='b', marker = 'o')
ax[3].set_ylim([1./8.*(1-1/100.), 1/8.*(1+1/100.)])
ax[3].yaxis.set_ticks( (1/8.,) )
ytext_inset = ax[3].set_ylabel('Flux')
xtext_inset = ax[3].set_xlabel('Energy Levels')
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
plt.savefig('kScalingKubo.eps')



