# -*- coding: utf-8 -*-
"""
Created on Sun Jan 24 15:58:32 2021

@author: maria
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp

# =============================================================================
#  import necessary data
# =============================================================================
galaxyCounts=np.loadtxt('./q = 2.5, min r = 3, local BG correction/variable_r galaxy counts.txt')
magnitudeError=np.loadtxt('./q = 2.5, min r = 3, local BG correction/Variable_r galaxy magnitude errors')
#galaxyRadii = np.loadtxt('./q = 2.5, min r = 3, local BG correction/variable_r radius lengths.txt')

ZP=2.530E+01 #zero point calbration taken from header
ZP_err = 2.0E-2 #taken from header of FITS file

magnitudes = []
magnitudeError2 = [] #errors corresponding to positive fluxes
for i in range(len(galaxyCounts)):
	if galaxyCounts[i]>0: #filter out negative fluxes 
		tempMag=ZP-2.5*np.log10(galaxyCounts[i])
		magnitudes.append(tempMag)
		magnitudeError2.append(magnitudeError[i])

bins1=list(np.arange(np.min(magnitudes), np.max(magnitudes), 0.2)) #bins for placing the magnitudes into
magnitudesZip = list(zip(magnitudes, magnitudeError2)) #keeping the magnitudes and their errors together
magnitudesSorted = sorted(magnitudesZip) #sorted according to magnitude
magValSorted=sorted(magnitudes) #sorted magnitudes, without associated error

magValsMinimum = [] #minimum magnitude, based on the error
magValsMaximum = [] #maximum magnitude, based on the error
for i in range(len(magnitudesSorted)):
	tempMin = magnitudesSorted[i][0]-magnitudesSorted[i][1]
	tempMax = magnitudesSorted[i][0]+magnitudesSorted[i][1]
	magValsMinimum.append(tempMin)
	magValsMaximum.append(tempMax)
magValsMinimum = sorted(magValsMinimum)
magValsMaximum = sorted(magValsMaximum)
	
numberSmallerMag = [] #cumulative number of objects with a given magnitude or lower
numberSmallerMag_max = [] #number in each bin if magnitude has lowest value
numberSmallerMag_min = [] #number in each bin if magnitude has highest value

for i in range(len(bins1)):
	a_max=np.where(magValsMaximum<=bins1[i])
	b_max=len(a_max[0])
	numberSmallerMag_max.append(b_max)
	a_min=np.where(magValsMinimum<=bins1[i])
	b_min=len(a_min[0])
	numberSmallerMag_min.append(b_min)

for i in range(len(bins1)):
	a=np.where(magValSorted<=bins1[i])
	b=len(a[0])
	numberSmallerMag.append(b)

binError = [] #difference in bin population when magnitude at highest vs lowest values
for i in range(len(numberSmallerMag)):
	errorTemp = numberSmallerMag_min[i] - numberSmallerMag_max[i]
	binError.append(0.434/errorTemp)

poissonError = [] #sqrt(N) error associated with counting pixels
for i in range(len(numberSmallerMag)):
	poissonError.append(0.434/np.sqrt(numberSmallerMag[i]))
	
radiusBins = []
for i in bins1:
	radiusBins.append((i-25.3)/(-2.5*0.434))

# =============================================================================
# ### plotting
# =============================================================================
#set 3rd argument to binError or poissonError
#set 1st argument to radiusBins to check radius behaviour
plt.errorbar(bins1, np.log10(numberSmallerMag), poissonError, marker='x', markersize=9,fmt=' ') #plotting binned data with errorbars
plt.grid()
plt.xlabel("Magnitude", fontsize=15)
plt.ylabel("log(N)", fontsize=15)
plt.tick_params(labelsize=12)


# =============================================================================
# ### linear function for fitting
# =============================================================================
def linear(x, m, c):
	return((m*x)+c)
 
#choose indices of linear region manually by examining plot
linearStartIndex = 7
linearEndIndex = 34


#set sigma to binError or poissonError
p,q=sp.optimize.curve_fit(linear, bins1[linearStartIndex:linearEndIndex], np.log10(numberSmallerMag[linearStartIndex:linearEndIndex]), sigma=(binError)[linearStartIndex:linearEndIndex], absolute_sigma = True)  #perform a linear fit

fitPoints = [] #y values for the linear fit plot
for i in range(linearStartIndex, linearEndIndex):
    fitPoints.append(p[0]*bins1[i]+p[1])
	
plt.plot(bins1[linearStartIndex:linearEndIndex],fitPoints, 'r',linewidth =2, label="N = "+str(p[0])+"m + "+str(p[1])) #plot the linear fit on top of data
#plt.legend()
plt.vlines(bins1[linearStartIndex],0,3.5,'r', 'dashed')
plt.vlines(bins1[linearEndIndex],0,3.5, 'r','dashed')
plt.xlim(8,20)
print("the gradient of the linear part is ",p[0], "+/-",np.sqrt(q[0,0]))
