# -*- coding: utf-8 -*-
"""
Created on Fri Jan 22 14:02:50 2021

@author: User
"""

#import relevant libraries
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp

#read data from the file
hdulist = fits.open("Mosaic.fits/mosaic.fits")
hdr=hdulist[0].header #meta data 
data=hdulist[0].data #astro image data
reduced_data = hdulist[0].data
hist_data=data.ravel() #made data file a 1d array of pixel values 
hist_data_sorted=sorted(hist_data) #sort data points into ascending order
background_data=[] #empty list - store the relevant background data
hdulist.close()

#number of rows on the image
imagey = np.shape(data)[0]
#number of columns on the image
imagex = np.shape(data)[1]

#defines the Gaussian function 
#x = position 
#a = peak height 
#mu = avg 
#sig = std dev
def gaussian_fit(x,a,mu,sig):
    gaussian = a*np.exp(-(x-mu)**2/(2*sig**2))
    return gaussian

#bin frequencies and bins from histogram data (one bin for each pixel value)
yhistogram, xhistogram = np.histogram(hist_data, bins=np.max(hist_data))
#popt has optimized variable values 
#fit the scipy function to the data
popt, pcov = sp.optimize.curve_fit(gaussian_fit, xhistogram[0:-1], yhistogram, [8000000, 3420, 20])

# plot the histogram 
plt.ylabel("Bin frequency", fontsize=17)
plt.xlabel("Pixel value", fontsize=17)
plt.tick_params(labelsize=15)
plt.hist(hist_data,bins=np.max(hist_data))
i = np.linspace(3000, 4000, 1001)
# plot the gaussian on top
plt.plot(i, gaussian_fit(i, *popt))
plt.xlim(3350,3550)
plt.ylim(0, 350000)
plt.grid()
plt.show()
#plt.savefig("gaussian_background.jpg")

#N is the number of noisy data points below the estimated mean
N = sum(yhistogram[0: np.abs(xhistogram - popt[1]).argmin()])
muindex = np.abs(xhistogram - popt[1]).argmin()
i = 0
#finds the value of xhistogram that is N datapoints above the mean 
while N > 0 : 
    N -= yhistogram[muindex - 1]
    i +=1
print('the pixel val N points away from the mean', xhistogram[muindex+i])
print('pixel val 1 sig from the mean', xhistogram[np.abs(xhistogram - (popt[1]+ popt[2])).argmin()])
print('pixel val 2 sig from the mean', xhistogram[np.abs(xhistogram - (popt[1]+ 2*popt[2])).argmin()])
print('pixel val 3 sig from the mean', xhistogram[np.abs(xhistogram - (popt[1]+ 3*popt[2])).argmin()])

print('therefore we will consider points only more than 2 sigma away from the mean')

#the results from this file will be quoted to speed up later image processing
print('the mean of our distribution is: ', popt[1], '+/-',pcov[1,1])
print('the standard deviation of our distribution is: ',  popt[2], '+/-', pcov[2,2])
print('the amplitude of our distribution is: ',  popt[0], '+/-',pcov[0,0])
#np.savetxt('image_parameters.txt', popt)
