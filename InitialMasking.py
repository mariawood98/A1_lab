# -*- coding: utf-8 -*-
"""
Created on Fri Jan 22 16:42:35 2021

@author: User
"""
#import relevant libraries
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

#read data from the file
hdulist = fits.open("Mosaic.fits/mosaic.fits")
#plt.imshow(data, cmap='gray')
hdr=hdulist[0].header #meta data 
data=hdulist[0].data #astro image data
hdulist.close()

#number of rows
imagey = np.shape(data)[0]
#number of columns
imagex = np.shape(data)[1]

popt = np.loadtxt('image_parameters.txt')
#consider 1s don't consider 0s
mask = np.ones((imagey,imagex)) #initiate a mask image 

#masking of foreground pixels including stars, bleeding and noisy edges
for i in range (imagey):
    for j in range(imagex):
        #masking the borders
        if i <= 200 or j <= 200 or j >= 2370 or i>= 4411: 
            mask[i,j] = 0
        # masking pixel values below 2 std devs from the mean 
        if data[i][j] <= popt[1] + 2*popt[2]:
            mask[i][j] = 0
        #mask oversaturated pixels
        if data[i][j] >= 50000:
            mask[i][j] = 0

#masking central vertical bloom
mask[:,1425:1447] = 0
#mask horizontal blooming
mask[415:480,1195:1655] = 0
mask[310:380,1020:1705] = 0
mask[200:285,1390:1480] = 0
#masking misc blooming:
mask[3370:7420,770:784] = 0
mask[3200:3279,770:784] = 0
#masking the stars (identified by eye)
mask[420:460,1025:1050] = 0
mask[2220:2240,900:910] = 0
mask[4320:4500, 200:400] = 0
mask[3740:3660,2120:2140] = 0
mask[420:440, 1100:1200] = 0
mask[2940:3900, 180:210] = 0
mask[4000:4030, 1450:1470] = 0
mask[2820:3000, 1400:1470] = 0

mask[2900:3525,1130:1755] = 0
mask[2225:2370,855:980] = 0
mask[2700:2855,915:1035] = 0
mask[3195:3425,710:830] = 0
mask[3710:3815,2095:2180] = 0
mask[3265:3340,2220:2305] = 0
mask[2265:2350,2090:2175] = 0
mask[1370:1465,2040:2140] = 0
mask[1745:1815,1375:1450] = 0
mask[4015:4065,1425:1490] = 0
mask[2250:2310,2280:2325] = 0
mask[3820:3875,2245:2305] = 0
mask[3815:3880,2240:2310] = 0
mask[2275:2330,425:475] = 0
mask[1465:1520,605:670] = 0
mask[540:605,1730:1815] = 0
mask[4065:4135,525:590] = 0
plt.figure()
plt.imshow(mask)
#plt.figure()
#plt.imshow(np.log10(data))
print('we have this many data points: ', np.count_nonzero(mask))
#np.savetxt('mask_2sigma.txt', mask)
