# -*- coding: utf-8 -*-
"""
Created on Sun Jan 24 19:13:24 2021

@author: User
"""
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

#read data from the file
hdulist = fits.open("Mosaic.fits/mosaic.fits")#plt.imshow(data, cmap='gray')
hdr=hdulist[0].header #meta data 
data=hdulist[0].data #astro image data
background_data=[] #empty list - store the relevant background data
hdulist.close()

datatest = data#[1500:1650,1900:2000] #[y,x]
data = data#[1500:1650,1900:2000] #[y,x]

#show the test image taken 
#plt.figure()
#plt.imshow(np.log10(datatest))
#plt.show()
#dimensions of the image 
testimagey = np.shape(datatest)[0]
testimagex = np.shape(datatest)[1]
#load in image parameters and the premade mask 
popt = np.loadtxt('image_parameters.txt')
mask = np.loadtxt('mask_2sigma.txt')
#this is the sample of the mask that corresponds to the image
masktest =  mask#[1500:1650,1900:2000] #[y,x]
objectmask = np.ones((testimagey, testimagex))
#displays the mask sample to show is corresponds to the image
#plt.figure()
#plt.imshow(masktest)
#plt.show()

datatest = datatest*masktest
#plt.figure()
#plt.imshow(np.log10(datatest))
#plt.show()
#make the image dataset 1D
datatest1d=datatest.ravel()
#sort the data into ascending order
datatest1d_sorted = sorted(datatest1d)
datatest1d_sorted = list(filter(lambda a: a!=0, datatest1d_sorted ))

#sums values of pixels inside 1st aperture, pixels of galaxy
def galaxyPhotons(x_centre,y_centre, radius):
    photon_vals=[];
    #take i and j to scan across the aprture index values
    for i in range(x_centre-radius,x_centre+radius):
        for j in range(y_centre-radius,y_centre+radius):
            #if the index is still inside our sample image and the mask is 1 (i.e. readable) 
            if i<testimagex and i>=0 and j<testimagey and j>=0 and objectmask[j,i] == 1: #only consider pixels within the image!
                #d is the distance fromthe image 
                d=np.sqrt((i-x_centre)**2+(j-y_centre)**2)
                #if d is smaller than the radius
                if d<radius:
                    #add this photon to the list of photon values
                    photon_vals.append(data[j,i])
                    #scanned pixels get marked off so they aren't rescanned
                    masktest[j,i] = 0 
                    objectmask[j,i] = 0
    #return a sum of all the photon values and the number of pixels scanned
    return(np.sum(photon_vals), len(photon_vals))

#averages values of pixels between 1st and 2nd aperture, local bakcground mean
def localBackground(x_centre,y_centre, initialRadius, secondaryRadius):
    bg_photon_vals=[];
    for i in range(x_centre-secondaryRadius,x_centre+secondaryRadius):#outer x
        for j in range(y_centre-secondaryRadius,y_centre+secondaryRadius):#outer y
            if i<testimagex and i>=0 and j<testimagey and j>=0: #only consider pixels within the image!
                #diameter of the secondary apeture
                d2=np.sqrt((i-x_centre)**2+(j-y_centre)**2)
                #take only points between the two rings 
                if initialRadius <= d2 and d2<=secondaryRadius and objectmask[j,i] == 1:
                    #the pixel is added to the list 
                    bg_photon_vals.append(data[j,i])
                    masktest[j,i] = 0 
    #return the average pixel value from the outer apeture
    return(np.mean(bg_photon_vals), np.std(bg_photon_vals))

def pixelPos(data, value):
    indices = np.where(datatest==value) #gives tuples of coordinates
    listIndices= list(zip(indices[0], indices[1])) # gives pairs of coordinates in form (x,y)
    #print(listIndices)
    return listIndices

def initialRad(counts):
    radius = np.log(counts - (popt[1]+4*popt[2]))*1 #set the value in the brackets based on the lower threshold
    return radius 


# initialRadius = 6
# secondaryRadius = initialRadius*2.5
galaxyCounts = []
galaxyCountsErr = []
galaxyLocation = []
initial_rad_list = []
#ratio between r1 and r2 
q = 2.5
for i in range(len(datatest1d_sorted)): 
    tempPixVal = datatest1d_sorted[-(1+i)] #find highest pixel value
    if tempPixVal > popt[1]+4*popt[2]:
        tempPixPos = pixelPos(datatest, tempPixVal) #find position of highest pixel value
        for j in range(len(tempPixPos)): #potentially multiple locations with same pixel value
            loc = tempPixPos[j]
            if masktest[loc[0],loc[1]] == 1 and objectmask[loc[0],loc[1]] == 1: #if pixel is available to use 
                initialRadius = int(initialRad(tempPixVal)) #calculate the radius based on the pixel val
                if initialRadius >= 3: #exclude radii smaller than 3
                    initial_rad_list.append(initialRadius)
                    #reads the flux sum and number of pixels in the galaxy
                    galaxyBrightness=galaxyPhotons(loc[1],loc[0],initialRadius)
                    secondaryRadius = int(initialRadius*q) #calculates the outer radius
                    #reads the avg background pixel val and the std dev of the background
                    galaxyBackground = localBackground(loc[1],loc[0],initialRadius, secondaryRadius)
                    #adds the corrected flux value of the galaxy to a list 
                    galaxyCounts.append(galaxyBrightness[0]-(galaxyBackground[0]*galaxyBrightness[1]))
                    #adds the error of the galaxy flux to a list
                    galaxyCountsErr.append(galaxyBackground[1]*galaxyBrightness[1])
                    galaxyLocation.append(loc) #stores the position of the galaxy
                    print(len(datatest1d_sorted) -i) #counter for n. iterations remaining

def Magnitude(galaxyflux, galaxyflux_err):
    m = []
    m_err = []
    ZPinst = 25.3
    ZPinst_err = 2e-2
    for i in range(len(galaxyflux)):
        m.append(ZPinst - 2.5*np.log10(galaxyflux[i]))
        merr = np.sqrt(ZPinst_err**2 + (2.5**2)*(0.434**2)*((galaxyflux_err[i]/galaxyflux[i])**2))
        m_err.append(merr)
    return(m, m_err)
    
galaxyMagnitude = Magnitude(galaxyCounts, galaxyCountsErr)



plt.figure()
plt.imshow(np.log10(data), cmap='jet')
p = 0		
for i in galaxyLocation:
    y = i[0]
    x = i[1]
    i = [x, y]
    mycircle = plt.Circle(i, initial_rad_list[p], color= 'r',fill=False)
    plt.gca().add_artist(mycircle)
    mycircle2 = plt.Circle(i, initial_rad_list[p]*q, color= 'r',fill=False)
    plt.gca().add_artist(mycircle2) 
    p+=1
plt.show()

plt.figure()
plt.imshow(np.log10(datatest), cmap='jet')
p = 0		
for i in galaxyLocation:
    y = i[0]
    x = i[1]
    i = [x, y]
    mycircle = plt.Circle(i, initial_rad_list[p], color= 'k',fill=False)
    plt.gca().add_artist(mycircle)
    mycircle2 = plt.Circle(i, initial_rad_list[p]*q, color= 'k',fill=False)
    plt.gca().add_artist(mycircle2) 
    p+=1
plt.show()
plt.figure()
plt.imshow(masktest)
plt.show()

##np.savetxt("variable_r galaxy count err.txt", galaxyCountsErr)
#np.savetxt("variable_r galaxy locations.txt", galaxyLocation)
#np.savetxt("variable_r radius lengths.txt", initial_rad_list)
#np.savetxt("variable_r galaxy magnitudes", galaxyMagnitude[0])
#np.savetxt("2sigma variable_r galaxy magnitude errors", galaxyMagnitude[1])