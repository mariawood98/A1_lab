# -*- coding: utf-8 -*-
"""
Created on Tue Jan 12 16:14:24 2021

@author: maria
"""
no_objects = 20 #how many galaxies
object_loc = []

image_size = 500
import numpy as np
import matplotlib.pyplot as plt
popt = np.loadtxt('image_parameters.txt')

radius1 = 10
basearray = np.zeros((image_size,image_size*2)) #test image, asymmetrical to check indices
#make background with gaussian distribution same as real image
for i in range(image_size):
	for j in range(image_size*2):
		basearray[i][j] = np.random.normal(popt[1],popt[2]) 
		
#assign location to galaxies, but ensuring they're not too close to edges
for i in range(no_objects):
	object_loc.append(((np.random.randint(radius1,(image_size)-radius1)),(np.random.randint(radius1,(image_size*2)-radius1))))
	i += 1
	
def gaussian2d(a, mux, muy, sigma,x, y):
	gauss = a*np.exp(-((x-mux)**2/(2*sigma**2))-((y-muy)**2/(2*sigma**2)))
	return gauss

for i in range(no_objects):
	loctemp=object_loc[i] 
	x_centre=loctemp[0]
	y_centre=loctemp[1]
#	print(radius1)
	for j in range(x_centre-radius1,x_centre+radius1): #6 is the radius, arbitrary can be changed
		for k in range(y_centre-radius1,y_centre+radius1):
#			print(j)
#			print(k)
			radius2 = np.sqrt((j-x_centre)**2+(k-y_centre)**2)
			if radius2<=radius1:
				basearray[j,k] = gaussian2d(5000, x_centre, y_centre, 10, j, k)
			  #makes the objects gaussian blobs, arbitrary amplitude and sigma, can be changed
 
plt.figure()
plt.imshow(basearray)
#plt.savetxt('simulatedimage.txt', basearray)

