# A1_lab
This project takes a CCD image of a Deep Galaxy Survey to find the relationship between galaxy magnitude and cumualtive number of galaxies with a given magnitude.
To run the files, the mosaic.fits data file is needed. It is also necessary to run the numbered files blow in the correct order, as each generates a .txt file used in the next.
## Files and their functions

1. The ImageStatistics file reads in the image data and determines the mean and standard deviation of the background by fitting a Gaussian distribution. It outputs a .txt file with the mean and standard deviation to be used in the other files.
2. The InitialMasking file produces a mask array of ones the same size as the data, with undesirable pixels set to a value of zero. It removes bleeding pixels, birght objects, noisy areas and background. The output is a .txt file with a 2D array of the inital mask, to be used in GalaxyCountsProgram.
3. The GalaxyCountsProgram file carries out the analysis of the image file. This is where the galaxies are identified, their fluxes counted, corrected to background and then masked. This outputs a number of .txt files with the counts, location, magnitude and errors for each image analysed.
4. The DataAnalysis file plots the final galaxy catalogue graph from the data collected in GalaxyCountsProgram and fits a straight line to the area where the lienar relation holds.
* The GaussianSimulator files generates Gaussian "blob" galaxies with random locations but fixed radii and flux to test the basic functionality of the code. It outputs a .txt 2D array image to be used in GalaxyCountsProgram for testing.
