# To justify why 1st and 58 frames are outliers
# A comparison of background electron count per pixel in different frames for a sample of data cubes

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches
import time
import os, sys
from astropy.io import fits
from astropy.stats import sigma_clip
from photutils import aperture_photometry
from photutils import CircularAperture
from numpy import std
import glob
import csv
path='/home/hema/Documents/a_photometry'
xn=1
for filename in glob.glob(os.path.join(path, '*bcd.fits')):
	print filename
	f=fits.open(filename)
	image_data=f[0].data
	x=np.ndarray( shape=(64,32,32), dtype=bool)
	xmask=np.ma.make_mask(x,copy=True, shrink=True, dtype=np.bool)
	xmask[:,:,:]= False
	xmask[:,13:18,13:18]=True
	masked= np.ma.masked_array(image_data, mask=xmask)
	convfact=f[0].header['GAIN']*f[0].header['EXPTIME']/f[0].header['FLUXCONV']
	image_data1=masked*convfact
	print convfact
	#print image_data1[2,10:18,10:18]
	bgsum = np.zeros(64)
	plt.xlabel('Frame number')
	plt.ylabel('Background electron count per pixel')
	for i in range (64):
		bgsum[i] = np.nansum(image_data1[i,:,:])/(32*32-25)
		#print i, bgsum[i]
	frameno=np.arange(1,65)
	plt.plot(frameno,bgsum)
	xn+=1
	plt.savefig(str(xn)+'.png')

