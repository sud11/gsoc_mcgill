'''
Program implements sigma clipping and background subtraction

2 Step process

Step 1-Sigma clipping

1.For each pixel,compare it to its values in , say 64 frames and calculate sigma ($) and median (m).
2.We remove all values ( set them to 0)  that are smaller or larger than m +- 4*$. Hence, by not comparing the pixel to the pixels of the same frame, we avoid the risk of losing our star's data during transit.

Step 2-Background subtraction for each frame of the image

1.Mask the central 16 pixels for every frame of the image.
2.Do steps 2 and 3 for frames 0 to 64 exluding the outlier frames (0 and 57)
2.Calculate the median of pixel values.This is the median background flux (bg_avg) of that frame of the image.
3.Subtract bg_avg from each pixel of that frame of image.

**PS: Remove bad frames before sigmaclipping (1st and 58th frame)
'''
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from astropy.utils.data import download_file
from astropy.io import fits
from astropy.stats import sigma_clip
from numpy import std

def sigmaclip(image_data):
	filtered_data = sigma_clip(image_data, sigma=4, iters=2, cenfunc=np.median,axis=0)
	return filtered_data
def bgsubtract(image_data):
	bgsubimg=image_data
	#Creating the mask for central 16 pixels for back
	x=np.ndarray ( shape=(64,32,32), dtype=bool)
	xmask=np.ma.make_mask(x,copy=True, shrink=True, dtype=np.bool)
	xmask[:,:,:]= False
	xmask[:,14:18,14:18]=True
	masked= np.ma.masked_array(image_data, mask=xmask)
	n=0
	#Background subtraction for each frame
	while(n<64):
		#1st frame and 58th frame are outliers
		if(n==0 or n==57):
			n+=1
			continue
		#Excludes NAN values in the calculations
		#bg_avg stores the average flux of the background
		bg_avg=np.median(masked[n,:,:])
		#Subtract the mean 
		bgsubimg[n]= bgsubimg[n,:,:] - bg_avg
		#Setting negative values to 0
	#	bgsubimg[n]= bgsubimg[n].clip(min = 0)
	#	plt.clf()
	#	plt.imshow(bgsubimg[n])
		n+=1
		#plt.savefig(str(n)) #To save the img
	return bgsubimg

# main
filename='SPITZER_I2_20150000_0000_0000_1_bcd'
filepath='SPITZER_I2_20150000_0000_0000_1_bcd.fits'
f=fits.open(filepath)
image_data=f[0].data
filtered_data=sigmaclip(image_data)
bgsubimg=bgsubtract(filtered_data)
#print filtered_data[5,14:18,14:18]
stacked=np.median( filtered_data, axis=0)
plt.imshow(stacked, interpolation='none',cmap='gray')
plt.title('Sigma clipped _'+filename)
plt.xlabel('x pixels')
plt.ylabel('y pixels')
plt.colorbar()
plt.savefig('Sigma clipped _'+filename)
f.close()
