# To justify why 1st and 58 frames are outliers
# Do an aperture photometry and verify.
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

def sigma_clipping(image_data):#,fname):
	#global badframetable
	#global tossed
	sig_clipped_data=sigma_clip(image_data, sigma=4, iters=2, cenfunc=np.median,axis=0)
	return sig_clipped_data

def bgsubtract(image_data):
	bgsubimg=image_data
	x=np.ndarray ( shape=(64,32,32), dtype=bool)
	xmask=np.ma.make_mask(x,copy=True, shrink=True, dtype=np.bool)
	xmask[:,:,:]= False
	xmask[:,14:18,14:18]=True
	#xmask[:,0:1,:]=True
	masked= np.ma.masked_array(image_data, mask=xmask)
	n=0
	#Background subtraction for each frame
	while(n<64):
		bg_avg=np.mean(masked[n,:,:])
		bgsubimg[n]= bgsubimg[n,:,:] - bg_avg
		n+=1
	return bgsubimg
def centroid(image_data):
	# Refer: Intra-Pixel Gain Variations and High-Precision Photometry with the Infrared Array Camera (IRAC)
	cx=np.zeros(64)
	cy=np.zeros(64)
	starbox = image_data[:, 13:18, 13:18]
	h,w = np.shape(starbox[0,:,:])
	x = np.arange(0,w)
	y = np.arange(0,h)
	X,Y = np.meshgrid(x,y)
	for i in range(64):
		cx[i]=(np.sum(X*starbox[i,:,:])/np.sum(starbox[i,:,:]))+13
		cy[i]=(np.sum(Y*starbox[i,:,:])/np.sum(starbox[i,:,:]))+13
	return cx,cy

def aperphot(image_data,ape_radius,cx,cy):
	ape_sum=np.zeros(64)
	for i in range(64):
		position=[cx[i],cy[i]]
		aperture=CircularAperture(position,r=ape_radius)
		phot_table=aperture_photometry(image_data[i,:,:],aperture)
		temp=phot_table['aperture_sum']
		ape_sum[i]=phot_table['aperture_sum']
	return ape_sum



path='/home/hema/Documents/mcgill/handy/datasample/'
#xn=1
#for filename in glob.glob(os.path.join(path, '*bcd.fits')):
ct=1
for filename in glob.glob(os.path.join(path, '*bcd.fits')):
	f=fits.open(filename)
	image_data0=f[0].data
	# convert MJy/str to electron count
	convfact=f[0].header['GAIN']*f[0].header['EXPTIME']/f[0].header['FLUXCONV']
	image_data1=image_data0#*convfact
	#sigma clip
	image_data2=sigma_clipping(image_data1)
	#bg subtract
	image_data3=bgsubtract(image_data2)

	#centroid
	xo, yo = centroid(image_data3)
	#aperture photmetry
	ape_sum=aperphot(image_data3,2.5,xo,yo)

	#print ape_sum
	#print len(ape_sum)
	#print ape_sum[0]
	#print ape_sum[57]

	x=np.ndarray( shape=(64,32,32), dtype=bool)
	xmask=np.ma.make_mask(x,copy=True, shrink=True, dtype=np.bool)
	xmask[:,:,:]= False
	xmask[:,13:18,13:18]=True
	masked= np.ma.masked_array(image_data3, mask=xmask)
	bgsum = np.zeros(64)
	for i in range (64):
		bgsum[i] = np.nansum(image_data1[i,:,:])/(32*32-25)
		#print i, bgsum[i]
	frameno=np.arange(0,64)
	fig, axes = plt.subplots(nrows=4, ncols=1, sharex=True)
	axes[0].plot(frameno,xo,color='k', mec ='r',marker='x', markevery=57)
	axes[0].set_ylabel('x0')


	axes[1].plot(frameno,yo,color='k' , mec ='r', marker='x', markevery=57)
	axes[1].set_ylabel('y0')

	axes[2].plot(frameno,bgsum,color='k', mec ='r', marker='x', markevery=57)
	axes[2].set_ylabel('b')

	axes[3].plot(frameno,ape_sum,color='k', mec ='r', marker='x', markevery=57)
	axes[3].set_ylabel('s')
	axes[3].set_xlabel('Frame number')
	fig.subplots_adjust(wspace=0)
	plt.savefig(str(ct)+'.png')
	ct+=1

'''
plt.xlabel('Frame number')
plt.ylabel('Background electron count per pixel')
for i in range (64):
	bgsum[i] = np.nansum(image_data1[i,:,:])/(32*32-25)
	#print i, bgsum[i]
frameno=np.arange(1,65)
plt.plot(frameno,bgsum)
xn+=1
plt.savefig(str(xn)+'.png')
'''
