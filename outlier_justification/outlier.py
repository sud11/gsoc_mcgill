# To justify why 1st and 58 frames are outliers
# Do an aperture photometry and verify.

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
import operator
#np.set_printoptions(threshold=np.nan)

def sigma_clipping(image_data):#,fname):
	#global badframetable
	#global tossed
	sig_clipped_data=sigma_clip(image_data, sigma=4, iters=2, cenfunc=np.mean,axis=0)
	return sig_clipped_data

def bgsubtract(image_data):
	bgsubimg=image_data
	x=np.ndarray ( shape=(64,32,32), dtype=bool)
	xmask=np.ma.make_mask(x,copy=True, shrink=True, dtype=np.bool)
	xmask[:,:,:]= False
	xmask[:,14:18,14:18]=True
	#xmask[:,0:1,:]=True
	masked= np.ma.masked_array(bgsubimg, mask=xmask)
	n=0
	#Background subtraction for each frame
	while(n<64):
		bg_avg=np.ma.median(masked[n])

		# np.nanmedian functions weirdly for some reason.
		# below code to calculate median by ignoring nan values
		'''
		copy=[]
		for i in range (32):
			for j in range (32):
				if(masked[n][i][j]!= np.nan):
					copy.append( masked[n,i,j])
		bg_avg=np.median(copy)
		'''
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

def bgnormalize(image_data,normbg):
	x=np.ndarray( shape=(64,32,32), dtype=bool)
	xmask=np.ma.make_mask(x,copy=True, shrink=True, dtype=np.bool)
	xmask[:,:,:]= False
	xmask[:,13:18,13:18]=True
	masked= np.ma.masked_array(image_data, mask=xmask)
	#print masked[5][12:19,12:19]
	bgsum = np.zeros(64)
	# Replace for loop with one line code
	for i in range (64):
		bgsum[i] = np.nanmean(masked[i]) #np.ma.mean
	#background average for the datecube
	bgdcbavg= np.nanmean(bgsum)
	#Normalize
	bgsum=bgsum/bgdcbavg
	normbg.extend(bgsum)
	
	#print " normal ", bgsum[5]
	#bg_avg = np.mean(bgsum)
	#bgsum=bgsum/
def normstar(ape_sum,normf):
	starmean=np.nanmean(ape_sum)
	ape_sum=ape_sum/starmean
	normf.extend(ape_sum)
	#print min(enumerate(normf), key=operator.itemgetter(1))
	
def normxycent(xo,yo,normx,normy):
	xo=xo/np.nanmean(xo)
	yo=yo/np.nanmean(yo)
	normx.extend(xo)
	normy.extend(yo)
	

def reshapelists(normf,normbg,normx,normy,ct):
	normf=np.reshape(normf,(ct,64))
	normbg= np.reshape(normbg,(ct,64))
	normx=np.reshape(normx,(ct,64))
	normy=np.reshape(normy,(ct,64))
	return normf,normbg,normx,normy

def stackit(normf,normbg,normx,normy):
	normf=np.mean(normf,axis=0)
	normbg=np.mean(normbg, axis=0)
	normx=np.mean(normx,axis=0)
	normy=np.mean(normy,axis=0)
	return normf,normbg,normx,normy

	
def plotcurve(xax,f,b,X,Y):
	fig, axes = plt.subplots(nrows=4, ncols=1, sharex=True)
	plt.grid()
	plt.minorticks_on()
	axes[0].plot(xax,f,color='k', mec ='r', marker='x', markevery=57)
	axes[0].set_ylabel('F')

	axes[1].plot(xax,b,color='k', mec ='r', marker='x', markevery=57)
	axes[1].set_ylabel('b')

	axes[2].plot(xax,X,color='k', mec ='r',marker='x', markevery=57)
	axes[2].set_ylabel('x0')
	

	axes[3].plot(xax,Y,color='k' , mec ='r', marker='x', markevery=57)
	axes[3].set_ylabel('y0')
	axes[3].set_xlabel('Frame number')

	fig.subplots_adjust(wspace=0)
	plt.savefig('preprocsdB-sim'+'.png')

#Normalised and stacked
normbg=[]
normf=[]
normx=[]
normy=[]
path='/home/hema/Desktop/Link to mcgill/handy/datasample/simulated'
#xn=1
#for filename in glob.glob(os.path.join(path, '*bcd.fits')):
ct=0
for filename in glob.glob(os.path.join(path, '*bcd.fits')):
	f=fits.open(filename)
	if(ct==250):
		break
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

	bgnormalize(image_data1,normbg)
	normstar(ape_sum,normf)
	normxycent(xo,yo,normx,normy)
	ct+=1
print ct
normf,normbg,normx,normy=reshapelists(normf,normbg,normx,normy,ct)
normf,normbg,normx,normy=stackit(normf,normbg,normx,normy)
frameno=np.arange(0,64)
plotcurve(frameno,normf,normbg,normx,normy)