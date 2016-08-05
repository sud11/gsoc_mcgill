# Connected component technique to separate noise
# Code to identify the best level to separate out just the noise
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
import photutils
from photutils.morphology import (centroid_com,centroid_1dg,centroid_2dg)

import astropy.units as u

def sigma_clipping(image_data):
	ptoss=set([])
	image_data=np.asarray(image_data)
	sig_clipped_data=sigma_clip(image_data, sigma=4, iters=2, cenfunc=np.ma.median,axis=0)
	for i in range (0,len(image_data)):
		oldstar=image_data[i,12:19,12:19]
		newstar=sig_clipped_data[i,12:19,12:19]
		truth= (newstar==oldstar)
		if(truth.sum() <truth.size):
			sig_clipped_data[i,:,:]=np.nan
			ptoss.add(i)
	return sig_clipped_data
	
def bgsubtract(image_data):
	bgsubimage_data=np.zeros((64,32,32))
	x=np.ndarray( shape=(64,32,32), dtype=bool)
	xmask=np.ma.make_mask(x,copy=True, shrink=True, dtype=np.bool)
	xmask[:,:,:]= False
	xmask[:,13:18,13:18]=True #removing starbox
	xmask[:,0:1,:]=True #Removing top defective rows 
	masked= np.ma.masked_array(image_data, mask=xmask)
	bg_err = np.zeros(64)	
	n=0
	while(n<64):
		bg_avg=np.ma.median(masked[n,:,:])
		bg_err[n] = np.ma.var(masked[n,:,:])
		bgsubimage_data[n]=image_data[n,:,:] - bg_avg
		n+=1
	return bgsubimage_data

def photometry(image_data, ape_radius, cx, cy):
	ape_sum=np.zeros(64)
	ape_sum_err=np.zeros(64)
	for i in range(64):
		position=[cx[i], cy[i]]
		aperture=CircularAperture(position, r=ape_radius)
		#http://photutils.readthedocs.io/en/latest/photutils/aperture.html#error-estimation
		#aperture_sum_error can be accessed only after passing a error parameter
		#data_error = bg_err[i]
		phot_table=aperture_photometry(image_data[i,:,:],aperture, pixelwise_error=False)
		ape_sum[i]=phot_table['aperture_sum']
		
	return ape_sum

# To find centroid after a Gaussian 2d fit
# This function takes a lot of run time
def centroidg2d(image_data):
	cx=np.zeros(64)
	cy=np.zeros(64)
	for i in range(64):
		cx[i],cy[i]=centroid_2dg(image_data3[i,13:18,13:18])
		cx[i]+=13
		cy[i]+=13
	return cx,cy

outerpath='/home/hema/Documents/mcgill/handy/corot2b/'
dirs=os.listdir(outerpath)
print dirs	
counter=0
ct=0
direc= dirs[0]
path=outerpath+direc+'/ch2/bcd'
for filename in glob.glob(os.path.join(path, '*bcd.fits')):
	if(ct==1):
		break
	f=fits.open(filename,mode='readonly')
	image_data0=f[0].data
	# As the image data is given in MJy/sr, but in the catalog the units are mJy, we have to convert the units
	factor = (u.MJy / u.sr * (1.2 * u.arcsec) ** 2 / u.pixel).to(u.mJy / u.pixel)
	image_data1= image_data0#*factor.value
	maskd = np.isfinite(image_data1)
	image_data2=sigma_clipping(image_data1)
	print 'pp',len(image_data2)
	#bg subtract
	image_data3=bgsubtract(image_data2)
	#get centroid
	xo_tmp, yo_tmp = centroidg2d(image_data3)

	flist=photometry(image_data3,2.5, xo_tmp, yo_tmp)
	print 'center weighed',flist[10]
	t1=flist[10]
	coords= zip(xo_tmp,yo_tmp)
	from photutils.psf import psf_photometry
	maskd=np.isfinite(image_data3)
	mx = np.ma.masked_array(image_data3, mask=maskd)
	
	from skimage import measure
	L = measure.label(image_data3[10])
	print "Number of components:", np.max(L)
	from matplotlib.backends.backend_pdf import PdfPages
	import matplotlib.gridspec as gridspec
	pdf =PdfPages( 'foo.pdf' )
	fig = plt.figure(figsize = (9,9))
	gs = gridspec.GridSpec(9, 9)#,wspace=0.025, hspace=0.05)
	i=0
	
	for k in np.arange(.5,10,.125):
		contours = measure.find_contours(mx[10], k)
		ax = plt.subplot(gs[i])
		ax.set_title(str(k),fontsize=5)
		ax.imshow(image_data3[10], interpolation='nearest', cmap=plt.cm.gray)
		ax.set_aspect('auto')
		for n, contour in enumerate(contours):
			ax.plot(contour[:, 1], contour[:, 0], linewidth=.1)
		ax.axis('image')
		ax.set_aspect('equal')
		ax.set_xticks([])
		ax.set_yticks([])
		i+=1
		#plt.figure(num=None, figsize=(.125, .125), dpi=80)
	plt.tight_layout()
	pdf.savefig(dpi=500)
	pdf.close()


	
	#fluxes_gaussian = photutils.psf.psf_photometry(mx[10],[coords[10]], psf_gaussian,fitshape=None)
	#print 'gaussian',fluxes_gaussian
	#t2=fluxes_gaussian
	#print 'fac',t2-t1
	ct+=1


'''
# photometry using a GaussianPSF https://github.com/astropy/photutils-datasets/blob/master/notebooks/PSFPhotometrySpitzer.ipynb
from photutils.psf import GaussianPSF
psf_gaussian = GaussianPSF(1)
fluxes_gaussian = psf_photometry(data, coords, psf_gaussian) #Fluxes are in mJy
data = fits.getdata('../data/spitzer_example_image.fits')
'''


