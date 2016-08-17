#--------------------------------------------------------------
#Author: Lisa Dang & Sudarsan 
#Created: 2016-06-09 7:35 AM EST
#Last Modified: 
#Title: Comparing Gaussian 2d method and flux weighed mean method for centroiding
#--------------------------------------------------------------

import numpy as np
import matplotlib
from matplotlib import rc
import matplotlib.pyplot as plt
import matplotlib.patches
import time
from matplotlib.ticker import MaxNLocator
import os, sys
from astropy.io import fits
from astropy.stats import sigma_clip
from photutils import aperture_photometry
from photutils import CircularAperture
from numpy import std
import glob
import csv
import operator
import matplotlib.ticker as mtick
from photutils.datasets import make_4gaussians_image
from photutils.morphology import (centroid_com,centroid_1dg,centroid_2dg)
from time import time
from warnings import warn
from scipy.linalg.fblas import dgemm
#np.set_printoptions(threshold=np.nan)

def sigma_clipping(image_data):#,fname):
	#x=np.ndarray( shape=(32,32), dtype=bool)
	#xmask=np.ma.make_mask(x,copy=True, shrink=True, dtype=np.bool)
	#xmask[:,:]= True
	sig_clipped_data=sigma_clip(image_data, sigma=4, iters=2, cenfunc=np.ma.median,axis=0)
	for i in range (1,64):
		oldstar=image_data[i,12:19,12:19]
		#med= np.ma.median(sig_clipped_data,axis=0)
		newstar=sig_clipped_data[i,12:19,12:19]
		truth= newstar==oldstar
		if(truth.sum() <truth.size):
			temp=np.column_stack(np.where(truth!=True))
			temp=temp+12
			for iter in temp:
				sig_clipped_data[i,iter[0],iter[1]]=np.nan#med[iter[0],iter[1]]
	return sig_clipped_data

def bgsubtract(image_data):
	bgsubimg=image_data
	x=np.ndarray ( shape=(64,32,32), dtype=bool)
	xmask=np.ma.make_mask(x,copy=True, shrink=True, dtype=np.bool)
	xmask[:,:,:]= False
	xmask[:,13:18,13:18]=True
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

# To find centroid after a Gaussian 2d fit
# This function takes a lot of run time
def centroidg2d(image_data):
	cx=np.zeros(64)
	cy=np.zeros(64)
	for i in range(64):
		cx[i],cy[i]=centroid_2dg(image_data[i,13:18,13:18])
		cx[i]+=13
		cy[i]+=13
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

def plotandsave(xax,yax,xlabl,ylabl,fname):
	plt.clf()
	y_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
	ax = plt.gca()

	lb=np.min(np.concatenate([xax,yax]))
	ub=np.max(np.concatenate([xax,yax]))
	plt.xlim(lb,ub)
	plt.ylim(lb,ub)
	ax.set_aspect('equal', adjustable='box')

	print lb,ub
	ax.set(xlim=(lb,ub), ylim=(lb,ub))
	#from scipy.optimize import curve_fit
	A = np.vstack([xax, np.ones(len(xax))]).T
	m,c= np.linalg.lstsq(A,yax)[0]
	#m, c = curve_fit(f, xax, yax)[0]
	plt.plot(xax, yax, 'ro', label='Original data')
	#Pretty one line code to fit a straight line
	#lt.plot(xax, np.poly1d(np.polyfit(xax, yax, 1))(xax))
	plt.plot(xax, m*xax + c, 'b', label='Fitted line (m='+str.format('{0:.4f}', m)+')',c="0.3")
	plt.legend(loc=2)

	#ax.plot(xax,yax,'ro')
	ax.plot(ax.get_xlim(), ax.get_ylim(), ls="--", c=".3")
	ax.get_xaxis().get_major_formatter().set_useOffset(False)
	ax.get_yaxis().get_major_formatter().set_useOffset(False)
	plt.title(fname)
	plt.xlabel(xlabl,fontsize=16)
	plt.ylabel(ylabl,fontsize=16)
	plt.savefig('test/'+fname+'.png',bbox_inches='tight',dpi=500)

#note that the outerpath must be complete with terminating '/'
outerpath='/home/hema/Documents/mcgill/handy/aorkeys-20-selected_AORs/'
dirs=os.listdir(outerpath)
print dirs
counter=0
for direc in dirs :
	print direc
	if(counter==1):
		break
	normbg=[]
	normf=[]
	normx=[]
	normy=[]
	normGx=[]
	normGy=[]
	normpsfwx=[]
	normpsfwy=[]
	normnpp=[]
	path=outerpath+direc+'/ch2/bcd'
	print path
	#xn=1
	ct=0
	for filename in glob.glob(os.path.join(path, '*bcd.fits')):
		if(ct==10):
			break
		#print filename
		t = time()
		f=fits.open(filename,mode='readonly')
		image_data0=f[0].data
		# convert MJy/str to electron count
		convfact=f[0].header['GAIN']*f[0].header['EXPTIME']/f[0].header['FLUXCONV']
		image_data1=image_data0*convfact
		
		#sigma clip
		image_data2=sigma_clipping(image_data1)
		#bg subtract
		image_data3=bgsubtract(image_data2)
		#centroid
		xo, yo = centroid(image_data3)
		#apply gaussian 2d fit and find the centroid
		gx, gy = centroidg2d(image_data3)
		#aperture photmetry
		#ape_sum=aperphot(image_data3,2.5,xo,yo)
		
		#psfwx,psfwy=psfwidth(image_data3,xo,yo,7)
		#npp=noisepixparam(image_data3,7)
		ct+=1
		print time()-t
		plotandsave(xo,gx,r'$x_0$',r'$G x_0$','x0_vs_Gxo_'+filename[filename.find('I2_')+3: len(filename)-5])
		plotandsave(yo,gy,r'$y_0$',r'$G y_0$','y0_vs_Gyo_'+filename[filename.find('I2_')+3: len(filename)-5])
		
	print ct
	counter+=1