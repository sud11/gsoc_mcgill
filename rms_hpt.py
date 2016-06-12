# Pass the data through a high pass filter.
# Plot for F RMS value for different aperture sizes
# For a high pass filter, do a boxcar smooth and subtract the smoothed data from raw data

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
from scipy.linalg.fblas import dgemm
from astropy.convolution import convolve, Box1DKernel

#np.set_printoptions(threshold=np.nan)

def sigma_clipping(image_data):#,fname):
	#x=np.ndarray( shape=(32,32), dtype=bool)
	#xmask=np.ma.make_mask(x,copy=True, shrink=True, dtype=np.bool)
	#xmask[:,:]= True
	sig_clipped_data=sigma_clip(image_data, sigma=4, iters=2, cenfunc=np.ma.median,axis=0)
	for i in range (1,len(image_data)):
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
	x=np.ndarray ( shape=(len(image_data),32,32), dtype=np.bool)
	xmask=np.ma.make_mask(x,copy=True, shrink=True,dtype=np.bool)
	xmask[:,:,:] = False
	xmask[:,14:18,14:18]=True
	#xmask[:,0:1,:]=True
	masked= np.ma.masked_array(bgsubimg, mask=xmask)
	n=0
	#Background subtraction for each frame
	while(n<len(image_data)):
		bg_avg=np.ma.median(masked[n])
		bgsubimg[n]= bgsubimg[n,:,:] - bg_avg
		n+=1
	return bgsubimg

def centroid(image_data):
	# Refer: Intra-Pixel Gain Variations and High-Precision Photometry with the Infrared Array Camera (IRAC)
	cx=np.zeros(len(image_data))
	cy=np.zeros(len(image_data))
	starbox = image_data[:, 13:18, 13:18]
	h,w = np.shape(starbox[0,:,:])
	x = np.arange(0,w)
	y = np.arange(0,h)
	X,Y = np.meshgrid(x,y)
	for i in range( len(image_data)):
		cx[i]=(np.sum(X*starbox[i,:,:])/np.sum(starbox[i,:,:]))+13
		cy[i]=(np.sum(Y*starbox[i,:,:])/np.sum(starbox[i,:,:]))+13
	return cx,cy

# To find centroid after a Gaussian 2d fit
# This function takes a lot of run time

def aperphot(image_data,ape_radius,cx,cy):
	ape_sum=np.zeros(len(image_data))
	print "Radius:",ape_radius
	for i in range(len(image_data)):
		position=[cx[i],cy[i]]
		aperture=CircularAperture(position,r=ape_radius)
		phot_table=aperture_photometry(image_data[i,:,:],aperture)
		temp=phot_table['aperture_sum']
		ape_sum[i]=phot_table['aperture_sum']
	return ape_sum

def normstar(ape_sum):
	starmean=np.nanmean(ape_sum)
	ape_sum=ape_sum/starmean
	return ape_sum

def getflist(allstackd):
	#Converting to numpy ndarray
	allstackd=np.asarray(allstackd)
	#print allstackd[4,13:18,13:18]
	#Boxcar Smooth
	temp=np.apply_along_axis(lambda m: convolve(m, Box1DKernel(50)), axis=0, arr=allstackd)
	#print temp[4,13:18,13:18]
	#Smoothening
	allstackd-=temp
	#allstackd1=np.subtract(allstackd,temp)

	image_data1= sigma_clipping(temp)
	image_data2= bgsubtract(image_data1)
	cx, cy = centroid(image_data2)

	#apertures
	xax=np.arange(1.5,4.75,step=0.25)
	flist=map(lambda t: aperphot(image_data2,t,cx,cy), xax)
	#print 'aperturephot',time()-clo
	flist = map(normstar,flist)
	flist = map( lambda t: t-1, flist)
	#print 'after norm',flist[0][10:30]

	'''
	print 'dimensions of flist',len(flist),len(flist[0])
	print 'step1',flist[0]
	tempo=np.square(flist[0])
	print 'square',tempo[10:30]
	print 'nanmean',np.nanmean(np.square(flist[0]))
	print 'rms',np.sqrt(np.nanmean(np.square(flist[0])))
	'''
	flist=map(lambda t: np.sqrt(np.nanmean(np.square(t))), flist)
	print 'after:'
	print 'dimensions of flist',len(flist)
	return flist

def tossoutframe(allstackd):
	#5 sigmaclip of x0,y0,F
	allstackd=np.asarray(allstackd)
	temp=np.apply_along_axis(lambda m: convolve(m, Box1DKernel(50)), axis=0, arr=allstackd)
	#print temp[4,13:18,13:18]
	#Smoothening
	allstackd-=temp
	image_data1= sigma_clipping(temp)
	image_data2= bgsubtract(image_data1)
	cx, cy = centroid(image_data2)
	cx=sigma_clip(cx,sigma=5)
	cy=sigma_clip(cy,sigma=5)
	xax=np.arange(1.5,4.75,step=0.25)
	flist=map(lambda t: sigma_clip(aperphot(image_data2,t,cx,cy),sigma=5), xax)
	flist=np.asarray(flist)
	flist = map(normstar,flist)
	flist = map( lambda t: t-1, flist)
	flist=map(lambda t: np.sqrt(np.nanmean(np.square(t))), flist)
	return flist


def plotandsave(xax,yax,xlabl,ylabl,fname):
	plt.clf()
	# y_formatter for pretty printing scale values
	y_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
	ax = plt.gca()
	plt.figure(figsize=(5,1))
	plt.tick_params(axis='x', labelsize=7)
	plt.tick_params(axis='y', labelsize=7)
	plt.plot(xax, yax, 'xb-')
	#plt.xlim([min(xax) - 0.25, max(xax) + 0.25])
	plt.xlabel(xlabl,fontsize=8)
	plt.ylabel(ylabl,fontsize=8)
	plt.savefig('paneltest/'+fname+'.png',bbox_inches='tight',dpi=200)

def plotcurve(xax,flist,ct):
	fig, axes = plt.subplots(nrows=3, ncols=1, sharex=True)
	plt.minorticks_on()
	fig.subplots_adjust(hspace = 0.001)
	plt.rc('font', family='serif')
	#fig.subplots_adjust(.15,.15,.9,.9,0,0)
	y_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
	
	axes[0].plot(xax,flist[0],'xb-')
	axes[0].set_ylabel(r'$raw RMS$',fontsize=16)
	axes[0].yaxis.set_major_formatter(y_formatter)
	axes[0].yaxis.set_major_locator(MaxNLocator(prune='both',nbins=5))

	axes[1].plot(xax,flist[1],'xb-')
	axes[1].set_ylabel(r'$frames RMS$',fontsize=16)
	axes[1].yaxis.set_major_formatter(y_formatter)
	axes[1].yaxis.set_major_locator(MaxNLocator(prune='both',nbins=5))

	axes[2].plot(xax,flist[2],'xb-')
	axes[2].set_ylabel(r'$\sigma-clipped$',fontsize=16)
	axes[2].set_xlabel(r'$Frame$ $number$',fontsize=16)
	axes[2].yaxis.set_major_formatter(y_formatter)
	axes[2].yaxis.set_major_locator(MaxNLocator(prune='both',nbins=5))
	plt.savefig('paneltest/rms3panesl'+str(ct)+'.png',bbox_inches='tight',dpi=200)

#note that the outerpath must be complete with terminating '/'
outerpath='/home/hema/Documents/mcgill/handy/aorkeys-20-selected_AORs/'
dirs=os.listdir(outerpath)
print dirs
allstackd1=[]
allstackd2=[]
flist=[]
counter=0
ct=0
for direc in dirs :
	print direc
	if(counter==1):
		break
	path=outerpath+direc+'/ch2/bcd'
	for filename in glob.glob(os.path.join(path, '*bcd.fits')):
		#print filename
		if(ct==5):
			break
		f=fits.open(filename,mode='readonly')
		image_data0=f[0].data
		# convert MJy/str to electron count
		convfact=f[0].header['GAIN']*f[0].header['EXPTIME']/f[0].header['FLUXCONV']
		image_data1=image_data0*convfact
		#smoothed_signal = convolve(noisy_signal, Box1DKernel(11))
		allstackd1.extend(image_data1)
		allstackd2.extend(image_data1[1:])
		#plt.imshow(smoothed_signal)
		#plt.show()
		ct+=1
		#plotandsave(xo,gx,r'$x_0$',r'$G x_0$','x0_vs_Gxo_'+filename[filename.find('I2_')+3: len(filename)-5])
		#plotandsave(yo,gy,r'$y_0$',r'$G y_0$','y0_vs_Gyo_'+filename[filename.find('I2_')+3: len(filename)-5])
	print 'ct',ct
	counter+=1


flist1= getflist(allstackd1)
flist2=getflist(allstackd2)
flist3=tossoutframe(allstackd2)

#flist4 = map(normstar,flist4)
#flist4 = map( lambda t: t-1, flist4)
#flist4=map(lambda t: np.sqrt(np.nanmean(np.square(t))), flist4)

flist.append(flist1)
flist.append(flist2)
flist.append(flist3)

xax=np.arange(1.5,4.75,step=0.25)
plotcurve(xax,flist,ct)
#plotandsave(xax,flist1,'Aperture radius (pixels) ','raw RMS','rmstest7_1aor_10f_p1')
#plotandsave(xax,flist2,'Aperture radius (pixels) ','frames RMS','rmstest7_1aor_10f_p2')
#plotandsave(xax,flist3,'Aperture radius (pixels) ','sig clipped','rmstest7_1aor_10f_p3')
#plotandsave(xax,flist4,'Aperture radius (pixels) ','panel4 ','rmstest7_1aor_10f_p4')




'''
F = aperphot(image_data2, 1.5,cx,cy)
print np.sqrt(np.nanmean(F**2))
F = aperphot(image_data2, 2,cx,cy)
print np.sqrt(np.nanmean(F**2))
F = aperphot(image_data2, 2.5,cx,cy)
print np.sqrt(np.nanmean(F**2))
F = aperphot(image_data2, 3,cx,cy)
print np.sqrt(np.nanmean(F**2))
F = aperphot(image_data2, 3.5,cx,cy)
print np.sqrt(np.nanmean(F**2))
F = aperphot(image_data2, 4,cx,cy)
print np.sqrt(np.nanmean(F**2))
F = aperphot(image_data2, 4.5,cx,cy)
print np.sqrt(np.nanmean(F**2))
'''