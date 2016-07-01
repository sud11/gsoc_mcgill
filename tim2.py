# CONVOLVING FLUXLIST and subtracting that from original flux list

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
import collections
#np.set_printoptions(threshold=np.nan)
ptoss = set([])
def sigma_clipping(image_data):#,fname):
	#x=np.ndarray( shape=(32,32), dtype=bool)
	#xmask=np.ma.make_mask(x,copy=True, shrink=True, dtype=np.bool)
	#xmask[:,:]= True
	sig_clipped_data=sigma_clip(image_data, sigma=4, iters=2, cenfunc=np.ma.median,axis=0)
	for i in range (0,len(image_data)):
		oldstar=image_data[i,12:19,12:19]
		newstar=sig_clipped_data[i,12:19,12:19]
		truth= (newstar==oldstar)
		if(truth.sum() <truth.size):
			sig_clipped_data[i,:,:]=np.nan
			ptoss.add(i)
	print toss
	return sig_clipped_data,ptoss

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
	print np.sum(image_data[20,13:18,13:18])
	#starbox=np.ndarray((len(image_data),5,5))
	#np.copyto(starbox,image_data[:,13:18,13:18])
	starbox = image_data[:, 13:18, 13:18]
	h,w = np.shape(starbox[0,:,:])
	x = np.arange(0,w)
	y = np.arange(0,h)
	X,Y = np.meshgrid(x,y)
	for i in range( len(image_data)):
		cx[i]=(np.sum(X*starbox[i,:,:])/np.sum(starbox[i,:,:]))+13
		cy[i]=(np.sum(Y*starbox[i,:,:])/np.sum(starbox[i,:,:]))+13
	return cx,cy

def aperphot(image_data,ape_radius,cx,cy,op):
	if(op==1):
		ape_sum=np.zeros(len(image_data))
		print "Radius:",ape_radius
		for i in range(len(image_data)):
			position=[cx[i],cy[i]]
			aperture=CircularAperture(position,r=ape_radius)
			phot_table=aperture_photometry(image_data[i,:,:],aperture)
			temp=phot_table['aperture_sum']
			ape_sum[i]=phot_table['aperture_sum']
	else:
		ape_sum=np.zeros(len(image_data))
		print "Radius:",ape_radius
		for i in range(len(image_data)):
			position=[cx[i],cy[i]]
			aperture=CircularAperture(position,r=ape_radius)
			phot_table=aperture_photometry(image_data[i,:,:],aperture,method='subpixel', subpixels=5)
			temp=phot_table['aperture_sum']
			ape_sum[i]=phot_table['aperture_sum']
	return ape_sum

def normstar(ape_sum):
	starmean=np.nanmean(ape_sum)
	ape_sum=ape_sum/starmean
	return ape_sum

def getflist(allstackd,op):
	'''
	allstackd=np.asarray(allstackd)
	temp=np.apply_along_axis(lambda m: convolve(m, Box1DKernel(50)), axis=0, arr=allstackd)
	'''
	image_data1,toss= sigma_clipping(allstackd)
	image_data2= bgsubtract(image_data1)
	cx, cy = centroid(image_data2)
	xax=[2.5]
	temp=0
	flist=aperphot(image_data2,2.5,cx,cy,op)
	return flist

# Flist after high pass filter
def highpassflist(allstackd,op):

	flist1=getflist(allstackd,op)
	print 'flist1',len(flist1)
	'''	
	# Binning version
	binstak=map(binning_data,flist1)
	binstak=np.asarray(binstak)
	print 'lengths',len(binstak),len(binstak[0])
	flist1=binstak
	'''
	flist2=convolve(flist1, Box1DKernel(50),boundary='extend')
	flist1=np.asarray(flist1)
	flist2=np.asarray(flist2)
	
	smooth=flist1-flist2
	
	'''
	sflist = map(normstar,smooth)
	sflist = map( lambda t: t-1, sflist)
	panel=map(lambda t: np.sqrt(np.nanmean(np.square(t))), sflist)
	'''
	return flist1,flist2,smooth

def binning_data(data):
	#reshaped_data = data.reshape((h,w))
	#binned_data=np.ma.median(data, axis=0)
	binnd=[]
	for i in range (0, len(data)/64):
		temp=np.ma.median(data[i*64:(i+1)*64])
		binnd.append(temp)
	return binnd
	#binned_data_std=np.std(reshaped_data, axis=1)
	return binned_data

def binning_time(data, h, w):
	reshaped_data = data.reshape((h,w))
	binned_data=np.ma.median(reshaped_data, axis=1)
	binned_data_std=np.std(reshaped_data, axis=1)
	return binned_data


def plot3panels(xax,p1,p2,p3,p4,toss,ct):
	toss=list(toss)
	p1=np.asarray(p1)
	p2=np.asarray(p2)
	p3=np.asarray(p3)
	p4=np.asarray(p4)
	xax=np.asarray(xax)
	fig, axes = plt.subplots(nrows=4, ncols=1, sharex=True)
	plt.minorticks_on()
	fig.subplots_adjust(hspace = 0.001)
	plt.rc('font', family='serif',serif='Times')
	y_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
	
	axes[0].plot(xax,p1,'-',c='k',mec='b',fillstyle='none')
	axes[0].set_ylabel(r'$original$',fontsize=13)
	axes[0].yaxis.set_major_formatter(y_formatter)
	axes[0].xaxis.set_major_formatter(y_formatter)
	axes[0].yaxis.set_major_locator(MaxNLocator(prune='both',nbins=5))

	axes[1].plot(xax,p2,'-',c='k',mec='b',marker='x', markevery=toss,fillstyle='none')
	axes[1].set_ylabel(r'$smoothed$',fontsize=13)
	axes[1].yaxis.set_major_formatter(y_formatter)
	axes[1].xaxis.set_major_formatter(y_formatter)
	axes[1].yaxis.set_major_locator(MaxNLocator(prune='both',nbins=5))

	axes[2].plot(xax,p3,'-',c='k',mec='b',marker='x', markevery=toss,fillstyle='none')
	axes[2].set_ylabel(r'$flattened$',fontsize=13)
	axes[2].yaxis.set_major_formatter(y_formatter)
	axes[2].xaxis.set_major_formatter(y_formatter)
	axes[2].yaxis.set_major_locator(MaxNLocator(prune='both',nbins=5))

	axes[3].plot(xax,p4,'-',c='k',mec='b',fillstyle='none',label='Soft')
	axes[3].set_ylabel(r'$\sigma$ $clipped$',fontsize=13)
	axes[3].set_xlabel(r'$time$ $(hours)$',fontsize=13)
	axes[3].yaxis.set_major_formatter(y_formatter)
	axes[3].xaxis.set_major_formatter(y_formatter)
	axes[3].yaxis.set_major_locator(MaxNLocator(prune='both',nbins=5))
	#axes[2].legend(numpoints=1)
	plt.savefig(str(ct)+'_4_lightcurve_tim2.png',bbox_inches='tight',dpi=200)

def plotcurve(xax,f1,f2,ct):
	fig, axes = plt.subplots(nrows=4, ncols=1, sharex=True)
	plt.minorticks_on()
	fig.subplots_adjust(hspace = 0.001)
	plt.rc('font', family='serif',serif='Times')
	y_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
	
	axes[0].plot(xax,f1[0],'D-',c='k',mec='b',fillstyle='none')
	axes[0].plot(xax,f2[0],'o-',c='g',mec='k',fillstyle='none')
	axes[0].set_ylabel(r'$raw$ $RMS$',fontsize=13)
	axes[0].yaxis.set_major_formatter(y_formatter)
	axes[0].yaxis.set_major_locator(MaxNLocator(prune='both',nbins=5))

	axes[1].plot(xax,f1[1],'D-',c='k',mec='b',fillstyle='none')
	axes[1].plot(xax,f2[1],'o-',c='g',mec='k',fillstyle='none')
	axes[1].set_ylabel(r'$frames$ $RMS$',fontsize=13)
	axes[1].yaxis.set_major_formatter(y_formatter)
	axes[1].yaxis.set_major_locator(MaxNLocator(prune='both',nbins=5))

	axes[2].plot(xax,f1[2],'D-',c='k',mec='b',fillstyle='none')
	axes[2].plot(xax,f2[2],'o-',c='g',mec='k',fillstyle='none')
	axes[2].set_ylabel(r'$\sigma-clipped$',fontsize=13)
	axes[2].yaxis.set_major_formatter(y_formatter)
	axes[2].yaxis.set_major_locator(MaxNLocator(prune='both',nbins=5))

	axes[3].plot(xax,f1[3],'D-',c='k',mec='b',fillstyle='none',label='Hard')
	axes[3].plot(xax,f2[3],'o-',c='g',mec='k',fillstyle='none',label='Soft')
	axes[3].set_ylabel(r'$\sigma$ $clipped$ $RMS$',fontsize=13)
	axes[3].set_xlabel(r'$aperture$ $(pixels)$',fontsize=13)
	axes[3].yaxis.set_major_formatter(y_formatter)
	axes[3].yaxis.set_major_locator(MaxNLocator(prune='both',nbins=5))
	axes[3].legend(numpoints=1)
	plt.savefig('paneltest/'+str(ct)+'updchanges.png',bbox_inches='tight',dpi=200)

def masterfunc1(op):
	xax=[2.5]
	global allstackd1
	allstackd1=np.asarray(allstackd1)
	print 'allstackd1',len(allstackd1)
	panel1,panel2,smooth1=highpassflist(allstackd1,op)


	'''
	panel1=map(normstar,flist1)
	panel1= map( lambda t: t-1, panel1)
	panel1=map(lambda t: np.sqrt(np.nanmean(np.square(t))), panel1)
	'''
	smooth1=np.asarray(smooth1)
	beforenan= np.argwhere(np.isnan(smooth1))
	beforenan=beforenan.flatten()
	panel3= sigma_clip(smooth1,sigma=5,cenfunc=np.nanmedian)
	p=np.ma.fix_invalid(panel3.data,panel3.mask,fill_value=0)
	tosstemp= [i for i, x in enumerate(panel3.mask) if x]
	panel3=np.asarray(panel3)
	for x in tosstemp:
		panel3[x]=np.nan
	afternan = np.ma.asarray(panel3)
	afternan = np.argwhere(np.isnan(afternan))
	afternan = afternan.flatten()
	'''
	smooth1=np.asarray(smooth1)
	beforenan= np.argwhere(np.isnan(smooth1))
	beforenan=beforenan.flatten()
	panel3= sigma_clip(smooth1,sigma=5)
	afternan = np.asarray(panel3)
	afternan = np.argwhere(np.isnan(afternan))
	afternan = afternan.flatten()
	'''

	toss = set(tosstemp) - set(beforenan)
	#beforenan= map(lambda t: beforenan[t], np.arange(0,len(beforenan)))
	#toss=np.asarray(afternan)-np.asarray(beforenan)
	#print toss
	
	'''
	panel3=map(normstar,panel3)

	panel3= map( lambda t: t-1, panel3)
	panel3=map(lambda t: np.sqrt(np.nanmean(np.square(t))), panel3)
	'''
	panel2=np.asarray(panel2)
	panels=[]
	panels.append(panel1)
	panels.append(panel2)
	panels.append(smooth1)
	panels.append(panel3)
	return panels,toss

def getfnames():
	outerpath='/home/hema/Documents/mcgill/handy/aorkeys-20-selected_AORs/'
	dirs=os.listdir(outerpath)
	print dirs
	counter=0
	ct=0
	nametime={}
	for direc in dirs :
		if(counter==1):
			break
		path=outerpath+direc+'/ch2/bcd'
		ct=0
		for filename in glob.glob(os.path.join(path, '*bcd.fits')):
			f=fits.open(filename,mode='readonly')
			nametime[filename] = f[0].header['AINTBEG']
			ct+=1
		counter+=1
	#sorting nametime based on time
	nametime = sorted(nametime.items(), key=operator.itemgetter(1))
	#fnames is the list of paths to filenames
	fnames= [x[0] for x in nametime]
	ftimes= [x[1] for x in nametime]
	return fnames,ftimes

allstackd1=[]
fnames,ftimes=getfnames()
ct=0
t=[]
for filename in fnames:
	#print filename
	if(ct==10):
		break
	f=fits.open(filename,mode='readonly')
	print filename
	image_data0=f[0].data
	t.extend( np.linspace(f[0].header['AINTBEG'],f[0].header['ATIMEEND'],64) )
	print f[0].header['AINTBEG'],f[0].header['ATIMEEND']
	# convert MJy/str to electron 
	convfact=f[0].header['GAIN']*f[0].header['EXPTIME']/f[0].header['FLUXCONV']
	image_data1=image_data0*convfact
	allstackd1.extend(image_data1)
	ct+=1

#convert time from sec to hours
t=np.asarray(t)
sectohrs=1.0/3600.0
time=t*sectohrs

print time[0:20]

panels,toss=masterfunc1(1)
'''
print 'panel1', panels[0][30:40]
print 'panel2', panels[1][30:40]
print 'panel3', panels[2][30:40]
print 'panel4', panels[3][30:40]

#panels=map(normstar,panels)
#lot3panels(time,panels[0],panels[1],panels[2],panels[3],toss,ct)

#plot3panels(time[30:-30],panels[0][30:-30],panels[1][30:-30],panels[2][30:-30],ct)
#xax=np.arange(1.5,4.75,step=0.25)
#plotcurve(xax,flist,ct)

'''


