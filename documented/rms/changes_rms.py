#--------------------------------------------------------------
#Author: Lisa Dang & Sudarsan 
#Created: 2016-06-28 4:13 PM EST
#Last Modified: 
#Title: Identifying optimal aperture radius for aperture photometry
#--------------------------------------------------------------

'''
Changes from multip_rms.py -
1. For the binned plots, I've median binned the flux list ( earlier I was median binning raw image data)
2. High pass filter the lightcurve (i.e flux list. Earlier, I was passing the raw image data through the high pass filter)
'''
#Test data - 3.6 AOR ['r46483456', 'r46482944', 'r46483200', 'r46482688']
#            4.5 AOR ['r46470912', 'r46469120','r46468608','r46470400']
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
from photutils import CircularAnnulus
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
from functools import partial

from multiprocessing import Pool
np.set_printoptions(threshold=np.nan)
# for np.1dunion as it takes only 2 args
from functools import reduce
import collections
def sigma_clipping(image_data):#,fname):
	#x=np.ndarray( shape=(32,32), dtype=bool)
	#xmask=np.ma.make_mask(x,copy=True, shrink=True, dtype=np.bool)
	#xmask[:,:]= True
	toss = set([])
	sig_clipped_data=sigma_clip(image_data, sigma=4, iters=2, cenfunc=np.ma.median,axis=0)
	for i in range (0,len(image_data)):
		oldstar=image_data[i,12:19,12:19]
		newstar=sig_clipped_data[i,12:19,12:19]
		truth= (newstar==oldstar)
		if(truth.sum() <truth.size):
			sig_clipped_data[i,:,:]=np.nan
			toss.add(i)
	return sig_clipped_data,toss

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
	xax=np.arange(1.5,4.75,step=0.25)
	temp=0
	flist=map(lambda t: aperphot(image_data2,t,cx,cy,op), xax)
	return flist

# Flist after high pass filter
def highpassflist(allstackd,op):

	flist1=getflist(allstackd,op)
	'''	
	# Binning version
	binstak=map(binning_data,flist1)
	binstak=np.asarray(binstak)
	print 'lengths',len(binstak),len(binstak[0])
	flist1=binstak
	'''
	flist2=np.apply_along_axis(lambda m: convolve(m, Box1DKernel(50)),axis=0,arr=flist1)
	flist1=np.asarray(flist1)
	flist2=np.asarray(flist2)
	smooth=flist1-flist2
	sflist = map(normstar,smooth)
	sflist = map( lambda t: t-1, sflist)
	panel=map(lambda t: np.sqrt(np.nanmean(np.square(t))), sflist)
	return panel,smooth

def tossoutframe(allstackd,op):
	#5 sigmaclip of x0,y0,F
	image_data1,toss1= sigma_clipping(allstackd)
	image_data2= bgsubtract(image_data1)
	cx, cy = centroid(image_data2)
	cx= sigma_clip(cx,sigma=5)
	cy= sigma_clip(cy,sigma=5)
	#bad_onlyx = cx.data[cx.mask]
	#bad_onlyy = cy.data[cy.mask]
	t1=np.argwhere(np.isnan(cx))
	t2=np.argwhere(np.isnan(cy))
	print 't1',t1
	print 't2',t2
	xax=np.arange(1.5,4.75,step=0.25)
	flist=map(lambda t:aperphot(image_data2,t,cx,cy,op), xax)

	beforenan= np.argwhere(np.isnan(flist))
	beforenan =  collections.Counter([x for (x,y) in beforenan])
	print 'beforenan',beforenan

	flist=map(lambda t:sigma_clip(t,sigma=5),flist)
	flist=np.asarray(flist)
	
	nanlist3= np.argwhere(np.isnan(flist))
	afternan = collections.Counter([x for (x,y) in nanlist3])
	
	print 'afternan',afternan
	temp1 = beforenan-afternan

	toss=np.empty(len(flist)); toss.fill(0)
	toss= map( lambda t: temp1[t], np.arange(0,13) )
	
	#tossed=reduce(np.union1d,(nanlist1,nanlist2,nanlist3))
	'''
	flist = map(normstar,flist)
	flist = map( lambda t: t-1, flist)
	'''
	return flist,toss

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
def plot3panels(xax,p1,p2,p3,p4,ct):
	fig, axes = plt.subplots(nrows=4, ncols=1, sharex=True)
	plt.minorticks_on()
	fig.subplots_adjust(hspace = 0.001)
	plt.rc('font', family='serif',serif='Times')
	y_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
	
	axes[0].plot(xax,p1,'D-',c='k',mec='b',fillstyle='none')
	axes[0].set_ylabel(r'$original$ $RMS$',fontsize=13)
	axes[0].yaxis.set_major_formatter(y_formatter)
	axes[0].yaxis.set_major_locator(MaxNLocator(prune='both',nbins=5))

	axes[1].plot(xax,p2,'D-',c='k',mec='b',fillstyle='none')
	axes[1].set_ylabel(r'$flattened$ $RMS$',fontsize=13)
	axes[1].yaxis.set_major_formatter(y_formatter)
	axes[1].yaxis.set_major_locator(MaxNLocator(prune='both',nbins=5))

	axes[2].plot(xax,p4,'D-',c='k',mec='b',fillstyle='none')
	axes[2].set_ylabel(r'$\sigma-clipped$',fontsize=13)
	axes[2].yaxis.set_major_formatter(y_formatter)
	axes[2].yaxis.set_major_locator(MaxNLocator(prune='both',nbins=5))
	
	axes[3].plot(xax,p3,'D-',c='k',mec='b',fillstyle='none',label='Soft')
	axes[3].set_ylabel(r'$\sigma$ $clipped$ $RMS$',fontsize=13)
	axes[3].set_xlabel(r'$aperture$ $(pixels)$',fontsize=13)
	axes[3].yaxis.set_major_formatter(y_formatter)
	axes[3].yaxis.set_major_locator(MaxNLocator(prune='both',nbins=5))
	axes[3].legend(numpoints=1)
	plt.savefig(str(ct)+'smoothlightS.png',bbox_inches='tight',dpi=200)

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
	plt.savefig(str(ct)+'updchanges.png',bbox_inches='tight',dpi=200)

def masterfunc1(op):
	xax=np.arange(1.5,4.75,step=0.25)
	global allstackd1
	global allstackd2
	allstackd1=np.asarray(allstackd1)
	allstackd2=np.asarray(allstackd2)
	
	panel1,smooth1=highpassflist(allstackd1,op)
	panel2,smooth2=highpassflist(allstackd2,op)

	beforenan= np.argwhere(np.isnan(smooth2))
	beforenan =  collections.Counter([x for (x,y) in beforenan])

	panel4= map(lambda t: sigma_clip(smooth2[t],sigma=5),np.arange(0,len(xax)))

	noisef1=np.ma.asarray(panel4)
	#panel3=np.ma.asarray(panel3)
	
	afternan = map( np.ma.count_masked,panel4)
	beforenan= map(lambda t: beforenan[t], np.arange(0,len(beforenan)))
	toss=np.asarray(afternan)-np.asarray(beforenan)

	panel4=map(normstar,panel4)
	panel4= map( lambda t: t-1, panel4)
	panel4=map(lambda t: np.sqrt(np.nanmean(np.square(t))), panel4)

	panel3=toss/(ct*64.0)
	panels=[]
	panels.append(panel1)
	panels.append(panel2)
	panels.append(panel3)
	panels.append(panel4)
	return panels



outerpath='/home/hema/Documents/mcgill/handy/aorkeys-20-selected_AORs/'
dirs=os.listdir(outerpath)
print dirs
allstackd1=[]
allstackd2=[]
counter=0
ct=0
for direc in ['r46470400'] :
	print direc
	if(counter ==1):
		break
	path=outerpath+direc+'/ch2/bcd'
	for filename in glob.glob(os.path.join(path, '*bcd.fits')):
		f=fits.open(filename,mode='readonly')
		image_data0=f[0].data
		# convert MJy/str to electron count
		convfact=f[0].header['GAIN']*f[0].header['EXPTIME']/f[0].header['FLUXCONV']
		image_data1=image_data0*convfact
		allstackd1.extend(image_data1)
		allstackd2.extend(image_data1[1:])
		
		ct+=1
	print 'ct',ct
	counter+=1

t=time()
pool1= Pool(processes=4)
f1,f2=pool1.map(masterfunc1,[1,2])

xax=np.arange(1.5,4.75,step=0.25)
plotcurve(xax,f1,f2,ct)
#plot3panels(xax,panel1,panel2,panel3,panel4,ct)


'''
binsmoo=[]
for i in range(0, len(smooth1)/64):
	binsmoo.append(binning_data(smooth1[i*64:(i+1)*64]))
'''
#binsmoo=np.asarray(binsmoo)
'''
panel2=getflist(binsmoo,1)


noisef1= map(lambda t: panel1[t]-panel2[t], np.arange(0,len(xax)))

beforenan= np.argwhere(np.isnan(noisef1))
beforenan =  collections.Counter([x for (x,y) in beforenan])

panel3= map(lambda t: sigma_clip(noisef1[t],sigma=5),np.arange(0,len(xax)))

noisef1=np.ma.asarray(noisef1)
#panel3=np.ma.asarray(panel3)
print 'count',np.ma.count_masked(panel3[0])

afternan = map( np.ma.count_masked,panel3)
beforenan= map(lambda t: beforenan[t], np.arange(0,len(beforenan)))
toss=np.asarray(afternan)-np.asarray(beforenan)

#toss= map( lambda t: temp1[t], np.arange(0,13))
toss=np.asarray(toss)
panel1=map(lambda t: np.sqrt(np.nanmean(np.square(t))), panel1)
panel2=map(lambda t: np.sqrt(np.nanmean(np.square(t))), noisef1)
panel3=map(lambda t: np.sqrt(np.nanmean(np.square(t))), panel3)
panel4=toss/float(ct)
plot3panels(xax,panel1,panel2,panel3,panel4,ct)
#plotcurve(xax,f1,f2,ct)
'''