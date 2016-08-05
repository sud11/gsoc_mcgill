# To justify why 1st and 58 frames are outliers
# Do an aperture photometry and verify.

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
	xmask[:,0:1,:]=True
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
		cx[i],cy[i]=centroid_2dg(image_data3[i,13:18,13:18])
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
#edg is the edge(in pixels) of the starbox
def psfwidth(image_data,xo,yo,edg):
	psfxw=np.zeros(64)
	psfyw=np.zeros(64)
	lb= (32-edg)/2
	ub= (32+edg)/2
	#ub+=1
	#lb+=1
	stx=np.ndarray((64,5,5))
	np.copyto(stx,image_data[:,lb:ub,lb:ub])
	for i in range(64):
		denom=0.0
		numerx=0.0
		numery=0.0
		for j in range(edg):
			for k in range(edg):
				f=stx[i][j][k]
				numerx+=f*(j-xo[i]+lb)*(j-xo[i]+lb)
				numery+=f*(k-yo[i]+lb)*(k-yo[i]+lb)
		denom=np.nansum(stx[i,:,:])
		widx=numerx/denom
		widy=numery/denom
		widx=np.sqrt(widx)
		widy=np.sqrt(widy)
		psfxw[i]=widx
		psfyw[i]=widy
	
	'''
		denom=0.0
		numerx=0.0
		numery=0.0
		X2= (X-xo[i])**2
		Y2= (Y-yo[i])**2

		if(i==0):
			print X2
			print Y2
		widx = (np.sum(X2*stx[i,:,:])/np.sum(stx[i,:,:]))
		widy = (np.sum(Y2*stx[i,:,:])/np.sum(stx[i,:,:]))
		widx = np.sqrt(widx)
		widy = np.sqrt(widy)
		psfxw[i]=widx
		psfyw[i]=widy
	'''
	return psfxw,psfyw


def psfwidth1(image_data, xo, yo, edg):
	psfxw=np.zeros(64)
	psfyw=np.zeros(64)
	starbox = image_data[:, 13:18, 13:18]
	h,w = np.shape(starbox[0,:,:])
	x = np.arange(0,w)
	y = np.arange(0,h)
	

	X,Y= np.meshgrid(x,y)
	X +=13
	Y +=13
	for i in range(64):
		X2 = (X - xo[i])**2
		Y2 = (Y - yo[i])**2
		
		widx = (np.sum(X2*starbox[i,:,:])/np.sum(starbox[i,:,:]))
		widy = (np.sum(Y2*starbox[i,:,:])/np.sum(starbox[i,:,:]))
		widx = np.sqrt(widx)
		widy = np.sqrt(widy)
		psfxw[i] = widx
		psfyw[i] = widy
	return psfxw, psfyw

# Noise pixel param
def noisepixparam(image_data,edg):
	lb= (32-edg)/2
	ub= (32+edg)/2
	npp=[]
	# Its better to operate on the copy of desired portion of image_data than on image_data itself.
	# This reduces the risk of modifying image_data accidently. Arguements are passed as pass-by-object-reference.
	stx=np.ndarray((64,7,7))
	np.copyto(stx,image_data[:,lb:ub,lb:ub])
	for img in stx:
		#To find noise pixel parameter for each frame. For eqn, refer Knutson et al. 2012
		numer= np.nansum(img)
		numer=np.square(numer)
		denom=0.0
		temp = np.square(img)
		denom = np.nansum(img)
		param= numer/denom
		npp.append(param)
	return npp

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
	normbg.append(bgsum)
	
	#print " normal ", bgsum[5]
	#bg_avg = np.mean(bgsum)
	#bgsum=bgsum/
def normstar(ape_sum,normf):
	starmean=np.nanmean(ape_sum)
	ape_sum=ape_sum/starmean
	normf.append(ape_sum)
	#print min(enumerate(normf), key=operator.itemgetter(1))
	
def normxycent(xo,yo,normx,normy):
	xo=xo/np.nanmean(xo)
	yo=yo/np.nanmean(yo)
	normx.append(xo)
	normy.append(yo)
def normGxycent(gx,gy,normGx,normGy):
	gx=gx/np.nanmean(gx)
	gy=gy/np.nanmean(gy)
	normGx.append(gx)
	normGy.append(gy)

def normpsfwidth(psfwx,psfwy,normpsfwx,normpsfwy):
	psfwx=psfwx/np.nanmean(psfwx)
	psfwy=psfwy/np.nanmean(psfwy)
	normpsfwx.append(psfwx)
	normpsfwy.append(psfwy)	

def normnoisepix(npp,normnpp):
	npp = npp/np.nanmean(npp)
	normnpp.append(npp)

def stackit(normf,normbg,normx,normy,normGx,normGy,normpsfwx,normpsfwy,normnpp):
	normf=np.nanmean(normf,axis=0)
	normbg=np.nanmean(normbg, axis=0)
	normx=np.nanmean(normx,axis=0)
	normy=np.nanmean(normy,axis=0)
	normGx=np.nanmean(normGx,axis=0)
	normGy=np.nanmean(normGy,axis=0)
	normpsfwx=np.nanmean(normpsfwx,axis=0)
	normpsfwy=np.nanmean(normpsfwy,axis=0)
	normnpp=np.nanmean(normnpp,axis=0)
	return normf,normbg,normx,normy,normGx,normGy,normpsfwx,normpsfwy,normnpp


def plotcurve(xax,f,b,X,Y,GX,GY,wx,wy,npp,direc,ct):
	devfactor=2
	fmed=np.nanmedian(f)
	fstdev=np.std(f)
	lb=fmed-devfactor*fstdev
	ub=fmed+devfactor*fstdev
	avoid=[]
	i=0
	for x in (0,57):
		if( f[x] <=lb or f[x]>=ub):
			avoid.append(x)
	print avoid
	fig, axes = plt.subplots(nrows=9, ncols=1, sharex=True)
	fig.set_figheight(8)
	plt.minorticks_on()
	fig.subplots_adjust(hspace = 0.001)
	plt.rc('font', family='serif')
	
	#fig.subplots_adjust(.15,.15,.9,.9,0,0)
	plt.xlim(0,64)
	y_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
	axes[0].plot(xax,f,color='k', mec ='r', marker='x', markevery=avoid,fillstyle='none')
	if 0 not in (avoid):
		axes[0].plot(xax,f,color='k', mec ='g', marker='s', markevery=[0],fillstyle='none')
	if 57 not in (avoid):
		axes[0].plot(xax,f,color='k', mec ='g', marker='s', markevery=[57],fillstyle='none')

	axes[0].plot(xax,f,color='k', mec ='r', marker='x', markevery=[0],fillstyle='none')
	axes[0].set_ylabel(r'$F$',fontsize=16)
	axes[0].yaxis.set_major_formatter(y_formatter)
	axes[0].yaxis.set_major_locator(MaxNLocator(prune='both',nbins=5))

	axes[1].plot(xax,b,color='k', mec ='r', marker='s', markevery=[57],fillstyle='none')
	axes[1].plot(xax,b,color='k', mec ='r', marker='x', markevery=[0],fillstyle='none')
	axes[1].set_ylabel(r'$b$',fontsize=16)
	axes[1].yaxis.set_major_formatter(y_formatter)
	axes[1].yaxis.set_major_locator(MaxNLocator(prune='both',nbins=5))

	axes[2].plot(xax,X,color='k', mec ='r',marker='s', markevery=[57],fillstyle='none')
	axes[2].plot(xax,X,color='k', mec ='r',marker='x', markevery=[0],fillstyle='none')
	axes[2].set_ylabel(r'$x_0$',fontsize=16)
	axes[2].yaxis.set_major_formatter(y_formatter)
	axes[2].yaxis.set_major_locator(MaxNLocator(prune='both',nbins=5))

	axes[3].plot(xax,Y,color='k' , mec ='r', marker='s', markevery=[57],fillstyle='none')
	axes[3].plot(xax,Y,color='k' , mec ='r', marker='x', markevery=[0],fillstyle='none')	
	axes[3].set_ylabel(r'$y_0$',fontsize=16)
	axes[3].yaxis.set_major_formatter(y_formatter)
	axes[3].yaxis.set_major_locator(MaxNLocator(prune='both',nbins=5))

	axes[4].plot(xax,GX,color='k', mec ='r',marker='s', markevery=[57],fillstyle='none')
	axes[4].plot(xax,GX,color='k', mec ='r',marker='x', markevery=[0],fillstyle='none')
	axes[4].set_ylabel(r'$G x_0$',fontsize=16)
	axes[4].yaxis.set_major_formatter(y_formatter)
	axes[4].yaxis.set_major_locator(MaxNLocator(prune='both',nbins=5))

	axes[5].plot(xax,GY,color='k' , mec ='r', marker='s', markevery=[57],fillstyle='none')
	axes[5].plot(xax,GY,color='k' , mec ='r', marker='x', markevery=[0],fillstyle='none')	
	axes[5].set_ylabel(r'$G y_0$',fontsize=16)
	axes[5].yaxis.set_major_formatter(y_formatter)
	axes[5].yaxis.set_major_locator(MaxNLocator(prune='both',nbins=5))

	axes[6].plot(xax,wx,color='k' , mec ='r', marker='s', markevery=[57],fillstyle='none')
	axes[6].plot(xax,wx,color='k' , mec ='r', marker='x', markevery=[0],fillstyle='none')
	axes[6].set_ylabel(r'$\sigma_x$',fontsize=16)
	axes[6].yaxis.set_major_formatter(y_formatter)
	axes[6].yaxis.set_major_locator(MaxNLocator(prune='both',nbins=5))	

	axes[7].plot(xax,wy,color='k' , mec ='r', marker='s', markevery=[57],fillstyle='none')
	axes[7].plot(xax,wy,color='k' , mec ='r', marker='x', markevery=[0],fillstyle='none')
	axes[7].set_ylabel(r'$\sigma_y$', fontsize=16)
	#axes[5].set_xlabel(r'$Frame$ $number$',fontsize=16)
	axes[7].yaxis.set_major_formatter(y_formatter)
	axes[7].yaxis.set_major_locator(MaxNLocator(prune='both',nbins=5))

	axes[8].plot(xax,npp,color='k' , mec ='r', marker='s', markevery=[57],fillstyle='none')
	axes[8].plot(xax,npp,color='k' , mec ='r', marker='x', markevery=[0],fillstyle='none')
	axes[8].set_ylabel(r'$\beta$', fontsize=16)
	axes[8].set_xlabel(r'$Frame$ $number$',fontsize=16)
	axes[8].yaxis.set_major_formatter(y_formatter)
	axes[8].yaxis.set_major_locator(MaxNLocator(prune='both',nbins=5))

	plt.savefig('test/4.5_'+str(ct)+'_'+direc+'.png',bbox_inches='tight',dpi=500)

#note that the outerpath must be complete with terminating '/'
outerpath='/home/hema/Documents/mcgill/handy/aorkeys-20-selected_AORs/'
dirs=os.listdir(outerpath)
print dirs
counter=0
for direc in dirs :
	print direc
	if(counter==10):
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
	print counter
	for filename in glob.glob(os.path.join(path, '*bcd.fits')):
		#print filename
		#if(ct==100):
		#	break
		print ct
		f=fits.open(filename,mode='readonly')
		image_data0=f[0].data
		# convert MJy/str to electron count
		convfact=f[0].header['GAIN']*f[0].header['EXPTIME']/f[0].header['FLUXCONV']
		image_data1=image_data0*convfact
		
		#sigma clip
		image_data2=sigma_clipping(image_data1)
		#for b
		bgnormalize(image_data2,normbg)
		#bg subtract
		image_data3=bgsubtract(image_data2)
		#centroid
		xo, yo = centroid(image_data3)
		#apply gaussian 2d fit and find the centroid
		gx, gy = centroidg2d(image_data3)
		#aperture photmetry
		ape_sum=aperphot(image_data3,2.5,xo,yo)

		
		psfwx,psfwy=psfwidth1(image_data3,xo,yo,5)
		
		npp=noisepixparam(image_data3,7)
		normstar(ape_sum,normf)
		normxycent(xo,yo,normx,normy)
		normGxycent(gx,gy,normGx,normGy)
		normpsfwidth(psfwx,psfwy,normpsfwx,normpsfwy)
		normnoisepix(npp,normnpp)
		ct+=1
		
	print ct

	#Since we are appending and not extending lists, we needn't reshape it
	#normf,normbg,normx,normy,normpsfwx,normpsfwy=reshapelists(normf,normbg,normx,normy,normpsfwx,normpsfwy,ct)
	normf,normbg,normx,normy,normGx,normGy,normpsfwx,normpsfwy,normnpp=stackit(normf,normbg,normx,normy,normGx,normGy,normpsfwx,normpsfwy,normnpp)
	frameno=np.arange(0,64)
	plotcurve(frameno,normf,normbg,normx,normy,normGx,normGy,normpsfwx,normpsfwy,normnpp,direc,ct)
	counter+=1