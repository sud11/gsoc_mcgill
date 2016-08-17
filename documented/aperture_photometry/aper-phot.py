#--------------------------------------------------------------
#Author: Lisa Dang & Sudarsan 
#Created: 2016-03-18 1:21 AM EST
#Last Modified: 
#Title: Observation and Photometry of XO-3b
#--------------------------------------------------------------

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import operator
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

def sigma_clipping(image_data,filenb,fname):
	global tossed
	global badframetable
	x=np.ndarray( shape=(32,32), dtype=bool)
	xmask=np.ma.make_mask(x,copy=True, shrink=True, dtype=np.bool)
	xmask[:,:]= True
	sig_clipped_data=sigma_clip(image_data, sigma=4, iters=2, cenfunc=np.median,axis=0)
	for i in range (1,64):
		# oldstar is the original central 7x7 pixels of frame i
		oldstar=image_data[i,12:19,12:19]
		# newstar is the central 7x7 pixels of the sigma-clipped frame i
		newstar=sig_clipped_data[i,12:19,12:19]
		# truth is the matrix which has 'False' value for the values that are masked (i.e. sigma-clipped pixels)
		truth= newstar==oldstar
		# truth.sum() returns the number of values which is 'True'
		if(truth.sum() <truth.size):
			#mask the frame which which has atleast one sigma-clipped pixel in the central 7x7 box
			sig_clipped_data[i,:,:]=np.ma.masked_array(sig_clipped_data[i,:,:], mask=xmask)
			#Update the table that keeps track of the bad frames
			badframetable.append([i,filenb,fname])
			tossed+=1
	return sig_clipped_data
	
def bgsubtract(image_data):
	bgsubimage_data=np.zeros((64,32,32))
	x=np.ndarray( shape=(64,32,32), dtype=bool)
	xmask=np.ma.make_mask(x,copy=True, shrink=True, dtype=np.bool)
	xmask[:,:,:]= False
	# Mask the central 5x5 pixels which has the star
	xmask[:,13:18,13:18]=True
	xmask[:,0:1,:]=True #Removing top defective rows 
	masked= np.ma.masked_array(image_data, mask=xmask)
	bg_err = np.zeros(64)	
	n=0
	while(n<64):
		#bg_avg stores the median background value of frame i
		bg_avg=np.ma.median(masked[n,:,:])
		bg_err[n] = np.ma.var(masked[n,:,:])
		bgsubimage_data[n]=image_data[n,:,:] - bg_avg
		n+=1
	return bgsubimage_data, bg_err

def printstar(image_data,n):
	print (image_data[n,13:18,13:18])

def centroid(image_data):
	#Weighted flux method to find the centroids
	cx=np.zeros(64)
	cy=np.zeros(64)
	starbox = image_data[:, 13:18, 13:18]
	h,w = np.shape(starbox[0,:,:])
	x = np.arange(0,w)
	y = np.arange(0,h)
	X,Y = np.meshgrid(x,y)
	for i in range(64):
		cx[i] = (np.sum(X*starbox[i,:,:])/np.sum(starbox[i,:,:])) + 13
		cy[i] = (np.sum(Y*starbox[i,:,:])/np.sum(starbox[i,:,:])) + 13
	return cx, cy

def badcentroid(image_data, filenb,fname, cx, cy):
	global tossed
	global badframetable
	filtered = image_data
	cx_median=np.median(cx)
	cx_std=np.std(cx)
	cy_median=np.median(cy)
	cy_std=np.std(cy)
	x=np.ndarray( shape=(32,32), dtype=bool)
	xmask=np.ma.make_mask(x,copy=True, shrink=True, dtype=np.bool)
	xmask[:,:]= True
	for i in range(len(cx)):
		if ((abs(cx_median-cx[i]) > 4*cx_std) | (abs(cy_median-cy[i]) > 4*cy_std)):
			filtered[i,:,:]=np.ma.masked_array(filtered[i,:,:], mask=xmask)
			badframetable.append([i,fname])
			tossed+=1
			#print(tossed)
	return filtered

def photometry(image_data, ape_radius, cx, cy, bg_err):
	ape_sum=np.zeros(64)
	ape_sum_err=np.zeros(64)
	for i in range(64):
		position=[cx[i], cy[i]]
		aperture=CircularAperture(position, r=ape_radius)
		#http://photutils.readthedocs.io/en/latest/photutils/aperture.html#error-estimation
		#aperture_sum_error can be accessed only after passing a error parameter
		data_error = bg_err[i]
		phot_table=aperture_photometry(image_data[i,:,:],aperture, error=data_error, pixelwise_error=False)
		ape_sum[i]=phot_table['aperture_sum']
		ape_sum_err[i]=phot_table['aperture_sum_err']
	#print (ape_sum.shape)
	return ape_sum, ape_sum_err

def binning_data(data, h, w):
	reshaped_data = data.reshape((h,w))
	binned_data=np.ma.median(reshaped_data, axis=1)
	binned_data_std=np.std(reshaped_data, axis=1)
	return binned_data, binned_data_std

def frame_discard(data, h, w):
	masked = np.zeros((h,w))
	a = np.ndarray(shape=(h,w), dtype=bool)
	mask0 = np.ma.make_mask(a,copy=True, shrink=True, dtype=np.bool)
	mask0[:,:] = False
	mask0[:,0] = True #masking 0th frame
	reshaped_data = data.reshape((h,w))
	masked = np.ma.masked_array(reshaped_data, mask=mask0)
	reshaped_masked = masked.reshape((1, h*w))
	return reshaped_masked

#makes a csv file for the badframes discarded in sigmaclipping
def badframes(badframetable):
	with open('badframes.csv', 'wb') as outcsv:
		writer = csv.writer(outcsv, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL, lineterminator='\n')
		writer.writerow(['Bad frame number', 'File name'])
		for item in badframetable:
			writer.writerow(item)		# Write item to outcsv
	return 

def getfnames(root="/home/hema/Documents/mcgill/handy/aorkeys-20-selected_AORs/",nfiles=1):
	dirs=os.listdir(root)
	print dirs
	counter=0
	ct=0
	# nametime is a dictionary that holds the path of file as key and its time as value
	# ie. nametime['home/fitsfiles/SPITZER_blah_bla_0001.fits'] = 2312323.34
	nametime={}
	for direc in dirs:
		path=root+direc+'/ch2/bcd'
		ct=0
		if(counter==1):
			break
		for filename in glob.glob(os.path.join(path, '*bcd.fits')):
			try:
				f=fits.open(filename,mode='readonly')	
			except IOError:
				continue
			nametime[filename] = f[0].header['AINTBEG']
			ct+=1
		counter+=1
	#sorting nametime based on time
	nametime=sorted(nametime.items(), key=operator.itemgetter(1))
	#fnames is the list of paths to filenames
	# to read files
	# for filename in fnames :
	# 	f=fits.open(filename)
	fnames= [x[0] for x in nametime]
	return fnames


tossed=0
badframetable=[]
t=[]
bg_err=[]
aperture_sum=[]
aperture_sum_err=[]
xo=[]
yo=[]
psfxw = []
psfyw = []

fnames=[]
'''
path="/home/hema/Documents/mcgill/handy/aorkeys-20-selected_AORs/r46471168/ch2/bcd"		# Resolving path to read file
#path="/home/hema/Documents/mcgill/XO-3_b_sim_ch2_150702/"+AOR_list[i]+"/ch2/bcd"
fnames.extend(  [filename for filename in glob.glob(os.path.join(path, '*bcd.fits'))])
#direc+=1

fnames.sort()
'''
fnames=getfnames("/home/hema/Documents/mcgill/handy/aorkeys-20-selected_AORs/")
fnames=np.asarray(fnames)
i=0
for filename in fnames :
	#get data
	'''
	print 'filename:',filename
	if(i==1):
		break
	'''
	print i
	hdu_list=fits.open(filename)
	image_data0=hdu_list[0].data
	#get time 
	t.extend(np.linspace(hdu_list[0].header['AINTBEG'],hdu_list[0].header['ATIMEEND'],64))
	# convert MJy/str to electron count
	convfact=hdu_list[0].header['GAIN']*hdu_list[0].header['EXPTIME']/hdu_list[0].header['FLUXCONV']
	# convert electron count to Mjy/str
	factor = - hdu_list[0].header['PXSCAL1']*hdu_list[0].header['PXSCAL2']*(1/convfact) 
	image_data1=convfact*image_data0
	#sigma clip
	image_data2=sigma_clipping(image_data1,i,filename[filename.find('ch2/bcd/')+8:])

	#bg subtract
	image_data3, bg_err_tmp=bgsubtract(image_data2)
	#print image_data3[10,13:17,13:17]

	bg_err.extend(factor*bg_err_tmp)
	#get centroid
	xo_tmp, yo_tmp = centroid(image_data3)
	'''
	print xo_tmp[0:10]
	print yo_tmp[0:10]
	'''

	xo.extend(xo_tmp)
	yo.extend(yo_tmp)
	#filter frame with bad centroid
	image_data4=badcentroid(image_data3,i,filename, xo_tmp, yo_tmp)
	#aperture_photometry and flux conversion to MJy/str
	ape_sum_tmp, ape_sum_err_tmp = photometry(image_data4, 2.5, xo_tmp, yo_tmp, bg_err_tmp)
	#sigma clipping flux outliers (entire frame "fishy")
	ape_sum_tmp2 = sigma_clip(ape_sum_tmp, sigma=4, iters=2, cenfunc=np.median)
	
	aperture_sum.extend(factor*ape_sum_tmp2)
	aperture_sum_err.extend(factor*ape_sum_err_tmp)
	#getting psf width
	#psfwxtmp, psfwytmp=psfwidth(image_data3,xo_tmp,yo_tmp,7)
	#psfxw.extend(psfwxtmp)
	#psfyw.extend(psfwytmp)
	i+=1

nfiles=len(fnames)
# Make a csv file for the badframes
#badframes(badframetable);
#badframetable=np.asarray(badframetable)
#np.savetxt('badframe.txt', badframetable,header="Frame in .fits (0 to 63), file nb, file name")

#convert time from sec to hours

sectohrs=1.0/3600.0

t=np.asarray(t)
time=t*sectohrs

aperture_sum=np.asarray(aperture_sum)
time=np.asarray(time)
xo=np.asarray(xo)
yo=np.asarray(yo)
#psfxw=np.asarray(psfxw)
#psfyw=np.asarray(psfyw)
aperture_sum_err=np.asarray(aperture_sum_err)
#tossed respective xo, yo, aperture sum, time of tossed out frames
#index_tossed=np.multiply(badframetable[0, :])

# tossing out 0th frame data
#time0 = frame_discard(time, nfiles, 64)
#bg_err0 = frame_discard(bg_err, nfiles, 64)
#aperture_sum0 = frame_discard(aperture_sum, nfiles, 64)
#aperture_sum_err0 = frame_discard(aperture_sum_err, nfiles, 64)
#xo0 = frame_discard(xo, nfiles, 64)
#yo0 = frame_discard(yo, nfiles, 64)

#binning for data analysis and data visualisation
binned_aperture_sum, binned_aperture_sum_std = binning_data(aperture_sum, nfiles, 64)
binned_time, binned_time_std = binning_data(time, nfiles, 64)
binned_xo, binned_xo_std = binning_data(xo, nfiles, 64)
binned_yo, binned_yo_std = binning_data(yo, nfiles, 64)
#binned_psfxw, binned_psfxw_std = binning_data(psfxw, nfiles, 64)
#binned_psfyw, binned_psfyw_std = binning_data(psfyw, nfiles, 64)

#save above data into file

#np.savetxt('ch2_datacube_full_AORs579.dat', np.c_[aperture_sum, aperture_sum_err, time, xo, yo],header="aperture_sum, aperture_sum_err, time, xo, yo")
#np.savetxt('ch2_datacube_binned_AORs579.dat', np.c_[binned_aperture_sum, binned_aperture_sum_std,binned_time, binned_time_std,binned_xo, binned_xo_std,binned_yo, binned_yo_std, binned_psfxw, binned_psfxw_std, binned_psfyw, binned_psfyw_std],header="binned_aperture_sum, binned_aperture_sum_std,binned_time, binned_time_std,binned_xo, binned_xo_std,binned_yo, binned_yo_std, binned_psfxw, binned_psfxw_std, binned_psfyw, binned_psdyw_std")

#np.savetxt('ch2_datacube_binned_AORs201.dat', )
print fnames[0:5]
formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)

#Plots
fig, axes = plt.subplots(nrows=3, ncols=1, sharex=True)
fig.suptitle("r46482688", fontsize="x-large")
#fig.suptitle("CoRoT-2b", fontsize="x-large")
axes[0].plot(binned_time, binned_aperture_sum,'k+')
axes[0].set_ylabel("Stellar Flux (MJy/pixel)")
axes[0].yaxis.set_major_formatter(formatter)

axes[1].plot(binned_time, binned_xo, '+')
axes[1].set_ylabel("$x_0$")
axes[1].yaxis.set_major_formatter(formatter)

axes[2].plot(binned_time, binned_yo, 'r+')
axes[2].set_xlabel("Time since IRAC turn-on (Hrs)")
axes[2].set_ylabel("$y_0$")
axes[2].yaxis.set_major_formatter(formatter)
axes[2].xaxis.set_major_formatter(formatter)
fig.subplots_adjust(wspace=0)
fig.savefig("4.5test/"+"r46482688",dpi=500)

