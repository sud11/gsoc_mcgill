#--------------------------------------------------------------
#Author: Lisa Dang & Sudarsan 
#Created: 2016-07-13 8:50 AM EST
#Last Modified: 
#Title: Identifying the best combination of boxcar window size
#       and aperture radius for aperture photometry.
#--------------------------------------------------------------


# contour plot

# CONVOLVING FLUXLIST and subtracting that from original flux list
# One boxcar width with all radii
# Next boxcar width with all radii
# etc


# Pass the data through a high pass filter.
# Plot for F RMS value for different aperture sizes
# For a high pass filter, do a boxcar smooth and subtract the smoothed data from raw data

import scipy.interpolate
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
def sigma_clipping(image_data):#,fname):
	fin=0
	tempo=[]
	toss=[]
	tot=0
	
	#Moving window of width 50
	'''
	tem,to=sigclip(image_data[0:10])
	tempo.extend(tem)
	toss.extend(np.asarray(list(to)))

	for i in range (len(image_data)-20):
		#print i,i+20,(i+i+20)/2
		temp,to=sigclip(image_data[i:i+20])
		if(np.isnan(temp[10,15,15])):
			toss.extend([i+10])
			tot+=1
		tempo.append(temp[10])
		fin=i

	tem,to=sigclip(image_data[-10:])
	tempo.extend(tem)
	to=np.asarray(list(to))+fin+10
	toss.extend(to)
	
	

	tempo=np.asarray(tempo)
	'''
	
	#Fixed window of 20
	for i in range (len(image_data)/20):
		temp,to=sigclip(image_data[i*20:(i+1)*20])
		tempo.extend(temp)
		tot+=len(to)
		if(len(to)>0):
			#toss.extend(np.asarray(np.asarray(to)+i*20)
			#toss.extend( t(list(to)+i*20))
			to=np.asarray(list(to))+i*20
			toss.extend(to)
		fin=i
	print len(allstackd1)
	print len(tempo)
	if( len(tempo) < len(allstackd1)):
		tem,to=sigclip(allstackd1[(fin+1)*20:])
		to=np.asarray(list(to))+(fin+1*20)
		toss.extend(to)
		tempo.extend(tem)
	tempo=np.ma.asarray(tempo)
	print len(toss),'toss len'
	
	return tempo,toss

def sigclip(image_data):
	#x=np.ndarray( shape=(32,32), dtype=bool)
	#xmask=np.ma.make_mask(x,copy=True, shrink=True, dtype=np.bool)
	#xmask[:,:]= True
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
	#print toss
	#print 'ptoss',ptoss
	#print len(ptoss)
	return sig_clipped_data,ptoss

def bgsubtract(image_data):
	bgsubimg=image_data
	x=np.ndarray ( shape=(len(image_data),32,32), dtype=np.bool)
	x[:,:,:]=False
	x[:,14:18,14:18]=True
	#xmask=np.ma.make_mask(x,copy=True, shrink=True,dtype=np.bool)
	#xmask[:,:,:] = False
	#xmask[:,14:18,14:18]=True
	#xmask[:,0:1,:]=True
	masked= np.ma.masked_array(bgsubimg, mask=x)
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

def getflist(allstackd,op,apr=2.5):
	'''
	allstackd=np.asarray(allstackd)
	temp=np.apply_along_axis(lambda m: convolve(m, Box1DKernel(50)), axis=0, arr=allstackd)
	'''
	image_data1,toss= sigma_clipping(allstackd)
	image_data2= bgsubtract(image_data1)
	cx, cy = centroid(image_data2)
	temp=0
	flist=aperphot(image_data2,apr,cx,cy,op)
	return flist

# Flist after high pass filter
def highpassflist(allstackd,op,apr=2.5,bxw=64):

	np.asarray(allstackd)
	flist1=getflist(allstackd,op,apr)
	#temp=np.apply_along_axis(lambda m: convolve(m, Box1DKernel(bxw),boundary='extend'), axis=0, arr=allstackd)
	flist2=convolve(flist1, Box1DKernel(bxw),boundary='extend')
	flist1=np.asarray(flist1)
	flist2=np.asarray(flist2)
	flist3= flist1-flist2	#raw-smooth
	'''	
	# Binning version
	binstak=map(binning_data,flist1)
	binstak=np.asarray(binstak)
	print 'lengths',len(binstak),len(binstak[0])
	flist1=binstak
	'''

	'''
	sflist = map(normstar,smooth)
	sflist = map( lambda t: t-1, sflist)
	panel=map(lambda t: np.sqrt(np.nanmean(np.square(t))), sflist)
	'''
	return flist1,flist2,flist3

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
	print 'testing plot3panels'
	print len(xax)
	print len(p1)

	fig, axes = plt.subplots(nrows=4, ncols=1, sharex=True)
	plt.minorticks_on()
	fig.subplots_adjust(hspace = 0.001)
	plt.rc('font', family='serif',serif='Times')
	y_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
	
	axes[0].plot(xax,p1,'-',c='k',mec='r',marker='x', markevery=toss,fillstyle='none')
	axes[0].set_ylabel(r'$original$',fontsize=13)
	axes[0].yaxis.set_major_formatter(y_formatter)
	axes[0].xaxis.set_major_formatter(y_formatter)
	axes[0].yaxis.set_major_locator(MaxNLocator(prune='both',nbins=5))

	axes[1].plot(xax,p2,'-',c='k',mec='r',marker='x', markevery=toss,fillstyle='none')
	axes[1].set_ylabel(r'$smoothed$',fontsize=13)
	axes[1].yaxis.set_major_formatter(y_formatter)
	axes[1].xaxis.set_major_formatter(y_formatter)
	axes[1].yaxis.set_major_locator(MaxNLocator(prune='both',nbins=5))

	axes[2].plot(xax,p3,'-',c='k',mec='r',marker='x', markevery=toss,fillstyle='none')
	axes[2].set_ylabel(r'$flattened$',fontsize=13)
	axes[2].yaxis.set_major_formatter(y_formatter)
	axes[2].xaxis.set_major_formatter(y_formatter)
	axes[2].yaxis.set_major_locator(MaxNLocator(prune='both',nbins=5))

	axes[3].plot(xax,p4,'-',c='k',mec='r',fillstyle='none',label='Soft')
	axes[3].set_ylabel(r'$\sigma$ $clipped$',fontsize=13)
	axes[3].set_xlabel(r'$time$ $(hours)$',fontsize=13)
	axes[3].yaxis.set_major_formatter(y_formatter)
	axes[3].xaxis.set_major_formatter(y_formatter)
	axes[3].yaxis.set_major_locator(MaxNLocator(prune='both',nbins=5))
	#axes[2].legend(numpoints=1)
	plt.savefig('test/'+str(ct)+'con_fig2lx.png',bbox_inches='tight',dpi=500)
	

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
	plt.savefig('test/'+str(ct)+'updchanges.png',bbox_inches='tight',dpi=500)

def masterplot(tim,p1,p2,p3,p4,toss,xax,yax,zax,ct):
	plot3panels(tim,p1,p2,p3,p4,toss,ct)
	plt.clf()
	plotandsave(xax,yax,zax,ct)


	
def plotandsave(xax,yax,zax,ct):

	'''
	plt.clf()
	#http://stackoverflow.com/questions/9008370/python-2d-contour-plot-from-3-lists-x-y-and-rho
	x=np.asarray(xax)
	y=np.asarray(yax)
	z=np.asarray(zax)
	xi, yi = np.linspace(x.min(), x.max(), 100), np.linspace(y.min(), y.max(), 100)
	xi, yi = np.meshgrid(xi, yi)
	# Interpolate
	rbf = scipy.interpolate.Rbf(x, y, z, function='linear')

	zi = rbf(xi, yi)
	print zi
	#plt.contour(zi)
	#line_colours=('BlueViolet', 'Crimson', 'ForestGreen','Indigo', 'Tomato', 'Maroon')
	CS=plt.contourf(zi,extend='neither',vmin=z.min(), vmax=z.max(), origin='lower',extent=[x.min(), x.max(), y.min(), y.max()],aspect='auto')
	'''

	'''
	Line contour
	#####
	import matplotlib.pyplot as plt
	from mpl_toolkits.mplot3d import Axes3D
	from matplotlib import cm

	fig = plt.figure()
	#ax = fig.add_subplot(111, projection='3d')
	cset = plt.contour(xax, yax, zax)#, cmap=cm.coolwarm)
	#####
	#plt.show()
	'''
	'''
	from mpl_toolkits.mplot3d import axes3d
	from matplotlib import cm
	print len(xax)
	print len(yax)
	X,Y=np.meshgrid(xax,yax)
	print X
	print Y
	Z=np.asarray(zax)
	plt.contourf(X,Y,Z)
	'''
	X,Y=np.meshgrid(xax,yax)
	#Below line for the normal contourf plots
	CS=plt.contourf(X,Y,zax,origin='lower')#,extent=[np.min(xax),np.max(xax), np.min(yax), np.max(yax)],aspect=None)

	plt.colorbar()
	#xax=list(xax)
	#yax=list(yax)
	#print xax
	#print yax

	#plt.scatter(xax,yax)
	#plt.plot(xax)
	#cs = plt.contour(scalar_field)
	#plt.scatter([5,56,57],[2.5,3.5,4.5],zorder=1)
	#plt.clabel(CS,inline=1, fontsize=10)
	#plt.scatter(x, y, c=z)
	#plt.show()
	#plt.xlim(np.min(xax),np.max(xax))
	
	#ax.get_xaxis().get_major_formatter().set_useOffset(False)
	#ax.get_yaxis().get_major_formatter().set_useOffset(False)
	#plt.xlabel('boxcar width',fontsize=16)
	#plt.ylabel('aperture radius',fontsize=16)
	#plt.title('flattened')
	plt.savefig('test2/con2_'+str(np.min(xax))+' to '+str(np.max(xax))+'_'+str(ct)+'ab_fig1lxT.png',dpi=500)


def masterfunc1(op,apr,bxw):
	xax=[2.5]
	global allstackd1
	allstackd1=np.asarray(allstackd1)
	print 'allstackd1',len(allstackd1)
	panel1,panel2,smooth1=highpassflist(allstackd1,op,apr,bxw)

	print "Ok till here,",bxw
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
	scl=panels[3]
	scl=np.asarray(scl)
	'''
	scl=map(abs,scl)
	
	scl=normstar(scl)
	scl-=1
	'''
	z=np.sqrt(np.nanmean(np.square(scl)))
	return z
	#return panels,toss


def supportplot(op):
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
	print 'tossed in suppo',toss
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
	for direc in dirs:
		if(counter==1):
			break
		path=outerpath+direc+'/ch2/bcd'
		print path
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

start_time=time()
allstackd1=[]
fnames,ftimes=getfnames()
ct=0
t=[]
print 'start time :',start_time
for filename in fnames:
	#print filename
	if(ct==220):
		break
	f=fits.open(filename,mode='readonly')
	#print filename
	image_data0=f[0].data
	t.extend( np.linspace(f[0].header['AINTBEG'],f[0].header['ATIMEEND'],64) )
	# convert MJy/str to electron 
	convfact=f[0].header['GAIN']*f[0].header['EXPTIME']/f[0].header['FLUXCONV']
	image_data1=image_data0*convfact
	allstackd1.extend(image_data1)
	ct+=1

#convert time from sec to hours
t=np.asarray(t)
sectohrs=1.0/3600.0
tim=t*sectohrs

aplist=[2,2.5,3,3.5,4,4.5]
#bxlist=np.arange(50,80,5)
#bxlist=[2,50,100,150,200,250,300,350,400,450,500,550,600,650,700,750,800,900,950,1000]
bxlist=[2,50,100,150,200,250,300,350,400,450,500,550,600]


panels,toss=supportplot(1)

zlist=map(lambda bx : map(lambda t: masterfunc1(1,aplist[t],bxlist[bx]),np.arange(len(aplist))) , np.arange(len(bxlist)))
np.save('tempzlist',zlist)
print 'FINE HERE TOO'
print 'zlist',zlist
print bxlist
print 'end_time', time()-start_time
#plot3panels(time,panels[0],panels[1],panels[2],panels[3],toss,ct)
#plotandsave(bxlist,aplist,zlist,ct)
masterplot(tim,panels[0],panels[1],panels[2],panels[3],toss,aplist,bxlist,zlist,ct)

#plot3panels(tim,panels[0],panels[1],panels[2],panels[3],toss,ct)