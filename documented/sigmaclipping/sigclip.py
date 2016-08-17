#--------------------------------------------------------------
#Author: Lisa Dang & Sudarsan 
#Created: 2016-03-18 1:21 AM EST
#Last Modified: 
#Title: To justify why 1st and 58 frames are outliers
#Test data - 3.6 AOR ['r46483456', 'r46482944', 'r46483200', 'r46482688']
#            4.5 AOR ['r46470912', 'r46469120','r46468608','r46470400']
#Tested by - Sudarsan
#--------------------------------------------------------------

import numpy as np
import os, sys
import glob
import operator
from astropy.io import fits
from astropy.stats import sigma_clip

def sigma_clipping(image_data,ch=0,wlen=50):#,fname):
	fin=0
	tempo=[]
	toss=[]
	tot=0
	if(ch==1):
		#Moving window
		mid=wlen/2
		tem,to=sigclip(image_data[0:mid])	
		tempo.extend(tem)
		toss.extend(np.asarray(list(to)))

		for i in range (len(image_data)-wlen):
			#print i,i+20,(i+i+20)/2
			print i
			temp,to=sigclip(image_data[i:i+wlen])
			if(np.isnan(temp[mid,15,15])):
				toss.extend([i+mid])
				tot+=1
			tempo.append(temp[mid])
			fin=i

		tem,to=sigclip(image_data[-mid:])
		tempo.extend(tem)
		to=np.asarray(list(to))+fin+mid
		toss.extend(to)

		tempo=np.asarray(tempo)
	
	else:
		#Fixed window of 20
		for i in range (len(image_data)/wlen):
			temp,to=sigclip(image_data[i*wlen:(i+1)*wlen])
			tempo.extend(temp)
			tot+=len(to)
			if(len(to)>0):
				to=np.asarray(list(to))+i*wlen
				toss.extend(to)
			fin=i
		if( len(tempo) < len(image_data)):
			tem,to=sigclip(image_data[(fin+1)*wlen:])
			to=np.asarray(list(to))+(fin+1*wlen)
			toss.extend(to)
			tempo.extend(tem)
		tempo=np.ma.asarray(tempo)
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
	return fnames

fnames = getfnames();
allstackd1=[]
for filename in fnames:
	f=fits.open(filename,mode='readonly')
	image_data0=f[0].data
	allstackd1.extend(image_data0)
print 'len all stacked:',len(allstackd1)

s1,t1=sigma_clipping(allstackd1,0)
print 'len all stacked:',len(allstackd1)
print 'len moving window', len(s1)
print 'tossed in moving',t1

s2,t2=sigma_clipping(allstackd1,1)

print 'len all stacked:',len(allstackd1)
print 'len moving window', len(s2)
print 'tossed in moving',t2

print 'len fixed window', len(s2)
print 'tossed in fixed window',t2
