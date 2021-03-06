# Perform aperture photometry and find stellar flux
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from astropy.utils.data import download_file
from astropy.io import fits
from astropy.stats import sigma_clip
from photutils import CircularAperture
from photutils import aperture_photometry
from numpy import std
import glob
import csv
import os
tossed = 0
def sigmaclip(image_data,fname):
	global badframetable
	global tossed
	sig_clipped_data=sigma_clip(image_data, sigma=4, iters=2, cenfunc=np.median,axis=0)
	for i in range (1,64):
		oldstar=image_data[i,13:19,13:19]
		newstar=sig_clipped_data[i,13:19,13:19]
		truth= newstar==oldstar
		if(truth.sum() <truth.size):
			sig_clipped_data[i,:,:]=np.nan
			badframetable.append([i,fname])
			tossed+=1
	return sig_clipped_data
def bgsubtract(image_data):
	bgsubimg=image_data
	#Creating the mask for central 16 pixels for back
	x=np.ndarray ( shape=(64,32,32), dtype=bool)
	xmask=np.ma.make_mask(x,copy=True, shrink=True, dtype=np.bool)
	xmask[:,:,:]= False
	xmask[:,13:18,13:18]=True
	masked= np.ma.masked_array(image_data, mask=xmask)
	n=0
	#Background subtraction for each frame
	while(n<64):
		#1st frame and 58th frame are outliers
		if(n==0 or n==57):
			n+=1
			continue
		#Excludes NAN values in the calculations
		#bg_avg stores the average flux of the background
		bg_avg=np.median(masked[n,:,:])
		#Subtract the mean 
		bgsubimg[n]= bgsubimg[n,:,:] - bg_avg
		#Setting negative values to 0
	#	bgsubimg[n]= bgsubimg[n].clip(min = 0)
	#	plt.clf()
	#	plt.imshow(bgsubimg[n])
		n+=1
		#plt.savefig(str(n)) #To save the img
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
def badframes(badframetable):
	with open('badframes.csv', 'wb') as outcsv:
		writer = csv.writer(outcsv, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL, lineterminator='\n')
		writer.writerow(['Bad frame number', 'File name'])
		for item in badframetable:
			writer.writerow(item)		# Write item to outcsv
	return 
def tocsv(ape_sum,fname):
	with open('aperture_photometry_auxfiles/stellarflux/'+fname[fname.find('ch2/bcd/')+8:fname.find('.fits')] + '.csv', 'wb') as outcsv:
		writer = csv.writer(outcsv, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL, lineterminator='\n')
		writer.writerow(['Frame number', 'Stellar Flux'])
		frameno= np.arange(1,65)
		temp= [ frameno,ape_sum]
		temp=zip(*temp)
		for item in temp:
			writer.writerow(item)

# Table to store the details of bad frames tossed out during sigma clipping
badframetable=[]
# Number of directories to scan for
ndirec=1

# Starting aor number
direc=20150000
# Stores the file names
fnames=[]
for i in range(ndirec):
	path="/home/hema/Documents/mcgill/XO-3_b_sim_ch2_150702/r"+str(direc)+"/ch2/bcd"		# Resolving path to read file
	fnames.extend(  [filename for filename in glob.glob(os.path.join(path, '*bcd.fits'))])
	direc+=1
# To get the file names in alphabetical order
fnames.sort()
# Number of files scanned
nfiles=len(fnames)

for i in range (nfiles):
	fname=fnames[i]
	f=fits.open(fname)
	image_data=f[0].data
	filtered_data=sigmaclip(image_data,fname[fname.find('ch2/bcd/')+8:])
	bgsubimg=bgsubtract(filtered_data)
	cx,cy= centroid(bgsubimg)
	ape_sum=aperphot(bgsubimg,2.5,cx,cy)
	tocsv(ape_sum,fname)	#Write the calculated stellar flux of 64 frames into csv file

# Make a csv file for the badframes
badframes(badframetable);
'''
stacked=np.median( ape_sum,axis=0)
plt.imshow(stacked, interpolation='none',cmap='gray')

'''

'''
stacked=np.median( filtered_data, axis=0)
plt.imshow(stacked, interpolation='none',cmap='gray')
plt.title('Sigma clipped _'+filename)
plt.xlabel('x pixels')
plt.ylabel('y pixels')
plt.colorbar()
plt.savefig('Sigma clipped _'+filename)
f.close()
'''