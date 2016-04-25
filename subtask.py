# Data used : XO-3_b_sim_ch2_150702 
# Simulated Data is used for this program. 
# Link to data : http://conference.ipac.caltech.edu/iracaas224/data_challenge/data/XO-3_b_sim_ch2_150702.tgz

import pyfits 
import numpy as np
import glob
import os 
import pandas
import time
from numpy import unravel_index
from operator import itemgetter
import csv
import dateutil.parser
import matplotlib.pyplot as plt
from matplotlib import pyplot
from datetime import datetime
import operator
from collections import defaultdict

#This list holds the values for making the CSV table
mtable=[]

# findbrightestpixel() takes the image_data and for each frame, finds the brightest pixel. The corresponding (x,y) of the brightest pixel is stored in a dictionary and its count is incremented. After processing the 64 frames, the coordinates (x,y) which has been the brightest in maximum number of frames is returned.

def findbrightestpixel(img_data,dict):

	n=0
	while( n < 64 ):
		# To get index which holds the brightest pixel. This propogates through NAN values	
		pixcoord= unravel_index(image_data[n].argmax(), image_data[n].shape)  
		xpix= pixcoord[0]																		  
		ypix= pixcoord[1]
		dict[(xpix,ypix)]+=1
		# To get the value of the brightest pixel.
		maxn=np.amax(image_data[n])
		n=n+1
	pixcoord= max(dict.iteritems(), key=operator.itemgetter(1))[0]
	return pixcoord
	
	
# updatetable() takes updates the mtable by creating a row with the parameters passed
def updatetable(xpix, ypix,flux,obstime):
	#To make it as a datetime object
	obsdatetime=dateutil.parser.parse(obstime)
	mtable.append([xpix,ypix,flux,obsdatetime])
	
	return

# makecsv() takes a numpy array and creates a CSV file for the same with proper headers
def makecsv(mtable):
	# To sort based on timestamp
	mtable= sorted(mtable,key=lambda y: y[3])	
	with open('brightpixeldata.csv', 'wb') as outcsv:
		#configure writer to write standard csv file
		writer = csv.writer(outcsv, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL, lineterminator='\n')
	    writer.writerow(['PIXEL_X', 'PIXEL_Y', 'FLUX', 'TIMESTAMP'])
	    for item in mtable:
	    	# Write item to outcsv
			writer.writerow(item)		

# readcsv() reads the CSV file and stores the data in lists
def readcsv(csvfilename):
	colnames = ['PIXEL_X','PIXEL_Y','FLUX','TIMESTAMP']
	path = "/home/hema/Documents/mcgill/main/"+csvfilename
	data = pandas.read_csv(path, names=colnames)
	# Read the CSV file data and store in lists
	x = data.PIXEL_X.tolist()
	y = data.PIXEL_Y.tolist()
	fluxlist= data.FLUX.tolist()
	timestmp = data.TIMESTAMP.tolist()
	i=1
	obsdate=[]
	count=[]				
	while ( i < len( timestmp) ):
		count.append(i)
		try:
			date_object = datetime.strptime(timestmp[i] ,'%Y-%j-%d %H:%M:%S.%f')	
		except ValueError:
			date_object = datetime.strptime(timestmp[i] ,'%Y-%j-%d %H:%M:%S')
		#Convert the timestamp to seconds.
		timestmp[i]=time.mktime(date_object.timetuple()) 
		obsdate.append( timestmp[i])
		i=i+1
	
	print "CSV file read"
	print "Plotting light curve"
	#Plot (flux vs. time) light curve
	plotlightcurve(obsdate,fluxlist[1:])

# plotlightcurbe() takes the axes lists for which a lightcurve will be plotted

def plotlightcurve(xaxis,yaxis):									
	pyplot.plot(xaxis,yaxis,'ro')
	# X axis : Time of observation ( taken from FITS header ) in seconds
	# Y axis : Flux in MJ/st
	pyplot.title('Light curve for brightest pixel in the sample')   
	pyplot.xlabel('Time of observation (s)')
	pyplot.ylabel('Flux (MJ/st)')
	pyplot.grid(True)
	pyplot.show()
	pyplot.close()

	       
# MAIN

#AOR directory name for including in path initialised with the starting value of directory
direc=20150000
# Number of directories to collect images from
ndirec=1			
ct =0;	

# Dictionary that holds the coordinates of potential brightest pixels and their count
dict = defaultdict(int)

for i in range(ndirec):
	# Resolving path to read file
	path = "/home/hema/Documents/mcgill/XO-3_b_sim_ch2_150702/r"+str(direc)+"/ch2/bcd"		
	for filename in glob.glob(os.path.join(path, '*bcd.fits')):
		f=pyfits.open(filename)
		image_data = f[0].data
		# bpx holds coordinates of the brightest pixel
		bpx= findbrightestpixel(image_data,dict)				
# x coordinate of brightest pixel
xpix=bpx[0]										
# y coordinate of brightest pixel
ypix=bpx[1]										
		
		
		
for i in range(ndirec):
	# Resolving path to read file
	path = "/home/hema/Documents/mcgill/XO-3_b_sim_ch2_150702/r"+str(direc)+"/ch2/bcd"		
	for filename in glob.glob(os.path.join(path, '*bcd.fits')):
		ct = ct +1
		f=pyfits.open(filename)
		image_data = f[0].data
		header= f[0].header	
		sumofmax=0.0
		n=0
		while( n < 64 ):
			# Intensity is the intensity of (xpix,ypix)
			maxn=image_data[n][xpix][ypix]											  
			sumofmax += maxn
			n=n+1
		# Finding the mean of the flux of brightest pixel for the 64 frames
		avgflux = sumofmax/64							 
		updatetable(xpix,ypix,avgflux,header['DATE_OBS'] )
	direc+=1
	# Printing the path of directories scanned for images
	print path					

print "Number of images processed:",ct
makecsv(mtable)
print "CSV file created"
readcsv("brightpixeldata.csv")

