# program to read series of images from the Spitzer Data Challenge, identify the brightest pixel in the images and produces a lightcurve (flux vs.time) for that pixel.
# The mean of the maximum flux per frame for 64 frames is taken
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

mtable=[]	#This list holds the values for making the CSV table

def updatetable(xpix, ypix,flux,obstime):
	obsdatetime=dateutil.parser.parse(obstime)  #To make it as a datetime object
	mtable.append([xpix,ypix,flux,obsdatetime])
	
	return
	
def makecsv(mtable):
	mtable= sorted(mtable,key=lambda y: y[3])	# To sort based on timestamp
	with open('brightpixeldata.csv', 'wb') as outcsv:   	    #configure writer to write standard csv file
	    writer = csv.writer(outcsv, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL, lineterminator='\n')
	    writer.writerow(['PIXEL_X', 'PIXEL_Y', 'FLUX', 'TIMESTAMP'])
	    for item in mtable:
	        writer.writerow(item)		# Write item to outcsv
	             
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
		date_object = datetime.strptime(timestmp[i] ,'%Y-%j-%d %H:%M:%S.%f')
		timestmp[i]=time.mktime(date_object.timetuple())		#Convert the timestamp to seconds. 
		obsdate.append( timestmp[i])
		i=i+1
	print "CSV file read"
	print "Plotting light curve"
	#plotlightcurve(count,fluxlist[1:])
	plotlightcurve(obsdate,fluxlist[1:])

							
def plotlightcurve(xaxis,yaxis):									#Plot (flux vs. time) light curve
	pyplot.plot(xaxis,yaxis,'ro')
	pyplot.title('Light curve for brightest pixel in the sample')
	pyplot.xlabel('Time of observation (s)')
	pyplot.ylabel('Flux (MJ/st)')
	pyplot.grid(True)
	pyplot.show()
	       
# MAIN
	
direc=20150005			#AOR directory name for including in path initialised with the starting value of directory
ndirec=1			# Number of directories to collect images from
ct =0;	
for i in range(ndirec):
	path = "/home/hema/Documents/mcgill/XO-3_b_sim_ch2_150702/r"+str(direc)+"/ch2/bcd"		# Resolving path to read file
	for filename in glob.glob(os.path.join(path, '*bcd.fits')):
		ct = ct +1
		f=pyfits.open(filename)
		image_data = f[0].data
		header= f[0].header	
		sumofmax=0.0
		n=0
		while( n < 64 ):
			xpix= 15											  # Hardcoding brightest pixel value based on manual observation
			ypix= 15
			maxn=image_data[n][xpix][ypix]											  # Intensity is the intensity of (xpix,ypix)
			sumofmax += maxn
			n=n+1
		avgflux = sumofmax/64							 # Finding the mean of the flux of brightest pixel for the 64 frames
		updatetable(xpix,ypix,avgflux,header['DATE_OBS'] )
	direc+=1
	print path					# Printing the path of directories scanned for images

print "Number of images processed:",ct
#print "UPD:", upd
makecsv(mtable)
print "CSV file created"
readcsv("brightpixeldata.csv")

