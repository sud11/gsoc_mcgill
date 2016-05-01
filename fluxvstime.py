#To generate flux vs time light curve for an AOR
import pyfits 
import numpy as np
import glob
import os 
import pandas
import time
from numpy import unravel_index
from operator import itemgetter
import dateutil.parser
import matplotlib.pyplot as plt
from matplotlib import pyplot
from datetime import datetime
import operator
from collections import defaultdict
mtable=[]
#Returns the average of the central 16 pixels for all 64 frames
def fluxofcentroid(image):
	return np.mean( image_data[:,14:18,14:18])	
#Updates mtable with flux of star and observation time
def updatetable(flux,obstime):
	#To make it as a datetime object
	obsdatetime=dateutil.parser.parse(obstime)
	mtable.append([flux,obsdatetime])
	return
def plotcurve(mtable,direc):
	mtable=sorted(mtable,key=lambda y: y[1])
	timstmp=[]
	#Each imagecube is taken after intervals of 128secs.
	for i in range(len(mtable)):
		timstmp.append(128*i)
	yaxis=[]
	yaxis =[x[0] for x in mtable] 
	plotlightcurve(timstmp,yaxis,direc)

def plotlightcurve(xaxis,yaxis,direc):									
	pyplot.plot(xaxis,yaxis,'ro')
	# X axis : Time of observation ( taken from FITS header ) in seconds
	# Y axis : Flux in MJ/st
	pyplot.title('Flux vs time lightcurve for AOR-r'+str(direc)) 
	pyplot.xlabel('Time of observation (s)')
	pyplot.ylabel('Flux (MJ/st)')
	pyplot.grid(True)
	pyplot.savefig('Flux vs time lightcurve for AOR-r'+str(direc)+'.pdf')
	pyplot.close()

ct = 0
#direc is the directory of images, the AOR number without 'r'
direc=20150013
path= "/home/hema/Documents/mcgill/XO-3_b_sim_ch2_150702/r"+str(direc)+"/ch2/bcd"
for filename in glob.glob(os.path.join(path, '*bcd.fits')):
	f=pyfits.open(filename)
	image_data=f[0].data
	header=f[0].header
	starflux=fluxofcentroid(image_data)
	updatetable(starflux,header['DATE_OBS'])
	ct +=1
plotcurve(mtable,direc)
