#--------------------------------------------------------------
#Author: Lisa Dang & Sudarsan S
#Created: 2016-07-28 09:02 PM IST
#Last Modified: 
#Title: Plot the fits image cube
#Test data:AOR r46470912- SPITZER_I2_46470912_0000_0000_2_bcd.fits
#--------------------------------------------------------------

import pyfits
from matplotlib.backends.backend_pdf import PdfPages
from pylab import *
import matplotlib.gridspec as gridspec
import numpy as np
import sys
import matplotlib.pyplot as plt

def plotimgcube(path,pdfname='dispimgcube'):
	#Open the FITS image
	try:
		f = pyfits.open(path,mode='readonly')
	except IOError:
		print("File does not exist")
		return
		
	image_data = f[0].data
	tam = image_data[0].shape
	print 'tam',tam
	#stacking 64 into 1 image by taking pixel wise median
	pdf =PdfPages(pdfname)

	stacked=np.nanmedian( image_data, axis=0)	
	plt.imshow(stacked,interpolation='none',cmap='gray')
	plt.colorbar()
	plt.title('Median image of 64 frames')
	plt.xlabel('x pixels')
	plt.ylabel('y pixels')
	plt.ylim([0, tam[0]])
	plt.xlim([0, tam[1]])
	pdf.savefig()
	#Plot image for each frame
	fig = plt.figure(figsize = (10,10))
	i=0
	gs = gridspec.GridSpec(8, 8)
	for frameno in np.arange(1,65):
		ax = plt.subplot(gs[i])
		ax.set_title('#'+ str(frameno),fontsize=5)
		ax.imshow(image_data[i], interpolation='none', cmap=plt.cm.gray)
		ax.set_aspect('auto')
		ax.axis('image')
		ax.set_aspect('equal')
		ax.set_xticks([])
		ax.set_yticks([])
		i+=1
	plt.tight_layout()
	pdf.savefig(dpi=500)

	# Set meta data for the pdf file
	d = pdf.infodict()
	d['Title'] = str()
	d['Author'] = u'Popeye'
	d['Subject'] = 'Median stacked image and frame wise plot for imagecube'
	pdf.close()


aorpath="/home/hema/Documents/mcgill/handy/aorkeys-20-selected_AORs/r46470912/ch2/bcd/"
filename="SPITZER_I2_46470912_0000_0000_2_bcd.fits"
path=aorpath+filename
plotimgcube(path,'r46470912_0.pdf')
