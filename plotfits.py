# Plot FITS images and save them in PDF
# Corrected to include mediam image of 64 frames
# Corrected to show images in pixelated form
import pylab
import pyfits
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

pdfname = 'imagespdf.pdf'
with PdfPages(pdfname) as pdf:
	#Open the FITS image
	f = pyfits.open('SPITZER_I2_20150000_0000_0000_1_bcd.fits')
	image_data = f[0].data
	tam = image_data[0].shape
	#stacking 64 into 1 image by taking pixel wise median
	stacked=np.median( image_data, axis=0)	
	pylab.clf()
	pylab.gray()
	pylab.imshow(stacked,interpolation='none',cmap='gray')
	pylab.title('Median image of 64 frames')
	pylab.xlabel('x pixels')
	pylab.ylabel('y pixels')
	pylab.ylim([0, tam[0]])
	pylab.xlim([0, tam[1]])
	pdf.savefig()
	#Plot image for each frame
	for frameno in range (1,64):
		plt.clf()
		plt.gray();
		plt.imshow(stacked, interpolation='none',cmap='gray')
		plt.xlabel('x pixels')
		plt.ylabel('y pixels')
		plt.title('frame number #'+ str(frameno)  )
		plt.ylim([0, tam[0]])
		plt.xlim([0, tam[1]])
		pdf.savefig()
