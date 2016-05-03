# Plot FITS images and save them in PDF
# Plot FITS images and save them in PDF
# Corrected to include mediam image of 64 frames
import pylab
import pyfits
from matplotlib.backends.backend_pdf import PdfPages
import pylab
import pyfits
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

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
	pylab.imshow(stacked)
	pylab.title('Median image of 64 frames')
	pylab.xlabel('x pixels')
	pylab.ylabel('y pixels')
	pylab.ylim([0, tam[0]])
	pylab.xlim([0, tam[1]])
	pdf.savefig()
	#Plot image for each frame
	for frameno in range (1,64):
		pylab.clf()
		pylab.gray();
		pylab.imshow(image_data[frameno])
		pylab.xlabel('x pixels')
		pylab.ylabel('y pixels')
		pylab.title('frame number #'+ str(frameno)  )
		pylab.ylim([0, tam[0]])
		pylab.xlim([0, tam[1]])
		pdf.savefig()
