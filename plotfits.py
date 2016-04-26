# Plot FITS images and save them in PDF
import pylab
import pyfits
from matplotlib.backends.backend_pdf import PdfPages

pdfname = 'image_pdf.pdf'
with PdfPages(pdfname) as pdf:
	#Open the FITS image
	f = pyfits.open('SPITZER_I2_46467072_0000_0000_2_bcd.fits')
	image_data = f[0].data
	tam = image_data[0].shape
	#Plot image for each frame
	# if loop is removed,one frame of image will be plotted
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
