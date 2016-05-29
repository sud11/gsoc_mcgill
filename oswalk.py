# Code snippet to read all the FITS files in the directories and subdirectories of the root and store the file names in a csv file
# Code tested in Ubuntu 15.04. Some tweaks might be required to make it work in Windows OS
import os
import glob
import sys
import csv
fnames=[]
outcsv = open("oswalk_filenames.csv",'wb')
writer = csv.writer(outcsv, dialect='excel',delimiter='\n')
writer.writerow(['File_name'])
root = '/home/hema/Documents/mcgill/handy'
for path, subdirs, files in os.walk(root):
   	fnames.extend(  [filename for filename in glob.glob(os.path.join(path, '*bcd.fits'))])
for i in range( len(fnames)):
	fnames[i]= fnames[i][fnames[i].find('SPITZER'):]
fnames.sort()
writer.writerow(fnames)
