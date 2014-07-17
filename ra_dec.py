import pyfits
import numpy as np

#try with just this for now, and simplify everything later (currently takes a few seconds to run)

#open the image FITS file
hdulist = pyfits.open("abell383/a383_mosaic_065mas_acs_wfc3ir_total_drz_20110825.fits.gz")
imdata = hdulist[0].data

#change all of the nonzero pixel values to 1
nonzero = np.nonzero(imdata)
imdata[nonzero[0], nonzero[1]] = 1
hdulist.writeto("abell383/weight.fits")
#can't save changes back to original fits file unless you unzip it, so just do this. But you can just go through and make the weights for everything. 
