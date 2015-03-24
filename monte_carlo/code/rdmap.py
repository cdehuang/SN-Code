#
#   chuang
#   edited--version from 03.24.2015
#

import pyfits
import numpy as np
import os.path

"""
   Makes a FITS file covering the ACS + WFC3 region of each cluster. The FITS file is 1 in the cluster region and 0 outside of the cluster region. This is just so that I make sure all of the RA and DEC points that the monte carlo generates are in the right location (the points get multiplied by the map and ones that aren't right get redistributed)
   output FITS files are approximately 100 MB in size...
"""

#all of the CLASH clusters (25) (need directories with these names)

clusters = ['abell1423', 'abell209', 'abell2261', 'abell611', 'macs0329', 'macs0416','macs0429', 'macs0647', 'macs0717', 'macs0744', 'macs1115', 'macs1149', 'macs1206', 'macs1311', 'macs1423', 'macs1720', 'macs1931', 'macs2129', 'ms2137', 'clj1226', 'rxj1347', 'rxj1532', 'rxj2129', 'rxj2248']

#all of the mosaics to use in making the fits files (I just grabbed the ones that were most recent at the time of writing).

fitsfiles = ['a1423_mosaic_065mas_acs_wfc3ir_total_drz_20130227.fits.gz', 'a209_mosaic_065mas_acs_wfc3ir_total_drz_20130119.fits.gz', 'a2261_mosaic_065mas_acs_wfc3ir_total_drz_20130123.fits.gz', 'a611_mosaic_065mas_acs_wfc3ir_total_drz_20120618.fits.gz', 'clj1226_mosaic_065mas_acs_f775w_drz_20130626.fits', 'macs0329_mosaic_065mas_acs_wfc3ir_total_drz_20111105.fits.gz', 'macs0416_mosaic_065mas_acs_wfc3ir_total_drz_20130119.fits.gz', 'macs0429_mosaic_065mas_acs_wfc3ir_total_drz_20130227.fits.gz', 'macs0647_mosaic_065mas_acs_wfc3ir_total_drz_20111202.fits.gz', 'macs0717_mosaic_065mas_acs_wfc3ir_total_drz_20111230.fits.gz', 'macs0744_mosaic_065mas_acs_wfc3ir_total_drz_20111230.fits.gz', 'macs1115_mosaic_065mas_acs_wfc3ir_total_drz_20120305.fits.gz', 'macs1149_mosaic_65mas_acs_wfc3ir_total_drz_20110825.fits.gz', 'macs1206_mosaic_065mas_acs_wfc3ir_total_drz_20110815.fits.gz', 'macs1311_mosaic_065mas_acs_wfc3ir_total_drz_20130709.fits', 'macs1423_mosaic_065mas_acs_wfc3ir_total_drz_20130313.fits.gz', 'macs1720_mosaic_065mas_acs_wfc3ir_total_drz_20120628.fits.gz', 'macs1931_mosaic_065mas_acs_wfc3ir_total_drz_20120628.fits.gz', 'macs2129_mosaic_065mas_acs_wfc3ir_total_drz_20120122.fits.gz', 'ms2137_mosaic_065mas_acs_f775w_drz_20130121.fits.gz', 'rxj1347_mosaic_065mas_acs_wfc3ir_total_drz_20110812.fits.gz', 'rxj1532_mosaic_065mas_acs_wfc3ir_total_drz_20120416.fits.gz', 'rxj2129_mosaic_065mas_acs_wfc3ir_total_drz_20121016.fits.gz', 'rxj2248_mosaic_065mas_acs_wfc3ir_total_drz_20130119.fits.gz']

clusterdict = dict(zip(clusters, fitsfiles))

for i,j in clusterdict.items():
    if os.path.isfile('../{}/weight.fits'.format(i)) is False:
        hdulist = pyfits.open('../{}/{}'.format(i,j))
        imdata = hdulist[0].data

        nonzero = np.nonzero(imdata)
        imdata[nonzero[0], nonzero[1]] = 1
        hdulist.writeto('../{}/weight.fits'.format(i))


#open the image FITS file
#hdulist = pyfits.open("abell383/a383_mosaic_065mas_acs_wfc3ir_total_drz_20110825.fits.gz")
#imdata = hdulist[0].data

#to get the shape
#imdata.shape
#imdata.dtype.shape

#change all of the nonzero pixel values to 1
