from scipy.stats import norm
import numpy as np
import math
import glob
import re
from scipy import integrate
import matplotlib.pylab as plt
import numpy.random
import pywcs

#Assume the star formation rates are given in the form of a probability distribution
def sim(num=100, ma=1, type='z15b', exf=1, SFR=1, mag_limit=30, A_v=1):

    #making the array that this returns
    length = num*16
    SNs = np.zeros((length, 5))

    #dealing with the individual RA an DECs for clusters
    hdulist = pyfits.open("abell383/weight.fits")
    ra_decdata = hdulist[0].data
    wcs = pywcs.WCS(hdulist[0].header)
    form = ra_decdata.shape
    max_rpix = form[0]-1
    max_dpix = form[1]-1
    pixcrd = np.array([[0,0], [max_rpix, max_dpix]])
    sky = wcs.wcs_pix2sky(pixcrd, 1)
    ra_max = sky[0,0]
    ra_min = sky[1,0]
    dec_max = sky[0,1]
    dec_min = sky[1,1]
    SNs[:,0] = np.random.uniform(ra_min, ra_max, length) #RA
    SNs[:,1] = np.random.uniform(dec_min, dec_max, length) #DEC

    #any supernovae with RA and DECs that fall within these values, but not in the ACS region get assigned a magnitude of 500. First part is matching the assigned supernovae ra and dec with the correct pixels. 
    pixcoordsdec = np.array([range(max_rpix), np.zeros(max_rpix)])
    pixcoordsdec = np.array([np.zeros(max_dpix), range(max_dpix)])
    ra_pix = pixcoordsr.transpose()
    dec_pix = pixcoordsdec.transpose()
    ra_sky = wcs.wcs_pix2sky(ra_pix, 1)
    dec_sky = wcs.wcs_pix2sky(dec_pix, 1)
    
    SNs[:,3] = np.random.uniform(0, 1000, size=length)
    for i in (np.arange(17)): #this is only going to make integer redshifts
        #SNs.append([0, 0, i, 0, 0])
        SNs[(i-1)*num:i*num,2] = i+4
    for q in range(length):
        idxr = (np.abs(ra_sky[:,0] - SNs[q,0])).argmin()
        idxd = (np.abs(dec_sky[:,1] - SNs[q,1])).argmin()
        if ra_decdata[idxr, idxd] == 0:
            SNs[q,4] = 500
        else:
            SNs[q,4] = SNmag(type, SNs[q,2], SNs[q,3])

     return SNs
