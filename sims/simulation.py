from scipy.stats import norm
import numpy as np
import math
import glob
import re
from scipy import integrate
import matplotlib.pylab as plt
import numpy.random
import pywcs
import mag_lim
import volfun
import SN_number #update the SN_number script that's uploaded onto github
#import SN_mag

zmin=5
zmax=6
type='z15g'
mag_lim=27
cluster= 'abell383'

def mumag(mag, mu):
    jsky = np.power(10, 23 -(mag+ 48.60)/2.5)
    mu_jsky = mu*jsky
    mu_mag = 2.5*(23 - np.log10(mu_jsky)) - 48.60
    return mu_mag

#changed this for doing this for one cluster at a time

#Since the elements of a np array all have to be the same, I will have to create an return a structured array if I want to be able to include information on the different types of SNs.
#Actually I decided to just get rid of the type. I can run it for a certain type one at a time. Or...I can have the type be a float and then map the float onto a string type at the end so I don't have to keep remaking record arrays
types = ['z15b', 'z15d', 'z15g', 'z25b', 'z25d', 'z25g', 'z40b', 'z40g']
type_num = types.index(type)

#seems like what you want to do is sort of dynamically create the array. Perhaps the array should be created given two criteria (the random values fall within the ACS regions) and the correct number for reach redshift based on Leonidas' code.

#Each column of the SN_array will contain the SN type, the RA and DEC coordinates, the redshift, the magnification, the time (in days after the explosion) and the magnitude of the object--> actually maybe I'll have two time columns. The first for the first observation and the second for the second observation. The second one will be two weeks after the first one, but we should adjust it so that the time after the first one varies.

#First pick an SN type and calculate the necessary magnification at various redshifts. First column is redshift, second column is the magnification necessary to see it. Needs mag_lim.py to run properly

spacing = zmax - zmin + 1
ztargets = np.linspace(zmin, zmax, spacing)
mag_array = np.zeros((1,2))
for x in ztargets:
    l_c = lightcurve(type, x)
    mu = mag_factor(l_c, mag_limit)
    mag_array = np.r_[mag_array, [[x, mu]]]
mag_array = mag_array[1:len(ztargets)+1, :]
#print mag_array.shape
print "mag_array", mag_array

#This first part needs the SN_number.py code to run. This is where we calculate the number of SNs that explode for each z and mu range. The first column of the SN_number output is the magnification and the second column is the number of SNs expected for the SFR and eta that were coded into the SN_number function.

#For the number, a lot of them are less than 1, so you can try doing somewhere you do a random draw and if the random draw is less than the number that you have, you get a SN if it is more, you don't. Do this first then feed the number array into the RA and DEC assigner
SN_array = np.zeros((1,8))
#couldn't figure out how to add them without adding in this useless first row...just take it out at the end I guess.
number_arr = np.zeros(0)
for x in range(len(ztargets)):
    #print ztargets[x]
    #print mag_array[x, 1]
    number_arr = number_SN(ztargets[x], mu_min = mag_array[x, 1])
    #this calculates the number at each magnitude, but note that it won't be integer numbers at this point
    number = np.c_[number_arr[0:len(number_arr)-1,0], np.abs(np.diff(number_arr[:,1]))]
    #now we turn the numbers into integer numbers to put into the rest of the code
    probs = number[:,1] - np.fix(number[:,1])
    #print "probs", probs
    number[:,1] = np.fix(number[:,1])
    #print "number", number
    for i in range(len(probs)):
        k = np.random.uniform(0, 1)
        if probs[i] > k:
            number[i, 1] += 1
    #now we have integer numbers, though they'll change each time, so now we are going to make that many supernovae at that redshift with those magnifications.
    for i in range(len(number)):
        SNs = np.zeros((number[i,1], 8))
        SNs[:,0] = type_num #type of SN
        SNs[:,3] = ztargets[x] #redshift
        SNs[:,4] = number[i,0] #magnification
        SNs[:,5] = np.random.uniform(1, 1000) #days after explosion
        SN_array = np.r_[SN_array, SNs]
SN_array = np.delete(SN_array, 0, 0)

#now to assign the correct RA and DEC coordinates
hdulist = pyfits.open("{}/weight.fits".format(cluster))
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

pixcoordsra = np.array([range(max_rpix), np.zeros(max_rpix)])
pixcoordsdec = np.array([np.zeros(max_dpix), range(max_dpix)])
ra_pix = pixcoordsra.transpose()
dec_pix = pixcoordsdec.transpose()
ra_sky = wcs.wcs_pix2sky(ra_pix, 1)
dec_sky = wcs.wcs_pix2sky(dec_pix, 1)

n = 0
while n < len(SN_array):
    RA = np.random.uniform(ra_min, ra_max)
    DEC = np.random.uniform(dec_min, dec_max)
    idxr = (np.abs(ra_sky[:,0] - RA)).argmin()
    idxd = (np.abs(dec_sky[:,1] - DEC)).argmin()
    if ra_decdata[idxr, idxd] == 1:
        SN_array[n, 1] = RA
        SN_array[n, 2] = DEC
        n += 1

#calculates the magnified magnitude this needs SN_mag.py and the SNmag function
#one way to make everything faster is to just calculate the lightcurve for that particular explosion type once at the beginning and then just reference that.

num_detection = 0
for n in range(len(SN_array)):
    day = SN_array[n,5]
    mag_1 = SNmag(type, SN_array[n,3], day)
    mag_1 = mumag(mag_1, SN_array[n,4])
    mag_2 = SNmag(type, SN_array[n,3], day + 14)
    mag_2 = mumag(mag_2, SN_array[n,4])
    SN_array[n, 6] = mag_1
    SN_array[n, 7] = mag_2
    if (SN_array[n,6] < 27) & (SN_array[n,7]<27):
        num_detection += 1

#if both of the last two are above the cutoff count it as a detection
print "number of detections", num_detection

#now we will pick two random times after the explosion (record the first time in the array) and then calculate the magnitudes at both times.

#turn it into a record array in the last step so that all of the columns and stuff will be labeled.
