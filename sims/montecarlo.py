from scipy.stats import norm
import numpy as np
import math
import glob
import asciitable
import astropy.io.ascii as ascii
import re
from scipy import integrate
import matplotlib.pylab as plt
import numpy.random
import pywcs
import maglim
import os
import volfun
import pyfits
import SN_number
import time
import SN_mag #update the SN_number script that's uploaded onto github

#   FUTURE CHANGES TO IMPLEMENT
#
#   noise: use the 1-sigma error maps provided by Adi and the average noise around the areas to get a good idea of the signal-to-noise (more accurate...but not sure if it is actually worth it)
#   High redshift star formation rates: This definitely has the biggest effects on the simulated numbers, but not sure if we really know this well enough
#   Change the files so that it outputs the SNs in a text document and it is also capable of reading in the volumes from a text document
#   make sure we've accounted for the duration of the survey correctly (check to see what the duration of the survey actually is)
#   Also add in timestamps like Steve's code and headers that are actually well documented
#   What is a good way to code in the SFR besides just using a function that is in the SN_number.py program?
#   Also it's probably a good idea to write in the option to read in light curves, or make it like Leonidas does and read in a file containing the names of files that we want.
#   Find a way to print out some kind of useful information about the simulation so that we know approximately how far along we've gotten
#   Scale the lightcurves up correctly (instead of just adding a magnitude or two)
#   Check that the cosmology really is consistent between all the parts...I fear it is not consistent with the output volumes.


#   OUTPUT
#
#   Each column of the SN_array will contain (in order):
#   the SN type: currently just using the brightest representative SN from each of the three differen masses to get an upper limit, column 0
#   the RA and DEC coordinates: columns 1 & 2
#   the redshift
#   the magnification
#   the time (in days after the explosion)
#   the magnitude of the object at the two epochs of observation.

def mumag(mag, mu):
    jsky = np.power(10, 23 -(mag+ 48.60)/2.5)
    mu_jsky = mu*jsky
    mu_mag = 2.5*(23 - np.log10(mu_jsky)) - 48.60
    return mu_mag

def find_nearest(array, value):
    idx = (np.abs(array-value)).argmin()
    return array[idx]

def find_nearest_idx(array, value):
    idx = (np.abs(array-value)).argmin()
    return idx

def type_mass(num):
    if num <=2:
        return 15
    elif (num > 2) & (num <= 5):
        return 25
    else:
        return 40

#type='z15g'
#type='z40g'
#type = ['z15g', 'z25g', 'z40b']
#type = ['z15g', 'z25g', 'z40g']
#mag_lim=26.8
#style = 10 #if anything but 0, then the program runs via the whalen method where you pick a minimum magnification (by inputting that number) and then calculate. If 0=, then the program runs by its usual method--it calculates the necessary minimum magnification
#cadence = 90 #in days
#   detectlim: if 0, then the detection limits are off (everything that is created is detected). If 1, then the detection limits are on (cadence & mag_lim become important)
#   once the lightcurves are fixed you can just have it read in light curves...
def simulation(zmin=5, zmax=12, type='z15g', mag_lim=26.8, style=10, cadence=90, detectlim=0):
    start = time.time()
    clusters = clist['clusterdir']
    
    types = ['z15b', 'z15d', 'z15g', 'z25b', 'z25d', 'z25g', 'z40b', 'z40g']
    type_num = types.index(type)
    mass = type_mass(type_num)
    #if type_num <= 2:
    #    mass = 15
    #elif (type_num > 2) & (type_num <=5):
    #    mass = 25
    #else:
    #    mass = 40

    #   Well, now we are always going to run the simulation for one representative type of each mass range, so we won't need to have the user input a specific type anymore. The only exception is that we might need to do two iterations (one for 'z40B' and one for 'z40G'

    spacing = zmax - zmin + 1
    ztargets = np.linspace(zmin, zmax, spacing)
    mag_array = np.zeros((1,2))

    #   Now we've changed it so that if you pick style = 0, that means that it will run the usual way. Assigning any other number to style causes it to pick that as the minimum magnification.
    if style == 0:
        for x in ztargets:
            l_c = maglim.lightcurve(type, x)
            mu = maglim.magfactor(l_c, mag_lim)
            mag_array = np.r_[mag_array, [[x, mu]]]
            print "[[x, mu]]", x, mu
        mag_array = mag_array[1:len(ztargets)+1, :]
        print "mag_array", mag_array
    else:
        mag_array = np.c_[ztargets, np.zeros(len(ztargets))+style]
        print "mag_array", mag_array

    SN_array = np.zeros((1,9))
    number_arr = np.zeros(0)
    for x in range(len(ztargets)):
        number_arr = SN_number.number_SN(ztargets[x], mass, mu_min = mag_array[x, 1])
        #this calculates the number at each magnitude, but note that it won't be integer numbers at this point

        nc = len(clusters)
        
        for y in np.linspace(1, nc, nc):
            number = np.c_[number_arr[0:len(number_arr)-1,0], np.abs(np.diff(number_arr[:,y]))]
            #now we turn the numbers into integer numbers to put into the rest of the code
            probs = number[:,1] - np.fix(number[:,1])
            number[:,1] = np.fix(number[:,1])
            for i in range(len(probs)):
                k = np.random.uniform(0, 1)
                if probs[i] > k:
                    number[i, 1] += 1
            #now we have integer numbers, though they'll change each time, so now we are going to make that many supernovae at that redshift with those magnifications
            for i in range(len(number)):
                SNs = np.zeros((number[i,1], 9))
                if SNs.size > 0:
                    SNs[:,0] = type_num #type of SN
                    SNs[:,3] = ztargets[x] #redshift
                    SNs[:,4] = number[i,0] #magnification
                    #SNs[:,5] = np.random.uniform(1, 1000) #days after explosion --> need to not assign the same day to every one of them! do this with the RA and DEC part
                    SNs[:,8] = y-1 #which cluster it is by number
                    hdulist = pyfits.open("../{}/weight.fits".format(clusters[y-1]))
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
                    while n < len(SNs):
                        RA = np.random.uniform(ra_min, ra_max)
                        DEC = np.random.uniform(dec_min, dec_max)
                        idxr = (np.abs(ra_sky[:,0] - RA)).argmin()
                        idxd = (np.abs(dec_sky[:,1] - DEC)).argmin()
                        if ra_decdata[idxr, idxd] == 1:
                            SNs[n, 1] = RA
                            SNs[n, 2] = DEC
                            SNs[n, 5] = np.random.uniform(1, 500)
                            n += 1

                SN_array = np.r_[SN_array, SNs]

        if detectlim == 0:
            num_detection = len(SN_array)
            detected_array = SN_array
            print("ran if statement")
        elif (detectlim != 0) & (os.path.isfile('../output/lightcurve_{}_{}.dat'.format(type, ztargets[x]))):
            start3 = time.time()
            num_detection = 0
            detected_array = np.zeros((1, 9))
            lightcurve = '../output/lightcurve_{}_{}.dat'.format(type, ztargets[x])
            for n in (np.arange(len(SN_array)-1) + 1):
                day = SN_array[n, 5]
                mu = SN_array[n, 4]
                redshift = SN_array[n, 3]
                dat = ascii.read(lightcurve)
                dat = np.asarray(dat)
                lc_dat = dat.view(np.float).reshape(dat.shape + (-1,))
                idx = find_nearest_idx(lc_dat[:,0], day)
                idx2 = find_nearest_idx(lc_dat[:,0], day + cadence)
                mag_1 = lc_dat[idx, 1] #unmagnified magnitude
                mag_1 = mumag(mag_1, mu)
                mag_2 = lc_dat[idx2, 1]
                mag_2 = mumag(mag_2, mu)
                SN_array[n, 6] = mag_1
                SN_array[n, 7] = mag_2
                flux_1 = np.power(10, 23 -(mag_1 + 48.60)/2.5)
                flux_2 = np.power(10, 23 -(mag_2 + 48.60)/2.5)
                flux_diff = np.abs(flux_1 - flux_2)
                flux_lim = np.power(10, 23 -(mag_lim + 48.60)/2.5)
                if flux_diff > flux_lim:
                    num_detection += 1
                    #print SN_array[n]
                    #print detected_array
                    detected_array = np.r_[detected_array, [SN_array[n]]]
            end3 = time.time()
            print("Completed one elif statement in %i seconds"%(end3-start3))
        #   change to suit the format of the lightcurve THIS PART IS INCOMPLETE. For now we will assume that the data is given in a format with days in one column (the first column) and the magnitudes in another column (the second column)
        else:
            num_detection = 0
            detected_array = np.zeros((1, 9))
            start2 = time.time()
            for n in (np.arange(len(SN_array)-1) + 1):
                day = SN_array[n,5]
                mag_1 = SN_mag.SNmag(type, SN_array[n,3], day) - 0.5
                mag_1 = mumag(mag_1, SN_array[n,4])
                mag_2 = SN_mag.SNmag(type, SN_array[n,3], day + cadence) - 0.5
                mag_2 = mumag(mag_2, SN_array[n,4])
                SN_array[n, 6] = mag_1
                SN_array[n, 7] = mag_2
                flux_1 = np.power(10, 23 -(mag_1 + 48.60)/2.5)
                flux_2 = np.power(10, 23 -(mag_2 + 48.60)/2.5)
                flux_diff = np.abs(flux_1 - flux_2)
                flux_lim = np.power(10, 23 -(mag_lim + 48.60)/2.5)
                if flux_diff > flux_lim:
                    num_detection += 1
                    #print SN_array[n]
                    detected_array = np.r_[detected_array, [SN_array[n]]]
            end2 = time.time()
            print("completed one else statement in %i seconds."%(end2-start2))

#   determining the number of detections. if detectlim = 0, then there are no detection limits and everything that is created is detected. for all others, we use the cadence and mag_lim to determine if it is a detection.

#   change this so that if some files are available of the lightcurves (or given I suppose) then it won't go to SN_mag, it will just directly look at the lightcurve files. Note that the lightcurve files will have to be files that are at the correct redshift for it to be faster
    SN_array = np.delete(SN_array, 0, 0)
    detected_array = np.delete(detected_array, 0, 0)
    print "number of detections", num_detection, np.count_nonzero(detected_array.size) != 0
    end = time.time()
    print("One run in %i seconds."%(end-start))
    if np.count_nonzero(detected_array) != 0:
        return num_detection, SN_array, detected_array
    else:
        return num_detection, SN_array, np.zeros((1,9))

#   now write something that runs this simulation many times.

lensmodels = '../data/vp_lensmodels.lis'
#   clist=asciitable.read(lensmodels, names=['alias','clusterdir','deflxfile','deflyfile','segfile','zcluster','zmodel'])
clist=asciitable.read(lensmodels, names=['alias','clusterdir','deflxfile','deflyfile','zcluster','zmodel'],guess=False)

clusters = clist['clusterdir']

runs=100
#supernovae = ['z15g', 'z40g']
supernovae = ['z15g', 'z25g', 'z40g']
#supernovae = ['z15g', 'z40g', 'z25g']
# supernovae = ['z15g', 'z25g', 'z40g']
detections = []
supernovae_info = np.zeros((1, 9))
detect_info = np.zeros((1, 9))
i = 0
while i < runs:
    totaldetect = 0
    totalinfo = np.zeros((1, 9))
    detectinfo = np.zeros((1, 9))
    for s in supernovae:
        outdat = simulation(zmin=5, zmax=12, type=s, detectlim=1)
        totaldetect = totaldetect + outdat[0]
        totalinfo = np.r_[totalinfo, outdat[1]]
        if np.count_nonzero(outdat[2]) != 0:
            detectinfo = np.r_[detectinfo, outdat[2]]
    totalinfo = np.delete(totalinfo, 0, 0)
    detectinfo = np.delete(detectinfo, 0, 0)
    detections.append(totaldetect)
    supernovae_info = np.r_[supernovae_info, totalinfo]
    detect_info = np.r_[detect_info, detectinfo]
    print("Finished run %i. %i more runs left."%((i+1), (runs-i-1)))
    i += 1
supernovae_info = np.delete(supernovae_info, 0, 0)
detect_info = np.delete(detect_info, 0, 0)
#write the output. First need to convert the supernovae info into a record array with the appropirate SN type and cluster in there
types = ['z15b', 'z15d', 'z15g', 'z25b', 'z25d', 'z25g', 'z40b', 'z40g']
rec_array = np.recarray((len(supernovae_info),), dtype=[('SN type', str, 4), ('RA', float), ('DEC', float), ('z', float), ('mu', float), ('epoch_days', float), ('ABmag_1', float), ('ABmag_2', float), ('cluster', str, 16)])
rec_array['SN type'] = [types[int(i)] for i in supernovae_info[:,0]]
rec_array['RA'] = supernovae_info[:,1]
rec_array['DEC'] = supernovae_info[:,2]
rec_array['z'] = supernovae_info[:,3]
rec_array['mu'] = supernovae_info[:,4]
rec_array['epoch_days'] = supernovae_info[:,5]
rec_array['ABmag_1'] = supernovae_info[:,6]
rec_array['ABmag_2'] = supernovae_info[:,7]
rec_array['cluster'] = [clusters[int(i)] for i in supernovae_info[:,8]]
if np.count_nonzero(detect_info) != 0:
    dect_array = np.recarray((len(detect_info),), dtype=[('SN type', str, 4), ('RA', float), ('DEC', float), ('z', float), ('mu', float), ('epoch_days', float), ('ABmag_1', float), ('ABmag_2', float), ('cluster', str, 16)])
    dect_array['SN type'] = [types[int(i)] for i in detect_info[:,0]]
    dect_array['RA'] = detect_info[:,1]
    dect_array['DEC'] = detect_info[:,2]
    dect_array['z'] = detect_info[:,3]
    dect_array['mu'] = detect_info[:,4]
    dect_array['epoch_days'] = detect_info[:,5]
    dect_array['ABmag_1'] = detect_info[:,6]
    dect_array['ABmag_2'] = detect_info[:,7]
    dect_array['cluster'] = [clusters[int(i)] for i in detect_info[:,8]]
    ascii.write(dect_array, '../output/mc_rftf_hsfrdz_{}.dat'.format(runs))
#colnames = ['SN type', 'RA', 'DEC', 'z', 'mu', 'epoch_days', 'ABmag_1', 'ABmag_2', 'cluster']
#ascii.write(rec_array, '../output/montecarlo_test_{}.dat'.format(runs))
ascii.write(rec_array, '../output/mc_rftf_hsfrdz_total_{}.dat'.format(runs))
plt.clf()
plt.hist(detections)
plt.savefig("../output/mac_rftf_lsfrdz_{}.png".format(runs))
#plt.savefig("mc_test.png")
plt.show()



