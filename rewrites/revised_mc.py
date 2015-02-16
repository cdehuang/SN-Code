#   chuang
#   01.27.15
#
#   Finally just decided to rewrite the original Monte Carlo code to fix the previous issues.
#   New one takes in a cadence, filter, limiting magnitude, and cluster name.
#

#   python packages to import

import numpy as np
import math
import glob
import asciitable
import astropy.io.ascii as ascii
import astropysics.phot
import pyfits
import time
import numpy.random
import os
import astropy.wcs as pywcs

#   modules I wrote that are necessary

import new_volfun
import new_SN_number
import new_SN_mag
import new_maglim

#   Takes in a magnitude and a magnification and returns the new magnitude after the flux has been magnified

def mumag(mag, mu):
    jsky = np.power(10, 23 -(mag+ 48.60)/2.5)
    mu_jsky = mu*jsky
    mu_mag = 2.5*(23 - np.log10(mu_jsky)) - 48.60
    return mu_mag

#   Find the value of the array that is nearest to the input value

def find_nearest(array, value):
    idx = (np.abs(array-value)).argmin()
    return array[idx]

#   Find the index of the array that is closest to the input value

def find_nearest_idx(array, value):
    idx = (np.abs(array-value)).argmin()
    return idx

#   Tells you how massive the supernovae is based on the type number (only useful for the Whalen stuff)

def type_mass(num):
    if num <=2:
        return 15
    elif (num > 2) & (num <= 5):
        return 25
    else:
        return 40

#   new class since record arrays can't change their size. But it's not that useful without tons of changes, so I think I won't use it (I'll just use the usual rec array)

class DynamicRecArray(object):
    def __init__(self, dtype):
        self.dtype = np.dtype(dtype)
        self.length = 0
        self.size = 10
        self._data = np.empty(self.size, dtype=self.dtype)
    
    def __len__(self):
        return self.length
    
    def append(self, rec):
        if self.length == self.size:
            self.size = int(1.5*self.size)
            self._data = np.resize(self._data, self.size)
        self._data[self.length] = rec
        self.length += 1
    
    def extend(self, recs):
        for rec in recs:
            self.append(rec)
    
    @property
    def data(self):
        return self._data[:self.length]

####
#
#   PARAMETERS
#
#   cluster name: name of the cluster (fully written out, no spaces, like abell2261--I don't think it matters, but noting is capitalized)
#   cadence: how long the cluster was observed (effectively)
#   filter: the reddest filter, written like F160W
#   magnitude limit: what was the limiting magnitude of the observations (note that this is differential)
#   zmin: the minimum redshift of SNe created
#   zmax: the maximum redshift of SNe created
#   detect_lim: if True, then the detection limits do matter, if False then the number of SNe is unaffected by the detection limits (i.e. everything created is detected)
#   type: you can choose the type of SNe you want to simulate
#   probably in the future this can be made so that the SNe light curve is an input parameter, or maybe it just reads in a file with all of the information...not sure which is easier
#
####

####
#
#   OUTPUT
#
#   An array of supernovae. The columns are, in order (and are floats unless otherwise specified):
#
#   0. The type of supernova explosion (string)
#   1. RA
#   2. DEC
#   3. redshift
#   4. magnification
#   5. date of first observation (relative to the explosion date, int)
#   6. apparent magnitude at 1st observation
#   7. apparent magnitude at 2nd observation
#   8. magnitude difference between the observations
#   9. the cluster in which it was located (string)
#
####

#   notes: looks like the magnitudes are getting assigned into the array in correctly. they seem to be right when they're being calculated though

#def simulation(parameter_file, filter, mag_lim, zmin=5, zmax=12, detect_lim=True, type='z15g', mass=15., observations_file='../../all_hst_surveys.csv', sfr='high', eta=1, IMF='Kroupa', maxvisibility=700, volume_file_dir='../../output/', volume_files_names='volumemag_whalen'):

#   these will be input parameters later. This code should probably read off the parameter file as well.
"""
type = 'z15g'
detect_lim=True
zmin=5
zmax=6
#mag_lim = 26.8
mag_lim=26.8
#filter='F140W'
filter='F160W'
cadence=87 #this should not be input, but should rather be a parameter from the observations_file
#   600 is the standard, but we probably want to add in something to just calculate this.
maxvisibility=700
#bandwidth = 3840
bandwidth=2683
mass = 15.
sfr='high'
IMF='Kroupa'

#   this parameter file is the list of the clusters, and should be the one that also gets fed into new_volfun so that they agree. It will be given to vol_fun.py through SN_number.py so we have to edit SN_number.py too.
parameter_file = '../../Debugging Things/vp_lensmodels.lis'

#   we may also want a second parameter file that has lists of the cadence, filters, and magnitude limits
observations_file = '../../output/file_26.8.dat'
observations = asciitable.read(observations_file)
cadences = observations['effective_cadence']
search = observations['note']
clusternames = observations['cluster']
"""
def simulate(configfile, mag_array, sn_type):
    c = {}
    execfile(configfile, c)
                                 
    #   start of the main code
    #   maybe pass in the magnification array and the parameter file
    observations_file = c['observations']
    
    lensmodels = c['lensmodels']
    bandwidth = c['bandwidth']
    observations_file = c['observations']
    observations = asciitable.read(observations_file)
    
    clist = asciitable.read(lensmodels, names=['alias', 'clusterdir', 'deflxfile', 'deflyfile', 'zcluster', 'zmodel'])
    
    clusters = clist['alias']
    
    start = time.time()
    
    spacing = zmax - zmin + 1
    ztargets = np.linspace(zmin, zmax, spacing)
    
    #   get the filter name
    
    filter_name = glob.glob('../../filters/WFC*{}*.dat'.format(filter))[0]
    SN_array = np.recarray((0,), dtype=[('sn_type', '|S10'), ('ra', float), ('dec',float), ('redshift', float), ('magnification', float), ('epoch', int), ('apparent_mag_1', float), ('apparent_mag_2', float), ('mag_difference', float) , ('cluster', '|S10'), ('cadence', int), ('note', '|S16')])
    
    for x in range(len(ztargets)):
        #   this calculates the number of SN we expect at each magnification (or greater), but note that it won't be integer numbers at this point.
        number_arr = new_SN_number.number_SN(ztargets[x], mass, configfile, mu_min=mag_array[x,1])
        just_numbers = new_SN_number.recarr_to_nparr(number_arr) #things are working up to here at least
        
        nc = len(observations)
        
        #   before we were iterating over clusters...now I think it would be a better idea to iterate over the searches...or rather each column in the file. There doesn't seem to be a good way to do it otherwise.
        
        for y in np.linspace(1, nc, nc):
            
            clustername = observations.dtype.names[int(y)]
            
            #   turns the number at that magnification or greater into something that's just the number at that magnification.
            number = np.c_[just_numbers[0:len(just_numbers)-1,0], np.abs(np.diff(just_numbers[:,y]))]
            
            #   the fractional SN are treated as probabilities.
            probs = number[:,1] - np.fix(number[:,1])
            
            #   now we turn the numbers into integer numbers to put into the rest of the code
            number[:,1] = np.fix(number[:,1])
            
            #   The Monte Carlo part: draw from a random uniform distribution between 0 and 1. If the drawn number is less than probs, then you get a SN, if it is greater then you don't.
            for i in range(len(probs)):
                k = np.random.uniform(0, 1)
                
                if probs[i] > k:
                    number[i, 1] += 1
            
            #   we have integer numbers, so now we are going to make that many supernovae at that redshift with those magnifications, and we can automatically assign them their type, redshift, cluster number, and magnification now.
            for i in range(len(number)):
                SNs = np.recarray((number[i,1],), dtype=[('sn_type', '|S10'), ('ra', float), ('dec',float), ('redshift', float), ('magnification', float), ('epoch', int), ('apparent_mag_1', float), ('apparent_mag_2', float), ('mag_difference', float) , ('cluster', '|S10'), ('cadence', int), ('search', '|S16')])
                if SNs.size > 0:
                    SNs['sn_type'] = type #mass of the SN
                    SNs['redshift'] = ztargets[x] #redshift
                    SNs['magnification'] = number[i,0] #magnification
                    SNs['cluster'] = clustername #which cluster it is by name
                    
                    #SNs[:,5] = np.random.uniform(1, 1000) #days after explosion --> need to not assign the same day to every one of them! do this with the RA and DEC part
                    
                    #   now assigning the RA and DEC of each supernova using the regions files. Weight.fits is a map where pixels in the cluster region get 1 and outside get 0.
                    hdulist = pyfits.open("../../{}/weight.fits".format(clusters[y-1]))
                    ra_decdata = hdulist[0].data
                    #   print hdulist[0].header
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
                            SNs[n][1] = RA
                            SNs[n][2] = DEC
                            
                            #   here it is assigning each supernova a random date of observation based on the observation window. Earliest date should be the number of days before a cadence. The latest should be after the brightest supernova with longest plateau/brightest magnification/highest redshift can be seen.
                            SNs[n][5] = np.random.uniform(-1*(cadence), maxvisibility)
                            n += 1
                # join the new SN with the master list of SN already created.
                SN_array = np.append(SN_array, SNs)
        
        #   Detections
        #
        #   basically here everything that is created is detected if we put in no detection limits
        if detect_lim == False:
            num_detection = len(SN_array)
            detected_array = np.recarray((0,), dtype=[('sn_type', '|S10'), ('ra', float), ('dec',float), ('redshift', float), ('magnification', float), ('epoch', int), ('apparent_mag_1', float), ('apparent_mag_2', float), ('mag_difference', float) , ('cluster', '|S10'), ('cadence', int), ('search', '|S16')])
            print("ran if statement")
        
        #   if we do have detection limits, then the program will first try to see if there is an output lightcurve file (since this saves a lot of time)
        elif (detect_lim == True) & (os.path.isfile('../../output/lightcurve_{}_{}_{}.dat'.format(type, ztargets[x], filter))):
            start3 = time.time()
            num_detection = 0
            detected_array = np.recarray((0,), dtype=[('sn_type', '|S10'), ('ra', float), ('dec',float), ('redshift', float), ('magnification', float), ('epoch', int), ('apparent_mag_1', float), ('apparent_mag_2', float), ('mag_difference', float) , ('cluster', '|S10'), ('cadence', int), ('search', 'S16')])
            lightcurve = '../../output/lightcurve_{}_{}_{}.dat'.format(type, ztargets[x], filter)
            
            #   now use the light curves to calculate the magnitude of the supernova at the date of observation by seeing how bright it is for a light curve at that day, and then magnifying it with mu
            for n in (np.arange(len(SN_array)-1) + 1):
                day = SN_array[n][5]
                mu = SN_array[n][4]
                #   print( len(SN_array))
                redshift = SN_array[n][3]
                dat = ascii.read(lightcurve)
                dat = np.asarray(dat)
                lc_dat = dat.view(np.float).reshape(dat.shape + (-1,))
                
                #   finding the magnitudes at date of first observation and date of follow-up
                idx = find_nearest_idx(lc_dat[:,0], day)
                idx2 = find_nearest_idx(lc_dat[:,0], day + cadence)
                
                #   just picked a large number to indicate that the SN hadn't occurred yet, since a matching index of zero indicates that it's before it exploded. Don't bother magnifying these, so the original and magnified magnitudes are both set to 500.
                if idx == 0:
                    mag_1 = 500.
                    mag_11 = 500.
                else:
                    mag_1 = lc_dat[idx, 1] #unmagnified magnitude
                    #print "unmagnified mag_1, day, idx, and mu", mag_1, day, idx, mu
                    mag_11 = mumag(mag_1, mu)
                    print mag_1, mag_11, mu, n
                #time.sleep(1)
                #print "magnified mag_1", mag_1
                if idx2 ==0:
                    mag_2 = 500.
                    mag_22 = 500.
                else:
                    mag_2 = lc_dat[idx2, 1]
                    mag_22 = mumag(mag_2, mu)
                    #time.sleep(1)
                    print mag_2, mag_22, mu, n
                
                #   enter the final magnified magnitudes into the array of SNe
                SN_array[n][6] = mag_11
                SN_array[n][7] = mag_22
                
                #   calculate the difference between the two observations to see if it met the magnitude limit
                flux_1 = astropysics.phot.mag_to_lum(mag_11)
                flux_2 = astropysics.phot.mag_to_lum(mag_22)
                flux_diff = np.abs(flux_1 - flux_2)
                mag_diff = astropysics.phot.lum_to_mag(flux_diff)
                SN_array[n][8] = mag_diff
                flux_lim = astropysics.phot.mag_to_lum(mag_lim)
                
                #   if greater than threshold, count as detection
                if flux_diff > flux_lim:
                    num_detection += 1
                    #print SN_array[n]
                    #print detected_array
                    detected_array = np.append(detected_array, SN_array[n])
            end3 = time.time()
            print("Completed one elif statement in %i seconds"%(end3-start3))
        
        #   change to suit the format of the lightcurve THIS PART IS INCOMPLETE. For now we will assume that the data is given in a format with days in one column (the first column) and the magnitudes in another column (the second column)
        else:
            num_detection = 0
            detected_array = np.recarray((0,), dtype=[('sn_type', '|S10'), ('ra', float), ('dec',float), ('redshift', float), ('magnification', float), ('epoch', int), ('apparent_mag_1', float), ('apparent_mag_2', float), ('mag_difference', float) , ('cluster', '|S10'), ('cadence', int)])
            start2 = time.time()
            for n in (np.arange(len(SN_array)-1) + 1):
                day = SN_array[n][5]
                mag_1 = new_SN_mag.SNmag(type, SN_array[n][3], day, fil=filter_name, bandwidth=bandwidth) - 0.5
                mag_11 = mumag(mag_1, SN_array[n][4])
                mag_2 = new_SN_mag.SNmag(type, SN_array[n][3], day + cadence, fil=filter_name, bandwidth=bandwidth) - 0.5
                mag_22 = mumag(mag_2, SN_array[n][4])
                SN_array[n][6] = mag_11
                SN_array[n][7] = mag_22
                flux_1 = astropysics.phot.mag_to_lum(mag_11)
                flux_2 = astropysics.phot.mag_to_lum(mag_22)
                flux_diff = np.abs(flux_1 - flux_2)
                mag_diff = astropysics.phot.lum_to_mag(flux_diff)
                SN_array[n][8] = mag_diff
                flux_lim = astropysics.phot.mag_to_lum(mag_lim)
                if flux_diff > flux_lim:
                    num_detection += 1
                    #print SN_array[n]
                    detected_array = np.append(detected_array, SN_array[n])
            end2 = time.time()
            print("completed one else statement in %i seconds."%(end2-start2))
    
    print SN_array
    #SN_array = np.delete(SN_array, 0, 0)
    #detected_array = np.delete(detected_array, 0, 0)
    print "number of detections", num_detection, np.count_nonzero(detected_array.size) != 0
    end = time.time()
    print("One run in %i seconds."%(end-start))

"""
if np.count_nonzero(detected_array) != 0:
    return num_detection, SN_array, detected_array
else:
    return num_detection, SN_array, np.zeros((1,9))
"""
if __name__=="__main__":
    c = {}
    execfile(configfile, c)
    
    min_mu = 18.9
    
    observations_file = c['observations']
    
    lensmodels = c['lensmodels']
    bandwidth = c['bandwidth']
    observations_file = c['observations']
    observations = asciitable.read(observations_file)
    
    clist = asciitable.read(lensmodels, names=['alias', 'clusterdir', 'deflxfile', 'deflyfile', 'zcluster', 'zmodel'])
    
    clusters = clist['alias']
    
    start = time.time()
    
    spacing = zmax - zmin + 1
    ztargets = np.linspace(zmin, zmax, spacing)
    
    #   get the filter name
    
    filter_name = glob.glob('../../filters/WFC*{}*.dat'.format(filter))[0]
    SN_array = np.recarray((0,), dtype=[('sn_type', '|S10'), ('ra', float), ('dec',float), ('redshift', float), ('magnification', float), ('epoch', int), ('apparent_mag_1', float), ('apparent_mag_2', float), ('mag_difference', float) , ('cluster', '|S10'), ('cadence', int), ('note', '|S16')])

    for x in range(len(ztargets)):
        #   this calculates the number of SN we expect at each magnification (or greater), but note that it won't be integer numbers at this point.
        number_arr = new_SN_number.number_SN(ztargets[x], mass, configfile, mu_min=mag_array[x,1])
        just_numbers = new_SN_number.recarr_to_nparr(number_arr) #things are working up to here at least
        
        nc = len(observations)
        
        #   before we were iterating over clusters...now I think it would be a better idea to iterate over the searches...or rather each column in the file. There doesn't seem to be a good way to do it otherwise.
        
        for y in np.linspace(1, nc, nc):
            
            clustername = observations.dtype.names[int(y)]
            
            #   turns the number at that magnification or greater into something that's just the number at that magnification.
            number = np.c_[just_numbers[0:len(just_numbers)-1,0], np.abs(np.diff(just_numbers[:,y]))]
            
            #   the fractional SN are treated as probabilities.
            probs = number[:,1] - np.fix(number[:,1])
            
            #   now we turn the numbers into integer numbers to put into the rest of the code
            number[:,1] = np.fix(number[:,1])
            
            #   The Monte Carlo part: draw from a random uniform distribution between 0 and 1. If the drawn number is less than probs, then you get a SN, if it is greater then you don't.
            for i in range(len(probs)):
                k = np.random.uniform(0, 1)
                
                if probs[i] > k:
                    number[i, 1] += 1
            
            #   we have integer numbers, so now we are going to make that many supernovae at that redshift with those magnifications, and we can automatically assign them their type, redshift, cluster number, and magnification now.
            for i in range(len(number)):
                SNs = np.recarray((number[i,1],), dtype=[('sn_type', '|S10'), ('ra', float), ('dec',float), ('redshift', float), ('magnification', float), ('epoch', int), ('apparent_mag_1', float), ('apparent_mag_2', float), ('mag_difference', float) , ('cluster', '|S10'), ('cadence', int), ('search', '|S16')])
                if SNs.size > 0:
                    SNs['sn_type'] = type #mass of the SN
                    SNs['redshift'] = ztargets[x] #redshift
                    SNs['magnification'] = number[i,0] #magnification
                    SNs['cluster'] = clustername #which cluster it is by name
                    
                    #SNs[:,5] = np.random.uniform(1, 1000) #days after explosion --> need to not assign the same day to every one of them! do this with the RA and DEC part
                    
                    #   now assigning the RA and DEC of each supernova using the regions files. Weight.fits is a map where pixels in the cluster region get 1 and outside get 0.
                    hdulist = pyfits.open("../../{}/weight.fits".format(clusters[y-1]))
                    ra_decdata = hdulist[0].data
                    #   print hdulist[0].header
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
                            SNs[n][1] = RA
                            SNs[n][2] = DEC
                            
                            #   here it is assigning each supernova a random date of observation based on the observation window. Earliest date should be the number of days before a cadence. The latest should be after the brightest supernova with longest plateau/brightest magnification/highest redshift can be seen.
                            SNs[n][5] = np.random.uniform(-1*(cadence), maxvisibility)
                            n += 1
                # join the new SN with the master list of SN already created.
                SN_array = np.append(SN_array, SNs)
        
        #   Detections
        #
        #   basically here everything that is created is detected if we put in no detection limits
        if detect_lim == False:
            num_detection = len(SN_array)
            detected_array = np.recarray((0,), dtype=[('sn_type', '|S10'), ('ra', float), ('dec',float), ('redshift', float), ('magnification', float), ('epoch', int), ('apparent_mag_1', float), ('apparent_mag_2', float), ('mag_difference', float) , ('cluster', '|S10'), ('cadence', int), ('search', '|S16')])
            print("ran if statement")
        
        #   if we do have detection limits, then the program will first try to see if there is an output lightcurve file (since this saves a lot of time)
        elif (detect_lim == True) & (os.path.isfile('../../output/lightcurve_{}_{}_{}.dat'.format(type, ztargets[x], filter))):
            start3 = time.time()
            num_detection = 0
            detected_array = np.recarray((0,), dtype=[('sn_type', '|S10'), ('ra', float), ('dec',float), ('redshift', float), ('magnification', float), ('epoch', int), ('apparent_mag_1', float), ('apparent_mag_2', float), ('mag_difference', float) , ('cluster', '|S10'), ('cadence', int), ('search', 'S16')])
            lightcurve = '../../output/lightcurve_{}_{}_{}.dat'.format(type, ztargets[x], filter)
            
            #   now use the light curves to calculate the magnitude of the supernova at the date of observation by seeing how bright it is for a light curve at that day, and then magnifying it with mu
            for n in (np.arange(len(SN_array)-1) + 1):
                day = SN_array[n][5]
                mu = SN_array[n][4]
                #   print( len(SN_array))
                redshift = SN_array[n][3]
                dat = ascii.read(lightcurve)
                dat = np.asarray(dat)
                lc_dat = dat.view(np.float).reshape(dat.shape + (-1,))
                
                #   finding the magnitudes at date of first observation and date of follow-up
                idx = find_nearest_idx(lc_dat[:,0], day)
                idx2 = find_nearest_idx(lc_dat[:,0], day + cadence)
                
                #   just picked a large number to indicate that the SN hadn't occurred yet, since a matching index of zero indicates that it's before it exploded. Don't bother magnifying these, so the original and magnified magnitudes are both set to 500.
                if idx == 0:
                    mag_1 = 500.
                    mag_11 = 500.
                else:
                    mag_1 = lc_dat[idx, 1] #unmagnified magnitude
                    #print "unmagnified mag_1, day, idx, and mu", mag_1, day, idx, mu
                    mag_11 = mumag(mag_1, mu)
                    print mag_1, mag_11, mu, n
                #time.sleep(1)
                #print "magnified mag_1", mag_1
                if idx2 ==0:
                    mag_2 = 500.
                    mag_22 = 500.
                else:
                    mag_2 = lc_dat[idx2, 1]
                    mag_22 = mumag(mag_2, mu)
                    #time.sleep(1)
                    print mag_2, mag_22, mu, n
                
                #   enter the final magnified magnitudes into the array of SNe
                SN_array[n][6] = mag_11
                SN_array[n][7] = mag_22
                
                #   calculate the difference between the two observations to see if it met the magnitude limit
                flux_1 = astropysics.phot.mag_to_lum(mag_11)
                flux_2 = astropysics.phot.mag_to_lum(mag_22)
                flux_diff = np.abs(flux_1 - flux_2)
                mag_diff = astropysics.phot.lum_to_mag(flux_diff)
                SN_array[n][8] = mag_diff
                flux_lim = astropysics.phot.mag_to_lum(mag_lim)
                
                #   if greater than threshold, count as detection
                if flux_diff > flux_lim:
                    num_detection += 1
                    #print SN_array[n]
                    #print detected_array
                    detected_array = np.append(detected_array, SN_array[n])
            end3 = time.time()
            print("Completed one elif statement in %i seconds"%(end3-start3))
        
        #   change to suit the format of the lightcurve THIS PART IS INCOMPLETE. For now we will assume that the data is given in a format with days in one column (the first column) and the magnitudes in another column (the second column)
        else:
            num_detection = 0
            detected_array = np.recarray((0,), dtype=[('sn_type', '|S10'), ('ra', float), ('dec',float), ('redshift', float), ('magnification', float), ('epoch', int), ('apparent_mag_1', float), ('apparent_mag_2', float), ('mag_difference', float) , ('cluster', '|S10'), ('cadence', int)])
            start2 = time.time()
            for n in (np.arange(len(SN_array)-1) + 1):
                day = SN_array[n][5]
                mag_1 = new_SN_mag.SNmag(type, SN_array[n][3], day, fil=filter_name, bandwidth=bandwidth) - 0.5
                mag_11 = mumag(mag_1, SN_array[n][4])
                mag_2 = new_SN_mag.SNmag(type, SN_array[n][3], day + cadence, fil=filter_name, bandwidth=bandwidth) - 0.5
                mag_22 = mumag(mag_2, SN_array[n][4])
                SN_array[n][6] = mag_11
                SN_array[n][7] = mag_22
                flux_1 = astropysics.phot.mag_to_lum(mag_11)
                flux_2 = astropysics.phot.mag_to_lum(mag_22)
                flux_diff = np.abs(flux_1 - flux_2)
                mag_diff = astropysics.phot.lum_to_mag(flux_diff)
                SN_array[n][8] = mag_diff
                flux_lim = astropysics.phot.mag_to_lum(mag_lim)
                if flux_diff > flux_lim:
                    num_detection += 1
                    #print SN_array[n]
                    detected_array = np.append(detected_array, SN_array[n])
            end2 = time.time()
            print("completed one else statement in %i seconds."%(end2-start2))

    print SN_array
    #SN_array = np.delete(SN_array, 0, 0)
    #detected_array = np.delete(detected_array, 0, 0)
    print "number of detections", num_detection, np.count_nonzero(detected_array.size) != 0
    end = time.time()
    print("One run in %i seconds."%(end-start))

#
#   Also write a function that is a bit more general than the one above (can take in any SN spectra and run it through SN Cosmo). The first one only uses the original spectra from Whalen et al.
#
