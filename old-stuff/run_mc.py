#
#   chuang
#
#   02.01.15
#
#
#   User-modifiable code to run the Monte Carlo simulation and print the data into an output file.
#
#
######

#from ConfigParser import SafeConfigParser
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

import conversion
import new_volfun
import new_SN_number
import new_SN_mag
import new_maglim

"""
Code to run the simulation. Need to fill out the parameter file in order for it to work correctly though.
"""

def mumag(mag, mu):
    """
        Takes in a magnitude and a magnification and returns a new magnitude after the flux has been magnified.
        Inputs: magnitude, magnification
    """
    jsky = np.power(10, 23 -(mag+ 48.60)/2.5)
    mu_jsky = mu*jsky
    mu_mag = 2.5*(23 - np.log10(mu_jsky)) - 48.60
    return mu_mag

def find_nearest(array, value):
    """
        Finds the value of the array that is nearest to the input value
        Inputs: the array, the value
    """
    idx = (np.abs(array-value)).argmin()
    return array[idx]

def find_nearest_idx(array, value):
    """
        Finds the index of the value in the array that is closest to the input value
        Inputs: the array, the value
    """
    idx = (np.abs(array-value)).argmin()
    return idx

#parser = SafeConfigParser()

if len(sys.argv)>1:
    configfile = sys.argv[1]
#    parser.read(configfile)
else:
    configfile = 'sn_params.cfg'
#    parser.read('sn_params.cfg')
if not (os.path.isfile(configfile)):
    print 'missing startup file'
    sys.exit()

c = {}
execfile(configfile, c)

#   separate into a few runs based on limiting magnitude
#   probably the easiest thing to do would be to run each line of the file--which basically means just one cluster at a time....

obs_array = np.recfromcsv(c['observations'])
lim_mag = c['limiting_magnitude']
new_files = conversion.make_observation_files(input_file=c['observations'], output_fname=c['fnames'])
for i in new_files:
    arr = np.recfromcsv(i)
    limiting_mag = arr[lim_mag]
    maglim = limiting_mag[0]
    detections = []
    #   need to find a way to have this change it's size, since the number of parameters that get output will depend on the number of
    supernovae_info = np.zeros((1, 10))
    detect_info = np.zeros((1, 10))
    all_sn_models = []
    detected_sn_models = []
    maglim =
    i = 0
    while i < runs:
        totaldetect = 0
        totalinfo = np.zeros((1, 10))
        detectinfo = np.zeros((1, 10))
        for s in range(len(sn_types)):
            outdat = rmc.simulation(parameter_file=lensmodels, mag_lim=maglim,filter=filter,  detectlim=detection_limits, observations_file=observations, zmin=zmin, zmax=zmax, type=sn_types[s], mass=sn_masses[s], sfr=star_formation_rate, eta=efficiency, IMF=IMF, maxvisibility=maximum_visibility, volume_file_dir=lensed_volumes, volume_files_names=volumes_fname)
            totaldetect = totaldetect + outdat[0]
            totalinfo = np.r_[totalinfo, outdat[1]]
            detect_sn = [sn_types[s]]*len(outdat[0])
            detected_sn_models = detected_sn_models.extend(detect_sn)
            all_sn = [sn_types[s]]*len(outdat[1])
            all_sn_models = all_sn_models.extend(all_sn)
            if np.count_nonzero(outdat[2]) != 0:
                detectinfo = np.r_[detectinfo, outdat[2]]
        totalinfo = np.delete(totalinfo, 0, 0)
        detectinfo = np.delete(detectinfo, 0, 0)
        detections.append(totaldetect)
        supernovae_info = np.r_[supernovae_info, totalinfo]
        detect_info  = np.r_[detect_info, detect_info]
        print("Finished run %i. %i more runs left"%((i+1), (runs-i-1)))
        i += 1
    supernovae_info = np.delete(supernovae_info, 0, 0)
    detect_info = np.delete(detect_info, 0, 0)
    rec_array = np.recarray((len(supernovae_info),), dtype=[('SN type', str, 4), ('RA', float), ('DEC', float), ('z', float), ('mu', float), ('epoch days', float), ('ABmag_1', float), ('ABmag_2', float), ('cluster', str, 16)])
    rec_array['SN type'] = all_sn_models
    rec_array['RA'] = supernovae_info[:,1]
    rec_array['DEC'] = supernovae_info[:,2]
    rec_array['z'] = supernovae_info[:,3]
    rec_array['mu'] = supernovae_info[:,4]
    rec_array['epoch days'] = supernovae_info[:,5]
    rec_array['ABmag_1'] = supernovae_info[:,6]
    rec_array['ABmag_2'] = supernovae_info[:,7]
    rec_array['cluster'] =

def run_mc(sn_params.cfg):
    print "working on it!"



#





