#
#    chuang
#    02.12.15
#
#    The final iteration of the Monte Carlo code for running everything. Only necessary inputs are the parameter file and the 
#
#

import numpy as np
import math
import glob
import os
import astropy.io.ascii as ascii
import sys
import time

import new_volfun
import new_SN_number
import new_SN_mag
import new_maglim
import revised_mc as rmc
#import revised_mc

#    necessary functions



#    start of the main body of code
if len(sys.argv)>1:
    configfile = sys.argv[1]
else:
    configfile = 'sn_params_test.cfg'
    #    fix this part of the code later
if not (os.path.isfile(configfile)):
    print("missing startup file")
    sys.exit()
    
c = {}
execfile(configfile, c)

#   First thing the code does is to process the observations. It then divides everything by limiting magnitude? Well presumably, the surveys should be run one at a time. In that case, we would assume that everything has the same limiting magnitude...but at any rate, it at least means that if you want to go through a whole survey at once, you'll need to go through a bunch of different clusters at all once. I guess the user can add together more complicated surveys with many magnitudes.
#   What I want to do here is get all of the observations, and if there are multiple magnitudes, break it up into multiple surveys.
#   Check that all of the files that we need are there

#   Loads the observations as a record array, returns an error if there isn't an observations file.
observations = c['observations']
if not (os.path.isfile(observations)):
    print("missing observations file")
    sys.exit()
olist = ascii.read(observations)
olist = np.asarray(olist)


#   Checks that the magnified volumes data files are available. Spits out code explaining which ones are missing if they aren't.
#   May change this so that later it will also make these files.
zmin = c['zmin']
zmax = c['zmax']
lensdir = c['lensed_volumes_directory']
vfname = c['volumes_fname']
target_z = np.arange(zmin, zmax + 1, dtype='float')
for z in target_z:
    volumes_filename = lensdir+vfname+'_{}.dat'.format(z)
    if os.path.isfile(volumes_filename):
        continue
    else:
        print("Missing {}.".format(volumes_filename))


#   Apply the detection limits
sn_types = c['sn_types']
sn_masses = c['sn_masses']
zipped_sns = zip(sn_types, sn_masses)
filter = c['filter']
fildir = c['filter_directory']
filter_name = glob.glob('{}WFC*{}*.dat'.format(fildir, filter))[0]
lim_mag = c['limiting_magnitude']

start = time.time()

datatypes = [('redshift', float)]
for i in sn_types:
    datatypes.append((i, float))
mag_array = np.recarray((len(target_z),), dtype= datatypes)
mag_array['redshift'] = target_z
i = 0
for s in sn_types:
    n = 0
    for z in target_z:
        if os.path.isfile('../../output/lightcurve_{}_{}_{}.dat'.format(s, z, filter)):
            lightcurve = '../../output/lightcurve_{}_{}_{}.dat'.format(s, z, filter)
            dat = ascii.read(lightcurve)
            dat = np.asarray(dat)
            lc_dat = dat.view(np.float).reshape(dat.shape + (-1,))
            mags = lc_dat[:,1]
            days = lc_dat[:,0]
            trimmed_mags, idx = new_maglim.remove_shock_breakout(mags)
            mu = new_maglim.magfactor(trimmed_mags, lim_mag)
        else:
             print("Missing lightcurve_{}_{}_{}.dat. Calculating light curve.".format(s, z, filter))
             l_c = new_maglim.whalen_lightcurve(s, z, fil=filter_name, bandwidth=bandwidth)
             mu = new_maglim.magfactor(l_c, lim_mag)
        mag_array[s][n] = mu
        n += 1

#   Pass the appropriate parts of the magnification array into the monte-carlo code?

        
#   Now calculate the number of SN expected per row of the observations. However note that this can only do it for one redshift at a time. Pass the output array into the revised SN code?
#   So it should calculate an array for each redshift where you have the magnifications (first column) then each cluster/survey in the other ones. It'll be a record array, so later we can grab the values if we know the cluster names and we don't need to do repeats.
#   need to find a way to have this change it's size, since the number of parameters that get output will depend on the number of
#   Iterate over all redshifts and all sn_types
sys.exit()
detections = []
supernovae_info = np.zeros((1, 10))
detect_info = np.zeros((1, 10))
all_sn_models = []
detected_sn_models = []
i = 0
while i < runs:
    totaldetect = 0
    totalinfo = np.zeros((1, 10))
    detectinfo = np.zeros((1, 10))
    for s in range(len(sn_types)):
        m_arr = np.c_[mag_array['redshift'], mag_array[s]]
        outdat = rmc.simulation(configfile, m_arr, sn_type)
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
rec_array['cluster'] = supernovae_info[:,8]


#for z in target_z:
#    for zz in zipped_sns:
#        minmu = mag_array[zz[0]][z-min(target_z)]
#        number_array = new_SN_number.number_SN(z, zz[1], configfile, mu_min=minmu)

