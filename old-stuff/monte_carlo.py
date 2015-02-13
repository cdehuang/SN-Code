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

#    First thing the code does is to process the observations. It then divides everything by limiting magnitude? Well presumably, the surveys should be run one at a time. In that case, we would assume that everything has the same limiting magnitude...but at any rate, it at least means that if you want to go through a whole survey at once, you'll need to go through a bunch of different clusters at all once. I guess the user can add together more complicated surveys with many magnitudes.
#    What I want to do here is get all of the observations, and if there are multiple magnitudes, break it up into multiple surveys.
#    Check that all of the files that we need are there

#    Loads the observations as a record array, returns an error if there isn't an observations file.
observations = c['observations']
if not (os.path.isfile(observations)):
    print("missing observations file")
    sys.exit()
olist = ascii.read(observations)
olist = np.asarray(olist)


#    Checks that the magnified volumes data files are available. Spits out code explaining which ones are missing if they aren't.
#    May change this so that later it will also make these files.
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


#    Apply the detection limits
sn_types = c['sn_types']
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

        
#   Now calculate the number of SN expected per row of the observations. However note that this can only do it for one redshift at a time. Pass the output array into the revised SN code?
#   So it should calculate an array for each redshift where you have the magnifications (first column) then each cluster/survey in the other ones. It'll be a record array, so later we can grab the values if we know the cluster names and we don't need to do repeats.
for z in target_z:
    for s in sn_type:
        minmu =
        new_SN_number.sn_num_cluster(z, mass, configfile, mu_min=mmu)

