#
#    chuang
#    02.12.15
#
#    The final iteration of the Monte Carlo code for running everything. Only necessary inputs are the parameter file
#
#

import numpy as np
import math
import glob
import os
import astropy.io.ascii as ascii
import sys
import time
import matplotlib.pyplot as plt
import numpy.lib.recfunctions as rfn

import new_volfun
import new_SN_number
import new_SN_mag
import new_maglim
import revised_mc as rmc
#import revised_mc

#    necessary functions

"""
   Runs the other parts of the code necessary for the monte carlo. 
   This is the part that runs the simulation a bunch of times for the various supernovae and redshifts.
"""

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

#   Get all of the observations, and if there are multiple magnitudes, break it up into multiple surveys.
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

#   start of the actual simulation

start = time.time()

datatypes = [('redshift', float)]
for i in sn_types:
    datatypes.append((i, float))
mag_array = np.recarray((len(target_z),), dtype= datatypes)
mag_array['redshift'] = target_z
lightcurve_dir = c['lightcurves_location']
i = 0
for s in sn_types:
    n = 0
    for z in target_z:
        lightcurve_fname = lightcurve_dir + 'lightcurve_{}_{}_{}.dat'.format(s, z, filter)
        if os.path.isfile(lightcurve_fname):
            lightcurve = lightcurve_fname
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

#   So it should calculate an array for each redshift where you have the magnifications (first column) then each cluster/survey in the other ones. It'll be a record array, so later we can grab the values if we know the cluster names and we don't need to do repeats.

runs = c['number_of_runs']

#rec_array_datatypes =[('sn_type', '|S10'), ('ra', float), ('dec',float), ('redshift', float), ('magnification', float), ('epoch', int), ('apparent_mag_1', float), ('apparent_mag_2', float), ('mag_difference', float) , ('cluster', '|S10'), ('cadence', int), ('id', '|S10')]
output_rec_array_dtypes = [('sn_type', '|S10'), ('ra', float), ('dec',float), ('redshift', float), ('magnification', float), ('epoch', int), ('apparent_mag_1', float), ('apparent_mag_2', float), ('mag_difference', float) , ('cluster', '|S10'), ('cadence', int), ('id', '|S10'), ('MC_run', int)]
detections = []
supernovae_info = np.zeros((0,), dtype=output_rec_array_dtypes)
detect_info = np.zeros((0,), dtype=output_rec_array_dtypes)
all_sn_models = []
detected_sn_models = []
i = 0
while i < runs:
    totaldetect = 0
    totalinfo = np.recarray((0,), dtype=output_rec_array_dtypes)
    detectinfo = np.recarray((0,), dtype=output_rec_array_dtypes)
    for s in zipped_sns:
        m_arr = np.c_[mag_array['redshift'], mag_array[s[0]]]
        outdat = rmc.simulate(configfile, m_arr, s[0], s[1], recarrdatypes=output_rec_array_dtypes)
        
        #   Inelegant code for adding in the run numbers.
        run1 = [i+1]*len(outdat[1])
        run2 = [i+1]*len(outdat[2])
        alloutput = outdat[1]
        detectedoutput = outdat[2]
        alloutput['MC_run'] = run1
        detectedoutput['MC_run'] = run2
        
        #   creating arrays of the detected supernovae and all of the supernovae
        totaldetect = totaldetect + outdat[0]
        totalinfo = rmc.concat((totalinfo, alloutput))
        detect_sn = [s[0]]*outdat[0]
        detected_sn_models = detected_sn_models + detect_sn
        all_sn = [s[0]]*len(alloutput)
        all_sn_models = all_sn_models + all_sn
        if np.count_nonzero(detectedoutput) != 0:
            detectinfo = rmc.concat([detectinfo, outdat[2]])
    detections.append(totaldetect)
    supernovae_info = rmc.concat((supernovae_info, totalinfo))
    detect_info  = rmc.concat((detect_info, detectinfo))
    print("Finished run %i. %i more runs left"%((i+1), (runs-i-1)))
    i += 1

output_dir = c['output_dir']
output_detected_fname = c['output_detected_fname']
output_all_sn_fname = c['output_all_sn_fname']

output_detected_filename = output_dir + output_detected_fname + '_{}_2.dat'.format(runs)
output_allsn_filename = output_dir + output_all_sn_fname + '_{}_2.dat'.format(runs)
if np.count_nonzero(detect_info) != 0:
    ascii.write(detect_info, output_detected_filename)
ascii.write(supernovae_info, output_allsn_filename)
print "Wrote output data for detected supernovae to", output_detected_filename
print "and output data for all the supernovae to", output_allsn_filename


end = time.time()
print "The Monte Carlo took %i seconds for %i simulations."%(end-start, runs)
plt.clf()
plt.hist(detections)
plt.show()

#for z in target_z:
#    for zz in zipped_sns:
#        minmu = mag_array[zz[0]][z-min(target_z)]
#        number_array = new_SN_number.number_SN(z, zz[1], configfile, mu_min=minmu)

