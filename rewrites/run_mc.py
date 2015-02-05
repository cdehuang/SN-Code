#
#   chuang
#
#   02.01.15
#
#
######

import numpy as np
import time

import conversion
import revised_mc

#
#
#   User-modifiable code to run the Monte Carlo simulation and print the data into an output file.
#
#
#
#   PARAMETERS
#
#   lensmodels: file containing information on the lens models necessary to run the simulation.
#   first column should be 'alias' (name of the cluster), the second column the directory in which the files are stored, the third column the name of the x deflection map, the fourth one column the name of the y deflection map, the fifth column is the redshift of the cluster, and the sixth column is the redshift of the model
#
#   observations: csv file containing the alias of all of the observed clusters, their effective cadence, reddest filter, and limiting magnitude
#   Doesn't matter how many columns you have or what order they're in as long as it's a csv in which the first row is headers and the cluster aliases are under the 'cluster' column, the effective cadences are under the 'effective_cadence' column and the reddest filters are under the 'reddest_filter' column. Limiting magnitudes should be in the 'limiting_mag_5sigma_ab' column...or change what you named it in the code
#
#   lightcurves_location: where the light curves data are located
#
#   lensed_volumes: where the lensed cluster volume data is located (filename and path)
#
#   volumes_fname: name of the file containing all of the volumes. new_SN_number.py will look for files with the format "lensed_volumes + volumes_fname + '_{}'.format(z)" where z is the redshift for which the volumes are being calculated. It also works best if the volumes file contains all of the volumes for all of the galaxy clusters that we want to calculate.
#
#   mc_output_loc: where to put the monte carlo outputs
#
#   mc_output_filename: what to call the output files of the monte carlo
#
#   star_formation_rate: the star formation rate as a function of redshift. to use the SFR in Whalen et al's paper, type 'high' for their high SFR and 'low' for their low SFR. Otherwise, input in an array of redshifts (first column) and associated star formation rates (given in terms of solar mass per year per Mpc^3 in the rest frame of the source--second column) and SN_number will just do a simple interpolation over those points
#
#   sn_types: which supernovae types you want to simulate
#
#   sn_masses: in solar masses
#
#   detection_limits: if the detection limits are set to True, then the limiting magnitude is used as the detection limit. If they are set to False, then you will just detect all of the SN that are created
#
#   zmin: the minimum redshift to simulate
#
#   zmax: the maximum redshift to simulate
#
#   filter: which filter you are using/location of the filter file (so, you can write 'F140W' or the full path + filename
#
#   maximum_visibility: the greatest amount of time for which the supernovae are visible, in days
#
#   bandwidth: the bandwidth of the filter you are observing in, in angstroms
#
#   number_of_runs: the number of times to run the monte carlo
#
#   efficiency: eta (the fraction of baryons converted into stars)
#
#   imf: the initial mass function to use (options are 'Kroupa', 'Salpeter' and 'KroupaSalpeter')

lensmodels = '../../data/vp_lensmodels.lis'
observations = '../../all_hst_surveys.csv'
lightcurves_location = '../../output/'
lensed_volumes = '../../output/'
volumes_fname = 'volumemag_whalen'
mc_output_loc = '../../output/'
mc_output_filename = 'mc_output'
star_formation_rate = 'high'
sn_types = ['z15g', 'z25g', 'z40g']
sn_masses = [15., 25., 40.]
detection_limits = True
zmin = 5
zmax = 12
filter = 'F140W'
maximum_visibility = 700
bandwidth = 3840
number_of_runs = 100
efficiency = 1
IMF = 'Kroupa'

#   separate into a few runs based on limiting magnitude

obs_array = np.recfromcsv(observations)
fnames = 'file'
new_files = conversion.make_observation_files(input_file=observations, output_fname=fnames)
for i in new_files:
    arr = np.recfromcsv(i)
    limiting_mag = arr['limiting_mag_5sigma_ab']
    maglim = limiting_mag[0]
    detections = []
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
            outdat = revised_mc.simulation(parameter_file=lensmodels, mag_lim=maglim,filter=filter,  detectlim=detection_limits, observations_file=observations, zmin=zmin, zmax=zmax, type=sn_types[s], mass=sn_masses[s], sfr=star_formation_rate, eta=efficiency, IMF=IMF, maxvisibility=maximum_visibility, volume_file_dir=lensed_volumes, volume_files_names=volumes_fname)
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





#





