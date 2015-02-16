import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib
import numpy as np
import matplotlib.pylab as plt
from scipy.interpolate import UnivariateSpline
from scipy.integrate import quad
from StringIO import StringIO
import pyfits
import os
import math
import new_volfun
import cosmolopy.distance as cd
import astropy.io.ascii as ascii
import scipy

def recarr_to_nparr(record_array):
    record_array = np.asarray(record_array)
    array = record_array.view(np.float).reshape(record_array.shape + (-1,))
    return array

#   since I find myself doing this all the time, this converts a .dat file into just a regular numpy array instead of a record array

def dat_to_nparr(filename):
    vlist = ascii.read(filename)
    vlist = np.asarray(vlist)
    array = vlist.view(np.float).reshape(vlist.shape + (-1,))
    return array

#this is the version of SN_number.py that is located in the "Debugging Things" directory

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx

def highSFR(x):
    #corresponds to the optimistic case in the Whalen paper
    if (x < 12):
        val = 0.2
    else:
        val =  np.power(10, -0.25*x + 2.3)
    return val

def lowSFR(x):
    #corresponds to the lower limit in the Whalen paper
    val = np.power(10, -0.25*x - 0.4489)
    return val

def Salpeter(x):
    return np.power(x, -2.35)

def Salpeter2(x):
    return np.power(x, -1.35)

def Kroupa(x):
    if x >= 1.:
        val = np.power(x, -7./3)
    elif (x < 1.) & (x >= 0.1):
        val = np.power(x, -4./3)
    else:
        val = np.power(x, -1/3.)
    return val

#Basically the Kroupa, but multiplied by the mass of the star as well so instead of number of stars you get the amoutn of stellar mass
def Kroupa2(x):
    if x >= 1:
        val = np.power(x, -4./3)
    elif (x < 1) & (x >= 0.1):
        val = np.power(x, -1./3)
    elif (x < 0.1) & (x >= 0.07):
        val = np.power(x, 2./3)
    else:
        val = 0
    return val

def KroupaSalpeter(x):
    if x >= 0.5:
        val = np.power(x, -2.35)
    elif (x < 0.5) & (x >= 0.08):
        val = np.power(x, -1.3)
    elif (x < 0.08) & (x >= 0.07):
        val = np.power(x, -0.3)
    else:
        val = 0
    return val

def KroupaSalpeter2(x):
    if x >= 0.5:
        val = np.power(x, -1.35)
    elif (x < 0.5) & (x >= 0.08):
        val = np.power(x, -0.3)
    elif (x < 0.08) & (x >= 0.07):
        val = np.power(x, .7)
    else:
        val = 0
    return val

def KroupaTest(x): #for finding the fraction of the mass that goes into SN
    val = np.power(x, -4./3)
    return val

#gives you a luminosity distance in Mpc
def lumdist(z, H_0=73, omega_M=0.308, omega_L=0.692, omega_K=0):
    pc_cm = 3.08567758E18
    c = 299792.458 #km/s
    D_H = c/H_0
    I = lambda z1: 1/np.sqrt(omega_M*(1+z1)**3 + omega_K*(1+z1)**2 + omega_L)
    Int, err = quad(I, 0, z)
    D_m = D_H*Int
    D_a = D_m/(1+z)
    D_l = ((1 + z)**2)*D_a
    lum_dist = D_l
    return lum_dist

def sphereshell(r):
    shell = 4*np.pi*np.power(r, 2)
    return shell

#as duncan noted, the distance between z=1 and z=2 should be a lot smaller than the distance between z=5 and z=6--> calc the comoving radial distance rather than the luminosity distance

def comovingradial(z, H_0=73, omega_M=0.308, omega_L=0.692, omega_K=0):
    c = 299792.458 #km/s
    D_H = c/H_0
    #I = lambda z1: 1/np.sqrt[
    comoving_dist = 10
    return comoving_dist


#   returns the number of SN expected at the given z for each mu. First column is the mu. Each of the later columns tells you how many SN are found in each cluster...seems like that's probably something that we would want to pass in in the future, or make a document that everything references since it's something I'll have to change later. Or, alternatively, I could make another function that does read that in...
#   mass is mass of star.
#   There needs to be something where it is only taking in the relevant clusters...maybe the full code does that, but it's not really a good idea. Perhaps it could read in the same input file as the other stuff.
#   maxvisibility is the longest ay of the SN could be seen at any redshift and magnification
#   need to remove the repetitive columns

def number_SN(z, mass, configfile, mu_min=1):
    """
        Calculates the number of SN expected per cluster per redshift per magnification.
    """
    #   Load in all the parameters
    c = {}
    execfile(configfile, c)
    
    if mu_min < 50:
        mu_max = 60
    else:
        mu_max = mu_min + 10
    
    efficiency = c['efficiency']
    
    #   cluster_file is the file that contains all the clusters, their reddest filters, and effective cadences
    cluster_array = ascii.read(c['observations'])
    
    clusters = cluster_array['cluster']
    cadences = cluster_array['effective_cadence']
    filter = cluster_array['reddest_filter']
    maxvisibility = c['maximum_visibility']
    search_number = cluster_array['search_number']
    sfr = c['star_formation_rate']
    
    IMF = eval(c['IMF'])
    IMF2 = eval(c['IMF2'])
    try:
        print "IMF found."
        all = quad(IMF2, 0.1, 100)
        allsn = all[0]*efficiency
        if mass <= 5:
            cc = quad(Kroupa, mass-0.5, mass+0.5)
        elif (5 < mass) & (mass <= 10):
            cc = quad(Kroupa, mass-1, mass+1)
        elif (10 < mass) & (mass <= 40):
            cc = quad(Kroupa, mass-5, mass+5)
        elif (mass > 40):
            cc = quad(Kroupa, mass-10, mass+10)
    except AttributeError:
        print("{} IMF unknown. Using Kroupa IMF instead...".format(IMF))
        all = quad(Kroupa2, 0.1, 100) #using Kroupa 2 to calculate the total stellar mass
        allsn = all[0]*efficiency
        if mass <= 5:
            cc = quad(Kroupa, mass-0.5, mass+0.5)
        elif (5 < mass) & (mass <= 10):
            cc = quad(Kroupa, mass-1, mass+1)
        elif (10 < mass) & (mass <= 40):
            cc = quad(Kroupa, mass-5, mass+5)
        elif (mass > 40):
            cc = quad(Kroupa, mass-10, mass+10)
    
    fraction = cc[0]/all[0] #the fraction of all the stellar mass created that goes into the stars we're care about.
    ccsn = fraction*efficiency
    
    #the actual amount of stellar mass per Mpc^3, also the interpolation part assumes that we get a bunch of points that are arranged in a 2-column array where the first column is redshift and the second one is the star formation rate I guess
    if sfr == 'high':
        I = quad(highSFR, z, z+1)
    elif sfr == 'low':
        I = quad(lowSFR, z, z+1)
    else:
        s = scipy.interpolate.interp1d(sfr[:,0], sfr[:,1])
        I = quad(s, z, z+1)
    
    numccsn = ccsn*I[0] #number of cc sn per Mpc^3 using their numbers it comes out to like 0.00024 SN per Mpc^3
    
    #for testing purposes, take out the dependence on Leonidas' code and see if we are still creating SNs. Just integrate the area
    if os.path.isfile('../../output/volumemag_whalen_{}.dat'.format(z)):
        #if os.stat('../output/volumemag_whalen_{}.dat'.format(z)).st_size != 0:
        magvolumes = '../../output/volumemag_whalen_{}.dat'.format(z)
        print "using the Whalen et al cosmology"
        vlist = ascii.read(magvolumes)
        vlist = np.asarray(vlist)
        volarray = vlist.view(np.float).reshape(vlist.shape + (-1,))
        volarray = volarray[mu_min-1: mu_max]
    elif os.path.isfile(volume_file_dir + volume_files_name + '_{}'.format(z)):
        magvolumes = volume_file_dir + volume_files_name + '_{}'.format(z)
        print "using the standard LCDM cosmology"
        vlist = ascii.read(magvolumes)
        vlist = np.asarray(vlist)
        volarray = vlist.view(np.float).reshape(vlist.shape + (-1,))
        volarray = volarray[mu_min-1: mu_max]
    elif os.path.isfile(volume_file_dir + volume_files_name + '_{}'.format(z)):
        magvolumes = '../../output/volumemag_{}.dat'.format(z)
        print "using user-input volumes and cosmology"
        vlist = ascii.read(magvolumes)
        vlist = np.asarray(vlist)
        volarray = vlist.view(np.float).reshape(vlist.shape + (-1,))
        volarray = volarray[mu_min-1: mu_max]
    else:
        print "Volume files not found. Calculating volumes..."
        volarray = new_volfun.volwmag(z, mumin=mu_min, mumax=mu_max, parameter_file=parameter_file)
    
    #   Runs through all of the clusters, gets their data, should return an array in which each column is the volumes of the clusters in order. first column is the magnifications
    
    all_cluster_vols = vlist['magnification']
    
    missing_clusters = []
    arr_columns = ['magnification']
    
    #   note that the macs clusters had Js in the name, but that had gotten dropped in earlier iterations of code, so now in order to match them properly i'm putting it back in.
    i = 0
    for cs in clusters:
        clust = cs.lower()
        if clust[0:5] == 'macsj':
            clust = ''.join([clust[0:4], clust[-4:]])
        try:
            clust_vols = vlist[clust]
            arr_columns.append(clust+'_s{}'.format(search_number[i]))
        except ValueError:
            missing_clusters.append(clust)
            continue
        all_cluster_vols = np.c_[all_cluster_vols, clust_vols]
        i += 1
    
    if len(missing_clusters) > 0:
        missing_clusters = set(missing_clusters)
        missing_cluster_indices = np.zeros(0)
        for i in missing_clusters:
            clusters = np.asarray([cs.lower() for cs in clusters])
            indices = np.where(clusters == i)
            missing_cluster_indices = np.r_[missing_cluster_indices, indices[0]]
    #   trim to the appropriate size
    all_cluster_vols = all_cluster_vols[mu_min-1:mu_max]
    arrshape = all_cluster_vols.shape
    
    #   note that this could be an array, depending on the input file
    duration = cadences + maxvisibility
    
    duration = duration/365.
    #print "duration", duration, duration.shape
    
    #time = (duration)/(1.+z) #using the survey duration in years in our frame
    time = (duration)/(1. + z)
    #print "time", time
    time = np.delete(time, missing_cluster_indices)
    
    #   multiple cluster version
    numarray = all_cluster_vols*1.
    clusts = len(all_cluster_vols[0]) -1
    for i in range(clusts):
        #numarray[:,i] = numarray[:,i]*numccsn
        numarray[:,i+1] = numarray[:,i+1]*numccsn*time[i] #for one of the cases
    #number_arr = np.c_[ number_arr, [vol*SFR[z]*eta[z]]]

    print "missing clusters:", missing_clusters
    #print "shape of the nparray", numarray.shape
    num_arr = np.core.records.fromarrays(numarray.transpose(), names=arr_columns)

    return num_arr

#   debugging
if __name__ =="__main__":
    z = 5.
    mass = 40.
    mu_min=50
    c = {}
    execfile(configfile, c)

    if mu_min < 50:
        mu_max = 60
    else:
        mu_max = mu_min + 10

    efficiency = c['efficiency']

    #   cluster_file is the file that contains all the clusters, their reddest filters, and effective cadences
    cluster_array = ascii.read(c['observations'])

    clusters = cluster_array['cluster']
    cadences = cluster_array['effective_cadence']
    filter = cluster_array['reddest_filter']
    maxvisibility = c['maximum_visibility']
    search_number = cluster_array['search_number']
    sfr = c['star_formation_rate']

    IMF = eval(c['IMF'])
    IMF2 = eval(c['IMF2'])
    try:
        all = quad(IMF2, 0.1, 100)
        allsn = all[0]*efficiency
        if mass <= 5:
            cc = quad(Kroupa, mass-0.5, mass+0.5)
        elif (5 < mass) & (mass <= 10):
            cc = quad(Kroupa, mass-1, mass+1)
        elif (10 < mass) & (mass <= 40):
            cc = quad(Kroupa, mass-5, mass+5)
        elif (mass > 40):
            cc = quad(Kroupa, mass-10, mass+10)
    except AttributeError:
        print("{} IMF unknown. Using Kroupa IMF instead...".format(IMF))
        all = quad(Kroupa2, 0.1, 100) #using Kroupa 2 to calculate the total stellar mass
        allsn = all[0]*efficiency
        if mass <= 5:
            cc = quad(Kroupa, mass-0.5, mass+0.5)
        elif (5 < mass) & (mass <= 10):
            cc = quad(Kroupa, mass-1, mass+1)
        elif (10 < mass) & (mass <= 40):
            cc = quad(Kroupa, mass-5, mass+5)
        elif (mass > 40):
            cc = quad(Kroupa, mass-10, mass+10)

    fraction = cc[0]/all[0] #the fraction of all the stellar mass created that goes into the stars we're care about.
    ccsn = fraction*efficiency

    #the actual amount of stellar mass per Mpc^3, also the interpolation part assumes that we get a bunch of points that are arranged in a 2-column array where the first column is redshift and the second one is the star formation rate I guess
    if sfr == 'high':
        I = quad(highSFR, z, z+1)
    elif sfr == 'low':
        I = quad(lowSFR, z, z+1)
    else:
        s = scipy.interpolate.interp1d(sfr[:,0], sfr[:,1])
        I = quad(s, z, z+1)

    numccsn = ccsn*I[0] #number of cc sn per Mpc^3 using their numbers it comes out to like 0.00024 SN per Mpc^3

    #for testing purposes, take out the dependence on Leonidas' code and see if we are still creating SNs. Just integrate the area
    if os.path.isfile('../../output/volumemag_whalen_{}.dat'.format(z)):
        #if os.stat('../output/volumemag_whalen_{}.dat'.format(z)).st_size != 0:
        magvolumes = '../../output/volumemag_whalen_{}.dat'.format(z)
        print "using the Whalen et al cosmology"
        vlist = ascii.read(magvolumes)
        vlist = np.asarray(vlist)
        volarray = vlist.view(np.float).reshape(vlist.shape + (-1,))
        volarray = volarray[mu_min-1: mu_max]
    #print volarray
    #print mu_min, mu_max
    #clusters = vlist.dtype.names
    elif os.path.isfile(volume_file_dir + volume_files_name + '_{}'.format(z)):
        magvolumes = volume_file_dir + volume_files_name + '_{}'.format(z)
        print "using the standard LCDM cosmology"
        vlist = ascii.read(magvolumes)
        vlist = np.asarray(vlist)
        volarray = vlist.view(np.float).reshape(vlist.shape + (-1,))
        volarray = volarray[mu_min-1: mu_max]
    elif os.path.isfile(volume_file_dir + volume_files_name + '_{}'.format(z)):
        magvolumes = '../../output/volumemag_{}.dat'.format(z)
        print "using user-input volumes and cosmology"
        vlist = ascii.read(magvolumes)
        vlist = np.asarray(vlist)
        volarray = vlist.view(np.float).reshape(vlist.shape + (-1,))
        volarray = volarray[mu_min-1: mu_max]
    #print mu_min, mu_max
    #clusters = vlist.dtype.names
    else:
        print "Volume files not found. Calculating volumes..."
        volarray = new_volfun.volwmag(z, mumin=mu_min, mumax=mu_max, parameter_file=parameter_file)

    #   Runs through all of the clusters, gets their data, should return an array in which each column is the volumes of the clusters in order. first column is the magnifications

    all_cluster_vols = vlist['magnification']

    missing_clusters = []
    arr_columns = ['magnification']

    #   note that the macs clusters had Js in the name, but that had gotten dropped in earlier iterations of code, so now in order to match them properly i'm putting it back in.
    i = 0
    for cs in clusters:
        clust = cs.lower()
        if clust[0:5] == 'macsj':
            clust = ''.join([clust[0:4], clust[-4:]])
        try:
            clust_vols = vlist[clust]
            arr_columns.append(clust+'_s{}'.format(search_number[i]))
        except ValueError:
            missing_clusters.append(clust)
            continue
        all_cluster_vols = np.c_[all_cluster_vols, clust_vols]
        i += 1

    if len(missing_clusters) > 0:
        missing_clusters = set(missing_clusters)
        missing_cluster_indices = np.zeros(0)
        for i in missing_clusters:
            clusters = np.asarray([cs.lower() for cs in clusters])
            indices = np.where(clusters == i)
            missing_cluster_indices = np.r_[missing_cluster_indices, indices[0]]
    #   trim to the appropriate size
    all_cluster_vols = all_cluster_vols[mu_min-1:mu_max]
    arrshape = all_cluster_vols.shape

    #   note that this could be an array, depending on the input file
    duration = cadences + maxvisibility

    duration = duration/365.
    #print "duration", duration, duration.shape

    #time = (duration)/(1.+z) #using the survey duration in years in our frame
    time = (duration)/(1. + z)
    #print "time", time
    time = np.delete(time, missing_cluster_indices)

    #   multiple cluster version
    numarray = all_cluster_vols*1.
    clusts = len(all_cluster_vols[0]) -1
    for i in range(clusts):
        #numarray[:,i] = numarray[:,i]*numccsn
        numarray[:,i+1] = numarray[:,i+1]*numccsn*time[i] #for one of the cases
    #number_arr = np.c_[ number_arr, [vol*SFR[z]*eta[z]]]

    #   now make a record array with all of the information in there
    print "number of columns", len(arr_columns)
    print "missing clusters:", missing_clusters
    print "shape of the nparray", numarray.shape
    num_arr = np.core.records.fromarrays(numarray.transpose(), names=arr_columns)
    print "num_arr: shape", num_arr.shape


#   The code needs a volumes file. Or it needs to calculate the volumes in some way. It also needs to know the name of the clusters. To calculate the volumes, in general we just need to give the volfun code a directory with the various files.
#   Calculates the number of SN for a given volume of space based on the input star formation rates and stuff
#def sn_num_cluster(z, mass, configfile, mu_min=1):



