import os
from os.path import join, getsize, expanduser
import re
import fileinput
import numpy as np
import matplotlib.pylab as plt
import scipy
from scipy.interpolate import UnivariateSpline
from scipy import integrate
from StringIO import StringIO
import pyfits

#   This is supposed to get the magnitude of a supernova at a given day in its lightcurve
def SNmag(sn_type, redshift, t, fil='../filters/WFC3_IR_F160W.dat', bandwidth=2683):
    #indent everything to go under this function
    filenames = []
    for root, dirs, files in os.walk('../../Lanl.ccsn.spectra/{}'.format(sn_type)):
        for name in files:
            filenames.append(os.path.join(root, name))
    fnames = [i for i in filenames if re.match('../../Lanl.ccsn.spectra/{}/................/spectra.out.6.1'.format(sn_type),i) is not None]
    fnames = [i for i in fnames if os.stat(i).st_size !=0]
    #print len(fnames)
    filter = fil
    filterdat = np.loadtxt(filter)
    filwv = filterdat[:,0]
    filamp = filterdat[:,1]
    filt = scipy.interpolate.UnivariateSpline(filwv, filamp)
    flen = len(filterdat)
    #Assume a cosmology to calculate a distance. Using Ned's calculator for now
    z = redshift

    #l_d=6.87604E4 for z = 7
    #calculate luminosity distance from the redshift using the values listed in the paper
    H_0 = 69.6 #not listed in the paper, but perhaps they used a different one?
    omega_M = 0.308
    omega_L = 0.692
    omega_K = 0
    pc_cm = 3.08567758E18
    c = 299792.458 #km/s
    D_H = c/H_0
    I = lambda z1: 1/np.sqrt(omega_M*(1+z1)**3 + omega_K*(1+z1)**2 + omega_L)
    Int, err = integrate.quad(I, 0, z)
    D_m = D_H*Int
    D_a = D_m/(1+z)
    D_l = ((1 + z)**2)*D_a
    lum_dist = D_l*1E6*pc_cm
    beta = 3E-13

    #fixes errors in the spectrum_dumps.list without making a new s_d.list file
    dirs = os.listdir('../../LANL.CCSN.SPECTRA/{}/'.format(sn_type))
    dirs = [i for i in dirs if re.match('....-.-dmp......',i) is not None]
    timedat = np.genfromtxt("../../LANL.CCSN.SPECTRA/{}/spectrum_dumps.list".format(sn_type), dtype=None, skip_header=1)
    i = 0
    n = 0
    while (i < len(timedat) -n):
        if ((timedat['f0'])[i] == dirs[i]):
            i += 1
        #print "ok"
        else:
            #print i
            n += 1
            timedat = np.delete(timedat, i-1)
    
    times = timedat['f1']
    dz_times = (times*(1+z))/86400 #converting from seconds to days and multiplying by redshift
    names = timedat['f0']
    ind = min(range(len(times)), key=lambda q:abs(dz_times[q]-t))
    if ind == 0:
        AB = 400
    else:
        filename = names[ind]

        dat = np.loadtxt('../../Lanl.ccsn.spectra/{}/{}/spectra.out.6.1'.format(sn_type, filename))
        wavelength = dat[:,0]
        z_dat = dat
        z_dat[:,0] *= (1+z)
        z_wvl = z_dat[:,0]

        ind = min(range(len(z_wvl)), key=lambda q:abs(z_wvl[q]-filwv[0]))
        ind_2 = min(range(len(z_wvl)), key=lambda q:abs(z_wvl[q]-filwv[flen-1]))
        #print ind, ind_2
        tz_dat = z_dat[ind_2:ind]
        tz_width = tz_dat[:,2]

        lum = 0
        for x in tz_dat:
            wt = filt(x[0])
            lum += x[4]*wt
        Flux = lum/(4*np.pi*np.power(lum_dist,2))
        F_l = Flux/bandwidth #wavelength of the bandpass
        F_v = F_l*(np.power(1.63,2))/(beta)
        #print "F_v", F_v
        #this last number is in janskys
        AB = 2.5*(23 - np.log10(F_v)) - 48.60 - 0.45
        #print "AB", AB
    return AB
