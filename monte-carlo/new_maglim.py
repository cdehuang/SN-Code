import os
from os.path import join, getsize, expanduser
import re
import fileinput
import numpy as np
import matplotlib.pylab as plt
from scipy.interpolate import UnivariateSpline
from scipy.integrate import quad
from StringIO import StringIO
import pyfits
import scipy

#   If you input a light curve, it figures out the minimum magnification needed for any part of it to be visible.
#   Tested and it looks like it is functioning correctly.
def magfactor(lc, lim):
    peak = min(lc)
    if lim < min:
        jsky = np.power(10, 23 -(peak + 48.60)/2.5)
        #print jsky
        lim = np.power(10, 23 -(lim + 48.60)/2.5)
        #print lim
        mu = lim/jsky
    else:
        mu = 0
    return mu

#   Should remove the shock breakout section of the light curve.
def remove_shock_breakout(mags):
    
    mag_diff = np.diff(mags)
    
    #   this is where the SN is getting brighter. Split into arrays where the light curve is continuously brightening.
    brighter = np.where(mag_diff < 0)[0]
    splitarray = np.array_split(brighter, np.where(np.diff(brighter)!=1)[0]+1)
    
    #   find where the lightcurve starts to rise for a while
    for i in splitarray:
        lengt = len(i)
        if lengt > 2:
            firstarr = i
            firstidx = i[0]
            break
    
    #   We could just remove everything that comes before the beginning index, but I think I'll just remove everything that comes before the dip instead. This may not be the best way to do it though...
    fidx = np.where(mags == np.max(mags[0:firstidx]))[0]
    if fidx == 0:
        fidx = firstidx
    
    #   Return the trimmed array
    return mags[fidx:len(mags)-1], fidx

#   Should return the lightcurve for the SN in a given redshift, SN_type, and filter
#   Note though that it currently only returns the magnitudes and not the dates, so it's not very useful (and you can't just straight up plot it and get something that makes sense)
def whalen_lightcurve(sn_type, redshift, steps=2, fil='../../filters/WFC3_IR_F160W.dat', bandwidth=2683):
    #indent everything to go under this function
    #sn_type = 'z15g'
    #redshift = 7
    #steps = 5
    #fil = 'WFC3_IR_F160W'
    filenames = []
    for root, dirs, files in os.walk('../../Lanl.ccsn.spectra/{}'.format(sn_type)):
        for name in files:
            filenames.append(os.path.join(root, name))
    fnames = [i for i in filenames if re.match('../../Lanl.ccsn.spectra/{}/................/spectra.out.6.1'.format(sn_type),i) is not None]
    fnames = [i for i in fnames if os.stat(i).st_size !=0]
    #uses only 1/10th of the files, need to modify the times also if you want to change this
    p = steps
    fnames = fnames[0::p]
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
    H_0 = 70 #not listed in the paper, but perhaps they used a different one? It doesn't start to look right until about 80.6
    omega_M = 0.308
    omega_L = 0.692
    omega_K = 0
    pc_cm = 3.08567758E18
    c = 299792.458 #km/s
    D_H = c/H_0
    I = lambda z1: 1/np.sqrt(omega_M*(1+z1)**3 + omega_K*(1+z1)**2 + omega_L)
    Int, err = quad(I, 0, z)
    D_m = D_H*Int
    D_a = D_m/(1+z)
    D_l = ((1 + z)**2)*D_a
    lum_dist = D_l*1E6*pc_cm
    
    flux_list = []
    m_list = []
    beta = 3E-13 #double checked with: http://www.stsci.edu/hst/nicmos/documents/handbooks/current_NEW/Appendix_B.14.3.html
    i = 0
    for y in fnames:
        dat = np.loadtxt(y)
        wavelength = dat[:,0]
        
        #Redshift
        z_dat = dat
        z_dat[:,0] *= (1+z)
        z_wvl = z_dat[:,0]
        
        #trim to the appropriate size
        ind = min(range(len(z_wvl)), key=lambda q:abs(z_wvl[q]-filwv[1]))
        ind_2 = min(range(len(z_wvl)), key=lambda q:abs(z_wvl[q]-filwv[flen-1]))
        #print ind, ind_2
        tz_dat = z_dat[ind_2:ind]
        tz_width = tz_dat[:,2]
        
        #integrate by summing
        #flux = 0
        lum = 0
        for x in tz_dat:
            wt = filt(x[0])
            lum += x[4]*wt
        #seg_lum = x[5]*wt*(1./(1+z))
        #seg_flux = seg_lum/(4*np.pi*np.power(lum_dist, 2))
        #lam = x[0]/1E4
        #f_v = seg_flux*(np.power(lam,2))/(beta)
        #flux += f_v
        
        #flux_list.append(flux)
        Flux = lum/(4*np.pi*np.power(lum_dist, 2))
        F_l = Flux/bandwidth
        F_v = F_l*(np.power(1.63, 2))/beta
        #AB = 2.5*(23-np.log10(flux))-48.60 -0.45
        AB = 2.5*(23 - np.log10(F_v)) - 48.60 #zeropoint #25.96
        m_list.append(AB)
    
    m_list = np.asarray(m_list)
    #return m_list
    tm_list, ind = remove_shock_breakout(m_list)
    return tm_list

if __name__ =="__main__":
    mag_diff = np.diff(mags)
    
    #   this is where the SN is getting brighter. Split into arrays where the light curve is continuously brightening.
    brighter = np.where(mag_diff < 0)[0]
    splitarray = np.array_split(brighter, np.where(np.diff(brighter)!=1)[0]+1)
    
    #   find where the lightcurve starts to rise for a while
    for i in splitarray:
        lengt = len(i)
        if lengt > 2:
            firstarr = i
            firstidx = i[0]
            break
    
    #   We could just remove everything that comes before the beginning index, but I think I'll just remove everything that comes before the dip instead. This may not be the best way to do it though...
    fidx = np.where(mags == np.max(mags[0:firstidx]))[0]
    if fidx == 0:
        fidx = firstidx
    
    #   Return the trimmed array
    #return mags[fidx:len(mags)-1], fidx
