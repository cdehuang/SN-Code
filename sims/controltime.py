# chuang
# 2014.10.01
"""
Calculates the control time. Takes in lightcurves, cadence, magnification limits as input. 

Should eventually hook this up with something else that will then go on to calculate the curves for analytical solutions to how often we expect to detect various SNs.
"""
import numpy as np
import time
from scipy.optimize import fsolve
from scipy import interpolate
from scipy.integrate import quad
import matplotlib.pyplot as plt
import astropy.io.ascii as ascii
import volfun
import os

def mag_to_flux(mag):
    flux = np.power(10, 23 - (mag + 48.60)/2.5)
    return( flux )

def mumag(mag, mu):
    jsky = np.power(10, 23 -(mag+ 48.60)/2.5)
    mu_jsky = mu*jsky
    mu_mag = 2.5*(23 - np.log10(mu_jsky)) - 48.60
    return( mu_mag )

def dat_to_array(filename):
    dat = ascii.read(filename)
    dat = np.asarray(dat)
    datarray = dat.view(np.float).reshape(dat.shape + (-1,))
    return( datarray )

def get_indexes(x, list):
    index_list = []
    for i in range(len(list)):
        if list[i] == x:
            index_list.append(i)
    return( index_list )

def find_nearest_idx(array, value):
    idx = (np.abs(array-value)).argmin()
    return( idx )

#   In this function, the control time is lesser of either the interval over which the SN remains above a detectable threshold or the time between consecutive observations. Note that lc is actually a filename
#   This is really bizarre and bad code....
#   If the shock breakout period is less than a day, this rounds it up to a day so that we don't *completely* miss it.
def control_time(lc, mu, cadence=90, maglim=26.8, obswindow=700):
    lc_dat = dat_to_array(lc)
    days = lc_dat[:, 0]
    #print days
    mags = lc_dat[:, 1]
    mmags = np.asarray([mumag(i, mu) for i in mags]) - maglim
    #if the difference is negative, then it is brighter than the limit and we can see it. Also not sure how to find the intersect, so I'll just approximate the number of days.
    l = interpolate.interp1d(days, mmags)
    dmin = min(days)
    dmax = obswindow
    daysint = dmax - dmin
    x = np.linspace(dmin, dmax, daysint)
    points = l(x)
    #print points
    negidxs = np.where(points<0)
    #print negidxs
    negidxs = negidxs[0]
    if negidxs.size == 0:
        ct = 0.0
    else:
        negidxdiff = np.diff(negidxs)
        endposidx = np.where(negidxdiff>1)
        endposidx = endposidx[0]
        endidx = negidxdiff[endposidx]
        endidxs = negidxs[endposidx]+1
        endidxs = np.append(endidxs, negidxs[-1])
        range = endidx - 1
        begidxs = endidxs[0:-1] + range
        begidxs = np.append(negidxs[0], begidxs)
        begdays = x[begidxs]
        #print begdays
        enddays = x[endidxs]
        #print enddays
        segments = enddays-begdays
        total_days = np.sum(segments)
        if len(negidxs) == 1:
            ct = 1.0
        elif total_days < 90:
            ct = total_days
        elif total_days >= 90:
            ct = 90.0
    return( ct )

def highSFR(x):
    if (x < 12):
        val = 0.2
    else:
        val =  np.power(10, -0.25*x + 2.3)
    return val

def Kroupa(x):
    if x >= 1.:
        val = np.power(x, -7./3)
    elif (x < 1.) & (x >= 0.1):
        val = np.power(x, -4./3)
    else:
        val = np.power(x, -1/3.)
    return val

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

#Calculates the total number of CCSN (for the given z and mass), basically sums over the mus
def number_SN(z, mass, sfr=1, eta=1, duration=700):

    mu_min = 1

    if mu_min < 50:
        mu_max = 60
    else:
        mu_max = mu_min + 10

    efficiency = eta

    all = quad(Kroupa2, 0.1, 100)
    allsn = all[0]*efficiency
    print allsn

    if mass < 40:
        cc = quad(Kroupa, mass-5, mass+5)
    elif (mass == 40):
        cc = quad(Kroupa, mass-10, mass+10)
    fraction = cc[0]/all[0] #the fraction of all the stellar mass created that goes into the stars we're care about.
    ccsn = fraction*efficiency

    if sfr == 1:
        I = quad(highSFR, z, z+1)
    else:
        s = UnivariateSpline(sfr[:,0], sfr[:,1])
        I = quad(s, z, z+1)

    #   number of cc sn per Mpc^3 per redshift per our mass range
    numccsn = ccsn*I[0]
    print numccsn

    lc = "../output/lightcurve_z{}G_{}.dat".format(mass, z)

    if os.path.isfile('../output/volumemag_whalen_{}.dat'.format(z)):
        #if os.stat('../output/volumemag_whalen_{}.dat'.format(z)).st_size != 0:
        magvolumes = '../output/volumemag_whalen_{}.dat'.format(z)
        print "using the Whalen et al cosmology"
        vlist = ascii.read(magvolumes)
        vlist = np.asarray(vlist)
        volarray = vlist.view(np.float).reshape(vlist.shape + (-1,))
        volarray = volarray[mu_min: mu_max]
        #print volarray
        #print mu_min, mu_max
        clusters = vlist.dtype.names
    elif os.path.isfile('../output/volumemag_{}.dat'.format(z)):
        magvolumes = '../output/volumemag_{}.dat'.format(z)
        print "using the standard LCDM cosmology"
        vlist = ascii.read(magvolumes)
        vlist = np.asarray(vlist)
        volarray = vlist.view(np.float).reshape(vlist.shape + (-1,))
        volarray = volarray[mu_min-1: mu_max]
        #print mu_min, mu_max
        clusters = vlist.dtype.names
    else:
        print "calculating volumes..."
        volarray = volfun.volwmag(z, mumin=mu_min, mumax=mu_max)

    #   get the total volume at each magnification
    volumes = volarray[:, 1:volarray[0].size]
    volumes = np.sum(volumes, 1)
    volumes = np.append(np.abs(np.diff(volumes)), volumes[-1])

    initialnumsn = volumes*numccsn
    #print initialnumsn

    control_times = np.asarray([control_time(lc, i) for i in np.linspace(mu_min, mu_max-1, mu_max - mu_min)])

    #control_times = np.zeros(60) + 690.
    #control_times[0:10] = 0

    control_times_years = control_times/365.
    control_times_years_z = control_times_years/(1. + z)
    #print control_times_years_z

    number_times = len(control_times) - 1
    indexes = np.linspace(0, number_times - 2, number_times)
    indexes = indexes.astype(int)

    numberdetections = [initialnumsn[i]*control_times_years_z[i] for i in indexes]
    totalsn = np.sum(numberdetections)

    return( totalsn )

#   okay now it's time to get the total number of SN by integrating over all masses

def detections_per_redshift(zmin=5.0, zmax=12.0):
    
    redshifts = np.linspace(zmin, zmax, zmax-zmin+1)
    masses = [15, 25, 40]
    SN_array = []
    for z in redshifts:
        SN_array_z = []
        for i in masses:
            n = number_SN(z, i)
            print "redshift, mass, number", z,i,n
            SN_array_z.append(n)
        SN_array_z = np.asarray(SN_array_z)
        SN_array.append(SN_array_z)
    SN_per_redshift = np.sum(SN_array, 1)

    return( SN_per_redshift )

def plotting():
    start = time.time()
    zmin = 5.0
    zmax = 12.0
    zs = np.linspace(zmin, zmax, zmax-zmin+1)
    sn_z = detections_per_redshift(zmin=zmin, zmax=zmax)
    sn_zc = sn_z/25.
    total_sn = np.sum(sn_z)
    s = interpolate.interp1d(zs, sn_zc)
    end = time.time()
    print "Finished calculations in {} seconds".format(end-start)
    print "Total {} SN detected".format(total_sn)
    print sn_z
    plt.clf()
    plt.semilogy(zs, s(zs))
    plt.xlabel("Redshift [z]")
    plt.ylabel("Detections per Cluster")
    plt.savefig("../output/controltime.png")
    plt.show()



