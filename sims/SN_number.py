import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib
import numpy as np
import matplotlib.pylab as plt
from scipy.interpolate import UnivariateSpline
from scipy.integrate import quad
from StringIO import StringIO
import pyfits
import volfun

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx

def SFRslope(x):
    #corresponds to the optimistic case in the Whalen paper
    if (x < 12):
        val = 0.2
    else:
        val =  np.power(10, -0.25*x + 2.3)
    return val

def Salpeter(x):
    return np.power(x, -2.35)

def Kroupa(x):
    if x >= 0.5:
        val = np.power(x, -2.3)
    elif (x < 0.5) & (x >= 0.08):
        val = np.power(x, -1.3)
    else:
        val = np.power(x, -0,3)
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

#Basically the Kroupa, but multiplied by the mass of the star as well so instead of number of stars you get the amount of stellar mass contained in some range of stars
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

def KroupaTest(x): #for finding the fraction of the mass that goes into SN, not used in the actual code
    val = np.power(x, -4./3)
    return val

#gives you a luminosity distance in Mpc, consider changing this to cosmolopy for more flexibility 
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

#returns the number of SN expected at the given z for each mu. First column is the mu. mass is mass of star in solar masses
def number_SN(z, mass, mu_min=1, eta=0.01):
    #setting the parameters to go into the volfun code 
    if mu_min < 50:
        mu_max = 60
    else:
        mu_max = mu_min + 10

    efficiency = eta

    #IMF--Salpeter IMF
    #all = quad(KroupaSalpeter, 0.1, 400)
    #note that this may or may not be the right way to handle it but
    all = quad(Kroupa2, 0.1, 100)
    allsn = all[0]*efficiency
    cc = quad(Kroupa2, mass-5, mass+5)
    fraction = cc[0]/all[0] #the fraction of all the stellar mass created that goes into the stars we're care about.
    ccsn = fraction*efficiency

    #the actual amount of stellar mass per Mpc^3
    I = quad(SFRslope, z, z+1)
    numccsn = ccsn*I[0] #number of cc sn per Mpc^3 using their numbers it comes out to like 0.00024 SN per Mpc^3, sounds believable

    #for testing purposes, take out the dependence on Leonidas' code and see if we are still creating SNs. Just integrate the area
    volarray = volfun.volwmag(z, mumin=mu_min, mumax=mu_max) #make sure this is the shell that you think it is...
    volarray[:,2] = volarray[:,2]*numccsn
    numarray = np.delete(volarray, 1, 1)
 
    return numarray




