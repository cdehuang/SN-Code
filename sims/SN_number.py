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

#returns the number of SN expected at the given z for each mu. First column is the mu.
def number_SN(z, mu_min=1, mu_max=100):
    #star formation rate, depends on z
    z_list = np.linspace(1, 30, 30)
    m = len(z_list)
    #SFR = np.zeros(m) + 1
    #SFR = np.zeros(m) + 2*z_list
    idx = find_nearest(z_list, 2)
    SFR = np.zeros(m)
    SFR[0:idx] = 0.5
    SFR[idx:m] = -2*z_list[idx:m]+ 60
    SFR[idx:m] = SFR[idx:m]/max(np.abs(SFR))

    #efficiency, also depends on z, but let's assume that it's flat for now
    idx = find_nearest(z_list, 5)
    eta = np.zeros(m)
    eta[idx:m] = z_list[idx:m]*0.75
    eta = eta/max(eta)
    
    volarray = volwmag(z_list[0], mumin=mu_min, mumax=mu_max)
    volarray[:,2] = volarray[:,2]*SFR[z]*eta[z]
    volarray = np.delete(volarray, 1, 1)
    #number_arr = np.c_[ number_arr, [vol*SFR[z]*eta[z]]]

    return volarray






