import os
from os.path import join, getsize, expanduser
import re
import fileinput
import numpy as np
import matplotlib.pylab as plt
from scipy.interpolate import UnivariateSpline
from scipy.integrate import quad
import scipy
from StringIO import StringIO
import pyfits

#change so that it is a function that takes in a supernova type (z25B, etc.), a redshift, and possible other things, and then returns a light curve. 
#checked the bolometeric luminosity..it summed up on the order of 10^44 erg/s which seems about right to me. For the bandpass, it's 10-^39, which seems believable.

def shock(b_list):
    f_list = b_list[0:len(b_list)-1]
    i_list = b_list[1:len(b_list)]
    diff = i_list - f_list
    ind_list = np.asarray(np.where(diff<0))
    f_ind = ind_list[0, 0:len(ind_list[0])-1]
    i_ind = ind_list[0, 1:len(ind_list[0])]
    ind_diff = i_ind - f_ind
    i = 0
    p = 0
    while (i < 4) & (p < len(ind_diff)):
        if (ind_diff[p] == 1):
            i += 1
            p += 1
            print i
        else:
            i = 0
            p += 1
            print i
    ind = p-3
    #ind = np.argmax(diff<0) #need to make something a bit more sophisticated, that will    iterate through for a few times...
    #diff_cut = diff[ind+1:len(diff)]
    #ind2 = np.argmax(diff_cut<0)
    #if (ind2-ind) > 1:
        #ind = ind2
    idx = ind_list[0,ind]
    tb_list = b_list[idx:len(b_list)]
    return tb_list, idx

#sys.exit()

filenames = []
sn_type = '1a'
redshift = 0.0
stepsize = 1
l_d = 6.87604E4 #default is luminosity distance for z=7
#for root, dirs, files in os.walk('Lanl.ccsn.spectra/{}'.format(sn_type)):
#    for name in files:
#        filenames.append(os.path.join(root, name))
#fnames = [i for i in filenames if re.match('Lanl.ccsn.spectra/{}/................/spectra.out.6.1'.format(sn_type),i) is not None]
#fnames = [i for i in fnames if os.stat(i).st_size !=0]
#fnames = fnames[0::stepsize]
filter = "../test/WFC3_IR_F160W.dat"
filterdat = np.loadtxt(filter)
filwv = filterdat[:,0]
filamp = filterdat[:,1]
filt = scipy.interpolate.interp1d(filwv, filamp)
flen = len(filterdat)
bandwidth = 2683 #looked this up
#Assume a cosmology to calculate a distance. Using Ned's calculator for now
z = redshift
pc_cm = 3.08567758E18 
m_list = []
beta = 3E-13 #double checked with: http://www.stsci.edu/hst/nicmos/documents/handbooks/current_NEW/Appendix_B.14.3.html
i = 0
H_0 = 73 #not listed in the paper, but perhaps they used a different one? It doesn't start to look right until about 80.6
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
fnames = "hsiao_template/snflux_1a.dat"
lfname = "hsiao_template/lc_template.dat"

#need to do something about this huge for loop

dat = np.loadtxt(fnames)
#Redshift
z_dat = dat
z_dat[:,1] *= (1+z)
z_wvl = z_dat[:,1] #the wavelengths are in angstroms, same as the files from whalen's group
#trim to the appropriate size
z_times = (z_dat[:,0])*(1+z)
z_dat[:,0] *= (1+z)
days = np.linspace(-20, 85, 106)*(1+z)
m_list = []
for x in np.arange(len(days)):
    flux = 0
    f_vlist = []
    if x < (len(days)-2):
        idx = (np.abs(z_times - days[x])).argmin()
        idx_2 = (np.abs(z_times - days[x+1])).argmin() - 1
        spectra = z_dat[idx:idx_2, :]
        #spectra[0:2]
        wvlnth = spectra[:,1]
        ind = min(range(len(wvlnth)), key=lambda x:abs(wvlnth[x]-filwv[1]))
        ind_2 = min(range(len(wvlnth)), key=lambda x:abs(wvlnth[x]-filwv[flen-1]))
        spectra_dat = spectra[ind:ind_2]
        tz_width = spectra_dat[2,1] - spectra_dat[1,1]#the width of each segment
        width = spectra[2,1] - spectra[1,1]
        lum = 0
        for x in spectra_dat:
            seg = x[2]*filt(x[1])
            lam = x[1]/1E4
            f_v = seg*(np.power(lam,2))/(beta)
            flux += f_v
    else:
        idx = (np.abs(z_times - days[x])).argmin()
        idx_2 = len(z_times) - 1
        spectra = z_dat[idx:idx_2, :]
        wvlnth = spectra[:,1]
        ind = min(range(len(wvlnth)), key=lambda x:abs(wvlnth[x]-filwv[1]))
        ind_2 = min(range(len(wvlnth)), key=lambda x:abs(wvlnth[x]-filwv[flen-1]))
        spectra_dat = spectra[ind:ind_2]
        tz_width = spectra_dat[2,1] - spectra_dat[1,1] #the width of each segment
        width = spectra[2,1] - spectra[1,1]
        lum = 0
        for x in spectra_dat:
            seg = x[2]*filt(x[1])
            lam = x[1]/1E4
            f_v = seg*(np.power(lam,2))/(beta)
            flux += f_v
            f_vlist.append(f_v)
    #Flux = lum/(4*np.pi*np.power(lum_dist,2))
    AB = 2.5*(23-np.log10(flux)) - 48.60
    m_list.append(AB)

m_list = np.asarray(m_list)

lfspectra = np.loadtxt(lfname, skiprows=1)
epoch = lfspectra[:,0]
H_band = lfspectra[:,8]


#plot light curve
plt.clf()
#plt.scatter(dzt_times, tm_list, color='green')
plt.plot(days, m_list)
plt.plot(epoch, H_band, color='green')
plt.gca().invert_yaxis()
#plt.xlim(0,1000)
#plt.ylim(34, 26)
plt.xlabel("time (days)")
plt.ylabel("AB Magnitudes (1.63 um)")
plt.savefig("{}_{}_comparison2.png".format(sn_type, redshift))
plt.show()
