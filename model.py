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

filenames = []
for root, dirs, files in os.walk('Lanl.ccsn.spectra/z25B'):
    for name in files:
        filenames.append(os.path.join(root, name))
fnames = [i for i in filenames if re.match('Lanl.ccsn.spectra/z25B/................/spectra.out.6.1',i) is not None]
filter = "test/WFC3_IR_F160W.dat" 
filterdat = np.loadtxt(filter)
filwv = filterdat[:,0]
filamp = filterdat[:,1]
flen = len(filterdat)
#Assume a cosmology to calculate a distance. Using Ned's calculator for now
z = 7
H_0 = 69.6 
omega_M = 0.308
omega_L = 0.692
lum_dist = 6.87604E4 #this is in Mpc
lpc_dist = lum_dist*10E6 #just changed to pc
pc_cm = 3.08567758E18 
m_list = []
beta = 3E-13 #double checked with: http://www.stsci.edu/hst/nicmos/documents/handbooks/current_NEW/Appendix_B.14.3.html
i = 0

#need to do something about this huge for loop
for x in fnames:
    dat = np.loadtxt(x)
    wavelength = dat[:,0]

    #Redshift
    z_dat = dat
    z_dat[:,0] *= (1+z)
    z_wvl = z_dat[:,0]

    #trim to the appropriate size
    ind = min(range(len(z_wvl)), key=lambda x:abs(z_wvl[x]-filwv[0]))
    ind_2 = min(range(len(z_wvl)), key=lambda x:abs(z_wvl[x]-filwv[flen-1]))
    #print ind, ind_2
    tz_dat = z_dat[ind_2:ind]
    tz_width = tz_dat[:,2]

    #integrate by summing--consider checking with their total bolometric luminosity, try checking the absolute magnitude of the source, somewhere in the range of
    #Type 1a - 19.5, so this is probably -20 or brighter, check the redshift behavior, try a fixed section around the peak fo the spectrum
    #k-correction extend filter to check
    lum = 0
    for x in tz_dat:
        val = min(range(len(filwv)), key = lambda i:abs(filwv[i] - x[0]))
        #print val
        seg = x[4]*filamp[val]
        lum += seg
    print lum

#this number was in erg s^-1, so now we want to get it into erg sec^-1 cm^-2 A^-1
    Flux = lum/(4*np.pi*np.power(pc_cm*lpc_dist,2))
    print "Flux", Flux
    F_l = Flux/16300 #wavelength of the bandpass
    print "F_l", F_l
    F_v = F_l*(np.power(1.63,2))/(beta)
    print "F_v", F_v
#this last number is in erg s^-1
    AB = -2.5*np.log10(F_v) - 48.60
    print "AB", AB
    m_list.append(AB)

#multiply the times by (1+z) also
timedat = np.genfromtxt("LANL.CCSN.SPECTRA/z25B/s_d.list", dtype=None)
times = timedat['f1']
dz_times = (times*(1+z))/86400

#check that there are no errors in the spectrum_dumps.list
i = 0
dirs = os.listdir('LANL.CCSN.SPECTRA/z25B/')
dirs = dirs[4:len(dirs)]
names = timedat['f0']
while (i < len(names)):
    if names[i] == dirs[i]:
	i += 1
    else:
        timedat = np.delete(timedat, i)
            
m_dat = np.asarray(m_list)
dz_dat = np.asarray(dz_times)
pyfits.writeto("z25B_7.fits", m_dat)
pyfits.append("z25B_7.fits", dz_dat)

#plot light curve
plt.plot(dz_times, m_list)
plt.gca().invert_yaxis()
plt.xlim(0,1000)
#plt.ylim(42, 38)
plt.xlabel("time (days)")
plt.ylabel("AB Magnitudes (1.63 um)")
#plt.savefig("Z25B")
plt.show()
