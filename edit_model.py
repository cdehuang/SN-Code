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
sn_type = 'z40g'
redshift = 7
stepsize = 1
l_d = 6.87604E4 #default is luminosity distance for z=7
for root, dirs, files in os.walk('Lanl.ccsn.spectra/{}'.format(sn_type)):
    for name in files:
        filenames.append(os.path.join(root, name))
fnames = [i for i in filenames if re.match('Lanl.ccsn.spectra/{}/................/spectra.out.6.1'.format(sn_type),i) is not None]
fnames = [i for i in fnames if os.stat(i).st_size !=0]
fnames = fnames[0::stepsize]
filter = "test/WFC3_IR_F160W.dat" 
filterdat = np.loadtxt(filter)
filwv = filterdat[:,0]
filamp = filterdat[:,1]
filt = scipy.interpolate.UnivariateSpline(filwv, filamp) 
flen = len(filterdat)
bandwidth = 2683 #looked this up
#Assume a cosmology to calculate a distance. Using Ned's calculator for now
z = redshift
H_0 = 69.6 
omega_M = 0.308
omega_L = 0.692
lum_dist = l_d #this is in Mpc
lpc_dist = lum_dist*1E6 #just changed to pc
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
        #val = min(range(len(filwv)), key = lambda i:abs(filwv[i] - x[0]))
        #print val
        #seg = x[4]*filamp[val]
        wt = filt(x[0])
        seg = x[4]*wt
        lum += seg
    #print lum

#this number was in erg s^-1, so now we want to get it into erg sec^-1 cm^-2 A^-1
    Flux = lum/(4*np.pi*np.power(pc_cm*lpc_dist,2))
    #print "Flux", Flux
    F_l = Flux/bandwidth 
    #print "F_l", F_l
    F_v = F_l*(np.power(1.63,2))/(beta)
    #print "F_v", F_v
#this last number is in jansky
    AB = 2.5*(23-np.log10(F_v)) - 48.60 #fixed this, sort of. I thought the original formula was wrong, but could not figure out why. I must have misread the wikipedia page. Still off by about 4 orders of magnitude 
    #print "AB", AB
    m_list.append(AB)

m_list = np.asarray(m_list)
tm_list, ind = shock(m_list)

#fixes errors in the spectrum_dumps.list without making a new s_d.list file
dirs = os.listdir('LANL.CCSN.SPECTRA/{}/'.format(sn_type))
dirs = [i for i in dirs if re.match('....-.-dmp......',i) is not None]
dirs = [i for i in dirs if os.stat('LANL.CCSN.SPECTRA/{}/{}/spectra.out.6.1'.format(sn_type,i)).st_size != 0] #corrects for the blank files in z40b and z40g

timedat = np.genfromtxt("LANL.CCSN.SPECTRA/{}/spectrum_dumps.list".format(sn_type), dtype=None, skip_header=1)
i = 0
n = 0
while (i < len(timedat) -n):
    if ((timedat['f0'])[i] == dirs[i]):
        i += 1
        #print "ok"
    else:
        print i
        n += 1
        timedat = np.delete(timedat, i-1)

times = timedat['f1']
times = times[0::stepsize]
t_times = times[ind:len(times)]
dzt_times = (t_times*(1+z))/86400
#dz_times = (times*(1+z))/86400 #converting from seconds to days and multiplying by redshift
names = timedat['f0']
            
#m_dat = np.asarray(m_list)
m_dat = np.asarray(tm_list)
#dz_dat = np.asarray(dz_times)
dz_dat = np.asarray(dzt_times)
print "first times:", m_dat[0:20], dz_dat[0:20]
#pyfits.writeto("{}_{}.fits".format(sn_type, redshift), m_dat)
#pyfits.append("{}_{}.fits".format(sn_type, redshift), dz_dat)

#plot light curve
plt.clf()
plt.scatter(dzt_times, tm_list)
#plt.scatter(dz_times, m_list)
plt.gca().invert_yaxis()
plt.xlim(0,1000)
plt.ylim(34, 26)
plt.xlabel("time (days)")
plt.ylabel("AB Magnitudes (1.63 um)")
plt.savefig("{}_{}_shock_scatter.png".format(sn_type, redshift))
plt.show()
