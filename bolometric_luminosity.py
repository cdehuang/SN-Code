#some more robust code to make a bolometric luminosity plot with the same scaling as: http://arxiv.org/pdf/1209.5459v4.pdf
import re
import fileinput
import numpy as np
import matplotlib.pylab as plt
import scipy
from scipy.interpolate import UnivariateSpline
from scipy.integrate import quad
from StringIO import StringIO
import pyfits

#indent everything to go under this function
filenames = []
sn_type = 'z40b'
redshift = 7
l_d = 6.87604E4
for root, dirs, files in os.walk('Lanl.ccsn.spectra/{}'.format(sn_type)):
    for name in files:
        filenames.append(os.path.join(root, name))
fnames = [i for i in filenames if re.match('Lanl.ccsn.spectra/{}/................/spectra.out.6.1'.format(sn_type),i) is not None]
fnames = [i for i in fnames if os.stat(i).st_size !=0]
    #Assume a cosmology to calculate a distance. Using Ned's calculator for now
z = 7
H_0 = 69.6 
omega_M = 0.308
omega_L = 0.692
lum_dist = l_d #this is in Mpc
lpc_dist = lum_dist*10E6 #just changed to pc
pc_cm = 3.08567758E18 
m_list = []
beta = 3E-13 #double checked with: http://www.stsci.edu/hst/nicmos/documents/handbooks/current_NEW/Appendix_B.14.3.html
i = 0
lum_list=[]
for y in fnames:
    #deal with empty files
    dat = np.loadtxt(y)
    wavelength = dat[:,0]

        #Redshift
    z_dat = dat
    z_dat[:,0] *= (1+z)
    z_wvl = z_dat[:,0]

        #integrate by summing checked this by replacing my original function weights with 1
    lum = np.sum(z_dat[:,4])
    lum_list.append(lum)

#fixes errors in the spectrum_dumps.list without making a new s_d.list file
dirs = os.listdir('LANL.CCSN.SPECTRA/{}/'.format(sn_type))
dirs = [i for i in dirs if re.match('....-.-dmp......',i) is not None]
dirs = [i for i in dirs if os.stat('LANL.CCSN.SPECTRA/{}/{}/spectra.out.6.1'.format(sn_type,i)).st_size != 0]
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
z_times = times*(1+z)
dz_times = (times*(1+z))/86400 #converting from seconds to days and multiplying by redshift
names = timedat['f0']
          
m_dat = np.asarray(m_list)
z_dat = np.asarray(z_times)
lum_dat = np.asarray(lum_list)
#pyfits.writeto("{}_{}_bol.fits".format(sn_type, redshift), m_dat)
pyfits.append("{}_{}_bol.fits".format(sn_type, redshift), z_dat)
pyfits.append("{}_{}_bol.fits".format(sn_type, redshift), lum_dat) #luminosities in erg/s

    #plot light curves
print len(dz_times)
print len(lum_dat)
plt.clf()
plt.loglog(z_times, lum_dat)
plt.xlim(10E3,10E6)
plt.ylim(10E39, 10E46)
plt.xlabel("time (s)")
plt.ylabel("Bolometric Luminosity (erg/s)")
plt.title("{}, z={}, bolometric light curves".format(sn_type, redshift))
plt.savefig("{}_{}_bolo.png".format(sn_type, redshift))
plt.show()
