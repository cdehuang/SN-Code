import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib
import numpy as np
import matplotlib.pylab as plt
from scipy.interpolate import UnivariateSpline
from scipy.integrate import quad
from StringIO import StringIO
import pyfits

#want to plot the number of SN for a given SFR

#star formation rate, depends on z
z_list = np.linspace(1, 10, 9)
m = len(z_list)
#SFR = np.zeros(m) + 1
SFR = np.zeros(m) + 2*z_list
SFR = SFR/(max(SFR))

#efficiency, also depends on z, but let's assume that it's flat for now
eta = np.zeros(m) + 2

#n is the number of clusters being simulated
#n=2
#voldat = np.genfromtxt("volumedata_a1689.dat", dtype=None, skip_header=n)

#by with seg we mean using a segmentation map, so smaller volume an vice-versa
#mus = voldat[:,0]
#withseg = voldat[:,2]
#noseg = voldat[:,1]

#need to do some 2D interpolation

#mulength = 33
#delta = 9/33.0
#X = np.arange(1, 10, delta)
#Y = np.arange(1, 60, 59/33.0)
#X, Y = np.meshgrid(X, Y)

volarray = volwmag(z_list[0])
Y = volarray[:,0]
X = z_list
vol = volarray[:,2]
vol = vol.reshape(len(vol),1)

for i in X[1:len(z_list)]:
    l = volwmag(i)
    l = (l[:,2]).reshape(len(l),1)
    vol = np.hstack((vol, l))

for k in range(len(SFR)):
    vol[:,k] = SFR[k]*vol[:,k]*eta[k]

Z = vol
Z = Z[0:5]
Y = Y[0:5]

#can try this interpolate thing later, for now it seems as though vol declines too quickly for it to want to do that
#f_z = interpolate.interp2d(X, Y, vol, kind ='cubic')
#xnew = np.linspace(1, 10, 0.5)
#ynew = np.linspace(1, 60, 1)
#znew = f_z(xnew, ynew)

plt.clf()
#plt.contour(xnew, ynew, znew)
plt.clabel(CS, inline=1, fontsize=10)
CS = plt.contour(X, Y, Z)
plt.xlabel("redshift")
plt.ylabel("$\mu$")
plt.title('N(z, $\mu$), for flat $\eta(z)$ and SFR')
plt.savefig("number_SN_lin.png")

plt.clf()
plt.plot(z_list, SFR)
plt.xlabel("redshift")
plt.ylabel("SFR")
plt.savefig("SFR_flat.png")
plt.show()



