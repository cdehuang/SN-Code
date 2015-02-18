"""
Calculation for the effective volumes accessible by galaxy clusters. 
L. Moustakas, JPL/Caltech
Begun 11 May 2012
After many hacks for the WeiOne JD1 paper and for a tailored calculation of Abell 1689,
generalizing for the volumes paper, September 2012. 
"""
#
# chuang
# 1.28.15
#
#
# the old version was getting to messy so this one is a bit cleaned up.
# this one seems to be working alright and it's able to accept a new parameter file, so I guess that's done then. 
#

import pyregion
import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate as si
from astropy.io import fits as pyfits
import os
import sys
import asciitable
from astropy.io import ascii
from astropy.table import Table
import cosmolopy as cp
import cosmolopy.distance as cd

# here are all the cosmology parameters

cdad = cd.angular_diameter_distance
cosmo = {'omega_M_0' : 0.3,
    'omega_lambda_0' : 0.7,
    'omega_k_0': 0.0,
    'h' : 0.7}
extras = {'w' : -1.0}
cosmo.update(extras)

#
# output is an array:
# first column is the magnification
# second column is the volume without the segementation map subtracted
# third column is the volume with the segmentation map subtracted
#
# NB: the volume in each column is the volume for that magnification and above. Not just at that magnification. This was done so that you could assign a magnification the entire volume (even if it was a lower limit).
#
#
# INPUT PARAMETERS:
#
# z: enter in the redshift for which you'd like to calculate the volumes
# mumin: the minimum magnification to calculate the volume. Default is 1 (no magnification)
# mumax: the maxmimum magnification to calculate the volume. Default is 100 (anything above 100 will thus be included in 100)
# parameter_file: this should give the location of all of the necessary files for this calculation.
# Columns of the file should read:
# 1-alias, 2-cluster directory, 3-deflection map x, 4-deflection map y, 5-redshift of the cluster, 6-redshift of the lensing model
# Default location is currently: '../"Debugging Things"/vp_lensmodels.lis'
#

def volwmag(z, mumin=1, mumax=100, parameter_file='../"Debugging Things"/vp_lensmodels.lis'):

    if (len(sys.argv)>1):
        lensmodels = sys.argv[1]
    else:
        lensmodels=parameter_file
        print lensmodels
    if not (os.path.isfile(lensmodels)):
        print 'missing startup file'
        sys.exit()

    clist = asciitable.read(lensmodels, names=['alias', 'clusterdir', 'deflxfile', 'deflyfile', 'zcluster', 'zmodel'])

    # the grid of magnifications
    mu = np.linspace(1, 100, 100)

    mspacing = mumax-mumin + 1
    magref = np.linspace(mumin, mumax, mspacing)

    TotalVolarrayWithSeg = np.zeros(mu.shape) 
    TotalVolarray = np.zeros(mu.shape)
    PartialVolarray = np.zeros((mu.size, clist.size))
    TotalSPVolArray = np.zeros(clist.shape)

    # target source redshift for rescaling model deflection maps
    ztarget = z

    # Plot stuff
    fig=plt.figure()
    MFig=fig.add_subplot(111)
    MFig.set_xlabel('Magnification [$\mu$]')
    MFig.set_ylabel('Effective Volume (>$\mu$) [Mpc$^3$]')
    MFig.set_xlim(0.2, 100)
    MFig.set_ylim(1.0, 2e5)
    MFig.set_xscale('log')
    MFig.set_yscale('log')

    # annotations -- check that these make sense
    MFig.text(30, 2e4, 'z=[9,10]', size=13)
    ytext = 1e5

    outdata = open('volumedata_a1689.dat','w')

    for ii in np.arange(clist.size):
        alias = clist['alias'][ii]
        clusterdir = clist['clusterdir'][ii]
        deflxfile = clist['deflxfile'][ii]
        deflyfile = clist['deflyfile'][ii]
        zcluster = clist['zcluster'][ii]
        zmodel = clist['zmodel'][ii]

    # note that the model redshifts are in names of Adi's files
    # cosmic SFR per year are for rest frame times

        rescale = \
            (cdad(ztarget, zcluster, **cosmo) /cdad(ztarget, **cosmo)) * \
            (cdad(zmodel, **cosmo) / cdad(zmodel, zcluster, **cosmo))

    # read in Adi's deflections
        axlist = pyfits.open('../{}/{}'.format(clusterdir, deflxfile))
        aylist = pyfits.open('../{}/{}'.format(clusterdir, deflyfile))
        ax, ay = rescale*axlist[0].data, rescale*aylist[0].data

    # read in segmentation map
        try:
            segfile = pyfits.open('../{}/{}'.format(clusterdir, segfile))
            segmap=segfile[0].data
            segflag = 1
        except:
            segflag = 0
            segmap=ax*0.0

    # do some initializations etc
        Unity = np.zeros(ax.shape)+1.0
        header = axlist[0].header
        try:
            cdelt1, cdelt2 = header['CDELT1'], header['CDELT2']
            header['CDELT1'], header['CDELT2'] = cdelt1, cdelt2
            pxsc = np.abs(cdelt1)*3600.0 # in arcseconds per linear pixel dimension
        except:
            cd1_1, cd2_2 = header['CD1_1'], header['CD2_2']
            pxsc = np.abs(cd1_1)*3600.0 # in arcseconds per linear pixel dimension

    # Calculate the magnification from Jacobian. Note that it first calculates the gradient with respect to the 'y' then with respect to the 'x' axis
        axy, axx = np.gradient(ax)
        ayy, ayx = np.gradient(ay)
        Jxx = Unity -axx
        Jxy =       -axy
        Jyx =       -ayx
        Jyy = Unity -ayy

        Mu    = Unity / (Jxx*Jyy - Jxy*Jyx)
        AbsMu = np.abs(Mu)

    # define the total wfc3ir outline mask
        regionfile = '../{}/{}_acs_wfc3ir_outline.reg'.format(clusterdir, clusterdir)
        rwcs=pyregion.open(regionfile)
        rcoo=pyregion.open(regionfile).as_imagecoord(header)
        maskboo=rwcs.get_mask(hdu=axlist[0])

    # diff_comoving_volume is differential comoving volume per unit redshift per unit solid angle, in units of Mpc**3 ster**-1. So convert this to a volume per (unlensed) pixel, for a source plane at ztarget.
        VolPix = (pxsc/206265.0)**2 * (cd.comoving_volume(ztarget+1, **cosmo) - cd.comoving_volume(ztarget, **cosmo))

    # The PixMap is now the Mpc**3 volume that each pixel corresponds to in the target source plane.
        PixMap = VolPix / AbsMu

    # Now calculate the source-plane volume for each of the lower-limit magnifications.
        TotalSPVol = np.sum(PixMap[((AbsMu>0) & (maskboo))])
        TotalSPVolArray[ii] = TotalSPVol
        VolarrayWithSeg = np.zeros(mu.shape)
        Volarray        = np.zeros(mu.shape)
        for jj in np.arange(mu.size):
            #        VolarrayWithSeg[jj] = np.sum( PixMap[((AbsMu>mu[jj]) & (rmask==1.0) & (segmap==0) )] )
            #        Volarray[jj]        = np.sum( PixMap[((AbsMu>mu[jj]) & (rmask==1.0))] )
            VolarrayWithSeg[jj] = np.sum( PixMap[((AbsMu>mu[jj]) & (maskboo) & (segmap==0) )] )
            Volarray[jj]        = np.sum( PixMap[((AbsMu>mu[jj]) & (maskboo))] )
        TotalVolarrayWithSeg += VolarrayWithSeg
        TotalVolarray += Volarray
        PartialVolarray[:,ii] = VolarrayWithSeg

        s=si.UnivariateSpline(np.log10(mu),np.log10(VolarrayWithSeg),s=5)
        # need to write the mu, Volarray data to a file.... 
        MFig.plot(mu,Volarray,label=alias)
        MFig.plot(mu,VolarrayWithSeg,'--')
        MFig.legend(loc=3,prop={'size':8})

    # called 'PartialMVolArray' because it doesn't give the total volume over all the clusters, but rather the volume in each individual cluster for that particular magnification and above. 
    PartialMVolArray = np.zeros((magref.size, clist.size))
    for x in np.arange(len(PartialVolarray[0])):
        sx = si.UnivariateSpline(np.log10(mu), np.log10(PartialVolarray[:,x]), s=5)
        partialvol = (np.asarray(10**sx(np.log10(magref))))
        PartialMVolArray[:,x] = partialvol
    # no error checking here yet...
    sws = si.UnivariateSpline(np.log10(mu),np.log10(TotalVolarrayWithSeg),s=5)
    s   = si.UnivariateSpline(np.log10(mu),np.log10(TotalVolarray),s=5)
    volwithseg = (np.asarray(10**sws(np.log10(magref)))).reshape(len(magref), 1)
    volnoseg = (np.asarray(10**s(np.log10(magref)))).reshape(len(magref),1)
    magref = (np.asarray(magref)).reshape(len(magref),1)
    MVolArray = np.hstack((magref, volnoseg, volwithseg))
    PartialMVolArray = np.c_[magref, PartialMVolArray]

    return PartialMVolArray

