#!/usr/bin/env python
"""
Calculation for the effective volume as a function of lower
magnification limit for CLASH clusters.
L. Moustakas, JPL/Caltech, 11
May 2012 CLASH discrete, do not circulate outside collaboration.
"""
import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate as si
import pyfits
import os
import asciitable
import cosmolopy as cp
import cosmolopy.distance as cd
from matplotlib import cm
import pyregion

cdad = cd.angular_diameter_distance
cosmo = {'omega_M_0' : 0.3,
         'omega_lambda_0' : 0.7,
         'omega_k_0': 0.0,
         'h' : 0.70}
extras = {'w' : -1.0}
cosmo.update(extras)

#This function returns an array that gives you the volume at various z, the first column tells you what those z's were

def volz(minmu, zmin=2.5, zmax=11.5):
    clist=asciitable.read('lensmodels.lis',
                          names=['alias','clusterdir','deflxfile','deflyfile','segfile', 'zcluster','zmodel'])

    # Observed F160W AB magnitude is 25.7. Its restframe B-band magnitude
    # is approximately -22.4, uncorrected for the magnification of ~15.5.
    # So the corrected absolute magnitude is -22.4+2.5*log(15.5) = -19.4.

    # target source redshift for rescaling model deflection maps
    spacing = zmax - zmin + 1
    ztarget = np.linspace(zmin,zmax,spacing)

    # a reasonable grid of magnifications
    mu=np.linspace(0.2,120,100)
    # fiducial threshold magnification
    mufiducial = minmu
    mustr = '%.1f' % mufiducial

    # set the anchor point of Mref <-> muref
    Mref, muref = -19.5, 9.0
    #Mref, muref = -19.4, 15.5
    Mlim = Mref+2.5*np.log10(mu/muref)
    def getmu(Mag):
        return muref*10**(0.4*(Mag-Mref))

    TotalVolarray = np.zeros(mu.shape)
    TotalRedarray = np.zeros(ztarget.shape)

    # Plot stuff
    fig=plt.figure()
    MFig=fig.add_subplot(111)
    MFig.set_xlabel('Redshift [z]')
    MFig.set_ylabel('Effective Volume ($\mu>$'+mustr+') [Mpc$^3$]')
    MFig.set_xlim(2.0,12)
    MFig.set_ylim(1.0,2e5)
    MFig.set_yscale('log')

    # # some annotations
    # MFig.text(30,2e4,'z=[9,10]',size=13)
    # ytext = 1e5
    # MFig.text(0.22,ytext,'Unlensed M$_{B}$ [AB] limit:',size=10)
    # plotmus=[-18.0,Mref,-21.0,-22.0]
    # for pp in plotmus: MFig.text(getmu(pp)/1.2,ytext,'%.1f' % pp,size=10)

    # outdata = open('redshiftdata.dat','w')

    for ii in np.arange(clist.size):

        alias      = clist['alias'][ii]
        clusterdir = clist['clusterdir'][ii]
        deflxfile  = clist['deflxfile'][ii]
        deflyfile  = clist['deflyfile'][ii]
        zcluster   = clist['zcluster'][ii]
        zmodel     = clist['zmodel'][ii]

    # read in the deflections from Adi's lens models 
        axlist = pyfits.open(clusterdir+'/'+deflxfile)
        aylist = pyfits.open(clusterdir+'/'+deflyfile)
    #    ax, ay = rescale*axlist[0].data, rescale*aylist[0].data
        ax, ay = axlist[0].data, aylist[0].data

    # do some initializations etc
        Unity  = np.zeros(ax.shape)+1.0
        header = axlist[0].header
        try:
            cdelt1, cdelt2 = header['CDELT1'], header['CDELT2']
            header['CDELT1'], header['CDELT2'] = cdelt1, cdelt2
            pxsc = np.abs(cdelt1)*3600.0 # in arcseconds per linear pixel dimension
        except:
            cd1_1, cd2_2 = header['CD1_1'], header['CD2_2']
            pxsc = np.abs(cd1_1)*3600.0 # in arcseconds per linear pixel dimension

        regionfile = clusterdir+'/'+'acs_outline.reg'
        rwcs=pyregion.open(regionfile)
        rcoo=pyregion.open(regionfile).as_imagecoord(header)
        maskboo=rwcs.get_mask(hdu=axlist[0])

        Redarray=np.zeros(ztarget.shape)
        for zi in np.arange(ztarget.size):
            #print alias, ztarget[zi]
            rescale = \
                (cdad(ztarget[zi], zcluster, **cosmo) / cdad(ztarget[zi], **cosmo)) * \
                (cdad(zmodel, **cosmo) / cdad(zmodel, zcluster, **cosmo))
            
    # Calculate the magnification from the Jacobian.  Note that it first
    # calculates the gradient with respect to the 'y' axis, and then wrt
    # the 'x' axis. 
            axy, axx = np.gradient(rescale*ax)
            ayy, ayx = np.gradient(rescale*ay)
            Jxx = Unity -axx
            Jxy =       -axy
            Jyx =       -ayx
            Jyy = Unity -ayy
        
            Mu    = Unity / (Jxx*Jyy - Jxy*Jyx)
            AbsMu = np.abs(Mu)

    # The diff_comoving_volume is the differential comoving volume per
    # unit redshift per unit solid angle, in units of Mpc**3 ster**-1.  So
    # we convert that to the Volume per (unlensed) pixel, for a source
    # plane at z=9-10.
            VolPix = (pxsc/206265.0)**2 * cd.diff_comoving_volume(ztarget[zi], **cosmo)

    # PixMap will now be the Mpc**3 volume that each pixel corresponds to
    # in the z=9-10 source plane.
            
            PixMap = VolPix / AbsMu 

            Redarray[zi]=np.sum(PixMap[((AbsMu>mufiducial) & (maskboo))])

        TotalRedarray += Redarray
        MFig.plot(ztarget,Redarray,label=alias)
        MFig.legend(loc=3,prop={'size':8})

    MFig.plot(ztarget,TotalRedarray,linewidth=4,color='b',label='CLASH 12')
    MFig.legend(loc=3,prop={'size':10})
    
    d1 = len(ztarget)
    VolZarray = ztarget.reshape(d1, 1)
    TotalRedarray = TotalRedarray.reshape(d1, 1)
    VolZarray = np.hstack((VolZarray, TotalRedarray))

    return VolZarray

minmu = np.linspace(1, 19, 10)

volarray = volz(minmu[0])
zarray = volarray[:,0]
vols = volarray[:,1]
vols = vols.reshape(len(vols), 1)

for i in minmu[1:len(minmu)]:
    l = volz(i)
    l = (l[:,1]).reshape(len(l),1)
    vols = np.hstack((vols, l))

plt.clf()
for x in range(len(minmu)):
    plt.plot(zarray, vols[:,x], color=cm.jet(1.*x/len(minmu)), label="$\mu$ = {}".format(minmu[x]))
plt.legend(loc=3,prop={'size':10})
plt.xlabel('Redshift [z]')
plt.ylabel('Effective Volume ($\mu>$'+mustr+') [Mpc$^3$]')
plt.xlim(2.0,12)
plt.ylim(1.0,2e5)
plt.yscale('log')
plt.savefig("vol_vs_redshift.png")
plt.show()




