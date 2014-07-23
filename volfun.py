#!/usr/bin/env python
"""
Calculation for the effective volumes accessible by galaxy clusters. 
L. Moustakas, JPL/Caltech
Begun 11 May 2012
After many hacks for the WeiOne JD1 paper and for a tailored calculation of Abell 1689,
generalizing for the volumes paper, September 2012. 
"""
import pyregion 
import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate as si
import pyfits
import os
import sys
import asciitable
import cosmolopy as cp
import cosmolopy.distance as cd

cdad = cd.angular_diameter_distance
cosmo = {'omega_M_0' : 0.3,
         'omega_lambda_0' : 0.7,
         'omega_k_0': 0.0,
         'h' : 0.7}
extras = {'w' : -1.0}
cosmo.update(extras)

def volwmag(z):
    #if __name__ == "__main__":
    if (len(sys.argv)>1):
        lensmodels = sys.argv[1]
    else:
        #    lensmodels='lensmodels.lis'
        lensmodels='vp_lensmodels.lis'
    if not (os.path.isfile(lensmodels)):
        print 'missing startup file'
        sys.exit()

    clist=asciitable.read(lensmodels,
                          names=['alias','clusterdir','deflxfile','deflyfile','segfile','zcluster','zmodel'])

    # Observed F160W AB magnitude is 25.7. Its restframe B-band magnitude
    # is approximately -22.4, uncorrected for the magnification of ~15.5.
    # So the corrected absolute magnitude is -22.4+2.5*log(15.5) = -19.4.

    # a reasonable grid of magnifications
    #mu=np.linspace(0.2,120,100)
    mu=np.linspace(1,100,100)

    # representative magnifications 
    ###magref=[5.4, 8.0, 13.5, 14.5, 18.7, 57.7]
    magref=[1, 2, 3, 4, 5, 5.4, 6, 6.6, 7, 7.5, 8.0, 9, 10, 11, 12, 13, 13.5, 14, 14.5, 15, 16, 17, 18, 18.7, 20, 25, 30, 35, 40, 45, 50, 57.7, 60]
    # set the anchor point of Mref <-> muref
    ###  Mref, muref = -19.5, 14.5
    Mref, muref = -19.5, 9.0
    Mlim = Mref+2.5*np.log10(mu/muref)
    def getmu(Mag):
        return muref*10**(0.4*(Mag-Mref))

    TotalVolarrayWithSeg = np.zeros(mu.shape)
    TotalVolarray        = np.zeros(mu.shape)

    # target source redshift for rescaling model deflection maps
    # *** 
    #ztarget = 9.5
    #ztarget = 9.6
    #ztarget = 5.5
    ztarget = z

    # Plot stuff
    fig=plt.figure()
    MFig=fig.add_subplot(111)
    MFig.set_xlabel('Magnification [$\mu$]')
    MFig.set_ylabel('Effective Volume (>$\mu$) [Mpc$^3$]')
    MFig.set_xlim(0.2,100)
    MFig.set_ylim(1.0,2e5)
    MFig.set_xscale('log')
    MFig.set_yscale('log')

    # some annotations
    # *** 
    MFig.text(30,2e4,'z=[9,10]',size=13)
    #MFig.text(30,2e4,'z=[5,6]',size=13) 
    ytext = 1e5
    # *** 
    #MFig.text(0.22,ytext,'Unlensed M$_{B}$ AB limit:',size=10)
    ## including some absolute magnitudes that correspond to magnifications 
    #plotmus=[-18.0,Mref,-21.0,-22.0]
    #for pp in plotmus: MFig.text(getmu(pp)/1.2,ytext,'%.1f' % pp,size=10)

    outdata = open('volumedata_a1689.dat','w')

    for ii in np.arange(clist.size):
    #for ii in [0]: 

        alias      = clist['alias'][ii]
        clusterdir = clist['clusterdir'][ii]
        deflxfile  = clist['deflxfile'][ii]
        deflyfile  = clist['deflyfile'][ii]
        segfile    = clist['segfile'][ii]
        zcluster   = clist['zcluster'][ii]
        zmodel     = clist['zmodel'][ii]

        rescale = \
            (cdad(ztarget, zcluster, **cosmo) / cdad(ztarget, **cosmo)) * \
            (cdad(zmodel, **cosmo) / cdad(zmodel, zcluster, **cosmo))
        
    # read in the deflections from Adi's lens models 
        axlist = pyfits.open(clusterdir+'/'+deflxfile)
        aylist = pyfits.open(clusterdir+'/'+deflyfile)
        ax, ay = rescale*axlist[0].data, rescale*aylist[0].data

    # read in the segmentation map, which we are implicitly assuming has
    # already been wregistered to the model
        try:
            segfile = pyfits.open(clusterdir+'/'+segfile)
            segmap=segfile[0].data
            segflag = 1
        except:
            segflag = 0
            segmap=ax*0.0 

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

    # Calculate the magnification from the Jacobian.  Note that it first
    # calculates the gradient with respect to the 'y' axis, and then wrt
    # the 'x' axis. 
        axy, axx = np.gradient(ax)
        ayy, ayx = np.gradient(ay)
        Jxx = Unity -axx
        Jxy =       -axy
        Jyx =       -ayx
        Jyy = Unity -ayy

        Mu    = Unity / (Jxx*Jyy - Jxy*Jyx)
        AbsMu = np.abs(Mu)

    # define the total wfc3ir outline mask
    #    regionfile = clusterdir+'/'+'wfc3ir_outline.reg'
        regionfile = clusterdir+'/'+'acs_outline.reg'
        rwcs=pyregion.open(regionfile)
        rcoo=pyregion.open(regionfile).as_imagecoord(header)
        maskboo=rwcs.get_mask(hdu=axlist[0])
        #    rmask=np.zeros(ax.shape)
        #    rmask[(maskboo)]=1.0

    # The diff_comoving_volume is the differential comoving volume per
    # unit redshift per unit solid angle, in units of Mpc**3 ster**-1.  So
    # we convert that to the Volume per (unlensed) pixel, for a source
    # plane at ztarget.
        VolPix = (pxsc/206265.0)**2 * cd.diff_comoving_volume(ztarget, **cosmo)

    # PixMap will now be the Mpc**3 volume that each pixel corresponds to
    # in the z=9-10 source plane.
        PixMap = VolPix / AbsMu 

    # Now let's zap the areas of the mosaic that are covered by objects
    # via the IR-detected segmentation map.
    # /Volumes/Archive/CLASH/archive.stsci.edu/pub/clash/outgoing/macs1149/HST/catalogs/mosaicdrizzle_image_pipeline/IR_detection/SExtractor

    # Now let us calculate the source-plane volume for each of the
    # magnification lower limits we established.
        VolarrayWithSeg = np.zeros(mu.shape)
        Volarray        = np.zeros(mu.shape)
        for jj in np.arange(mu.size):
            #        VolarrayWithSeg[jj] = np.sum( PixMap[((AbsMu>mu[jj]) & (rmask==1.0) & (segmap==0) )] )
            #        Volarray[jj]        = np.sum( PixMap[((AbsMu>mu[jj]) & (rmask==1.0))] )
            VolarrayWithSeg[jj] = np.sum( PixMap[((AbsMu>mu[jj]) & (maskboo) & (segmap==0) )] )
            Volarray[jj]        = np.sum( PixMap[((AbsMu>mu[jj]) & (maskboo))] )
        TotalVolarrayWithSeg += VolarrayWithSeg
        TotalVolarray += Volarray

        s=si.UnivariateSpline(np.log10(mu),np.log10(VolarrayWithSeg),s=5)
        print alias,clusterdir,10**s(np.log10(muref))
        outstring = '%s %s %.1f ' % (alias,clusterdir,10**s(np.log10(muref)))
        outdata.write(outstring+'\n')
        # need to write the mu, Volarray data to a file.... 
        MFig.plot(mu,Volarray,label=alias)
        MFig.plot(mu,VolarrayWithSeg,'--')
        MFig.legend(loc=3,prop={'size':8})

    # no error checking here yet...
    sws = si.UnivariateSpline(np.log10(mu),np.log10(TotalVolarrayWithSeg),s=5)
    s   = si.UnivariateSpline(np.log10(mu),np.log10(TotalVolarray),s=5)
    volwithseg = (np.asarray(10**sws(np.log10(magref)))).reshape(len(magref), 1)
    volnoseg = (np.asarray(10**s(np.log10(magref)))).reshape(len(magref),1)
    magref = (np.asarray(magref)).reshape(len(magref),1)
    MVolArray = np.hstack((magref, volnoseg, volwithseg))

    return MVolArray


# # in pyraf:
# ### --> wregister alphaX_11495_z9p7LAST.fits detectionImage_SEGM.fits alphaX_11495_z9p7LAST_segmatch.fits fluxconserve- boundary=constant constant=0.0
# """
# wregister detectionImage_SEGM.fits alphaX_11495_z9p7LAST.fits detectionImage_SEGM_defmatch.fits fluxconserve- boundary=constant constant=0.0
# """
# # ok, magnification maps generated from deflection maps that have been wregistered to the segmentation map size are bedunked.  perhaps this isn't surprising, but still.
# # julian's maps are super coarse.  not going to use those right now.
# # so, two choices for exploring. first, i could generate the magnification map, and match that to the higher res seg map. second, i could simply match the seg map to the deflection map with conserve=no and 
# 
