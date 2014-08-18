__author__ = 'rodney'
"""
2014.08.17
S.Rodney

Demo for using sncosmo and sncosmohst to generate
pop III SN light curves from the Whalen models.

This assumes you have installed sncosmo (and all prerequisites) following
the instructions here :
   http://sncosmo.readthedocs.org/en/latest/install.html

and that you have cloned SR's github repo for sncosmost and put it on
your python path :
   git clone git@github.com:srodney/sncosmost.git

"""
import numpy as np
from matplotlib import pyplot as pl
from matplotlib import ticker
import sncosmo
from sncosmost import hstbandpasses, pop3models
from astropy.cosmology import FlatLambdaCDM

cosmo = FlatLambdaCDM( 70., 0.3 )

def lcfig( z=7 ):
    """  Generate a light curve figure that mimics Whalen et al 2013
    :return:
    """

    for modelname in ['z15G','z25G','z40G'] :

        # initialize a supernova model :
        snmodel = sncosmo.Model(source=modelname)

        # Fix the redshift for this instantiation of the model
        # (NOTE: this does not apply any cosmological dimming. It only
        #   shifts the wavelengths)
        snmodel.set( z = z )

        # generate the H band light curve
        tobs = np.arange( 0, 1000, 10 )
        M160 = snmodel.bandmag( 'f160w', 'ab', tobs ) # Absolute Magnitude
        m160 = M160 + cosmo.distmod( z ).value # apparent magnitude

        pl.plot( tobs, m160 )
    ax = pl.gca()
    ax.invert_yaxis()
    ax.set_xlabel('Time (observer-frame days)')
    ax.set_ylabel('Apparent Magnitude in F160W')
    ax.set_xlim( 0, 1000 )
    ax.set_ylim( 36, 28 )

    pl.draw()

