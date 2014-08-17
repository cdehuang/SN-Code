# S.Rodney
# 2014.08.15
"""
Reformat a pile of spectra from Dan Whalen to make them accessible 
for modeling light curves with Kyle Barbary's sncosmo package.

NOTE: the SED files generated with this program have units of flux density
  (in erg/sec/cm^2/Angstrom) through a shell of radius 10 pc.  Thus, when you
  load these as a model in sncosmo and use amplitude=1,  the bandmag() function
  will deliver absolute magnitudes for the specified band.
"""
import os
from astropy.io import ascii
from glob import glob
import time
import numpy as np

# surface area of the 'SPECTRUM' radiative transfer grid used by Whalen et al
# for computing luminosities  (r=10^18 cm)
# A_GRIDSURFACE  = 4 * np.pi * 1e36   # cm^2

# surface area of a shell with radius = 10 pc  (r=3.08567758e19 cm)
A_10pc = 4 * np.pi * (3.08567758e19)**2   # cm^2


def mksedfile( modelname, rename=True, fixtypos=True, combine=True,
               minphasestep=0.5,  Nwaveskip=10 ):
    start = time.time()
    if rename :
        print("Renaming spectra.out files to more useful file names.")
        filelist, phaselist = renamedir( modelname, minphasestep=minphasestep )
    else :
        filelist = glob( '%s/%s*.dat'%(modelname,modelname))

    if fixtypos :
        print("Fixing E-100 typos in %i files"%len(filelist))
        filelist = fixtypos_list( filelist )
    else :
        filelist = glob( '%s/%s*.dat'%(modelname,modelname))

    if combine :
        print("Constructing an SED file from %i files"%len(filelist))
        sedfile = writesedfile( modelname, filelist, Nwaveskip=Nwaveskip )
    end = time.time()

    print( "    Total time = %i sec"%( end-start ) )
    return( sedfile )

def renamedir( modelname, minphasestep=0.5 ):
    """ Rename and relocate the spectra.out files to give them names that encode
    the phase relative to explosion.
    We only process files that are separated in time by at least <minphasestep>
    days from the previous file.
    """
    start = time.time()

    dumplist = os.path.abspath( os.path.join( modelname, 'spectrum_dumps.list' ) )
    dumpdat = ascii.read( dumplist, data_start=1 )
    dirnamelist = dumpdat['col1']
    tsec = dumpdat['col2']
    phase = tsec / 3600. / 24.

    filelist = []
    phaselist = []
    lastphase = -minphasestep
    for i in range( len(dirnamelist) ):
        subdirname = os.path.join(modelname, dirnamelist[i] )
        oldnamelist = glob( subdirname+'/spectra.out*')
        if len(oldnamelist) and (phase[i]-lastphase) >= minphasestep :
            oldname = oldnamelist[0]
            newname = os.path.join(modelname, modelname.strip('/') + '_%07.3f.dat'%phase[i])
            os.rename( oldname, newname )
            filelist.append( newname )
            phaselist.append( phase[i])
            lastphase = phase[i]
            print( newname )
        # remove the old subdirectory and all contents
        for oldfile in oldnamelist :
            if os.path.exists( oldfile) : os.remove( oldfile )
        if os.path.isdir( subdirname ) :
            os.removedirs( subdirname )
    end = time.time()
    print( "    Renamed %i files in %i seconds. (%.3e sec per file)"%(len(filelist), (end-start), (end-start)/len(filelist)))
    return(filelist, phaselist)

def fixtypos_list( filelist ):
    """  fix the E-100 typo on a list of files
    """
    start = time.time()
    for file in filelist:
        fixtypo1( file )
    end = time.time()
    print( "    Typo fixing :  %i seconds elapsed (%.3e sec per file)"%((end-start), (end-start)/len(filelist)))
    return( filelist )


def fixtypo1( filename ) :
    """  fix the E-100 typo in a single file
    """
    import re
    fin = open(filename,'r')
    contents = fin.read()
    fin.close()
    contents = re.sub(r'(?<!E)-(?=\d{3})', 'E-', contents)
    fout = open(filename,'w')
    fout.write( contents )
    fout.close()

def trimlist( filelist, Nwaveskip=5 ):
    """
    Trim all files in filelist using trim1, keeping only every Nth wavelength
    """
    start = time.time()
    for filename in filelist :
        trim1( filename, Nwaveskip=Nwaveskip )
    end = time.time()
    print( "    Renamed %i files in %i seconds. (%.3e sec per file)"%(len(filelist), (end-start), (end-start)/len(filelist)))
    return(filelist)


def trim1( filename, Nwaveskip=5, wavemax=20000 ):
    """
    Extract the wavelength and flambda values from a spectra.out file,
    keeping only every Nth wavelength step.
    :param filename:
    :return:
    """
    fin = open(filename,'r')
    linelist = fin.readlines()[1:]
    fin.close()
    wave = np.array([ float(line.split()[0]) for line in linelist ])
    flux = np.array([ float(line.split()[-1]) for line in linelist ])[wave.argsort()]
    wave.sort()
    itrim = np.arange(0,len(np.where(wave<wavemax)[0]),Nwaveskip)
    #fout = open(filename+'.fasttrim','w')
    # for w,f in zip( wave[itrim], flux[itrim] ):
    #    print>>fout,'%15.5f    %15.6e'%(w,f)
    #fout.close()
    return( wave[itrim], flux[itrim] )

def writesedfile( modelname, inputfilelist, outsedfile=None, Nwaveskip=5 ):
    """ Write out a single SED file for the given <modelname>.  The output
    is an ascii text file containing an informative header and three columns
    giving PHASE (days since explosion), WAVELENGTH (angstroms) and
    FLAMBDA (i.e. flux density, in erg/s/cm^2/Angstrom).
    """
    if outsedfile is None :
        outsedfile = 'popIII-%s.sed.restframe10pc.dat'%modelname

    if modelname.startswith('z') : progenitor ='Red Supergiant'
    elif modelname.startswith('u') : progenitor ='Blue Supergiant'
    else : progenitor = 'unknown'

    try :
        mass = int(modelname[1:3])
        ecode = modelname[-1]
    except :
        mass, ecode = 0, None

    if ecode == 'B' : energy = 0.6
    elif ecode == 'D' : energy = 1.2
    elif ecode == 'G' : energy = 2.4
    else : energy = 0

    header = """### %s ########################
#
# Pop III Core Collapse SN model from
#   Whalen, D.J. et al. 2013 ApJ, 768, 95
#
#  MODEL NAME : %s
#  PROGENITOR TYPE: %s
#  PROGENITOR MASS: %i Msun
#  EXPLOSION ENERGY: %.1f x 10^51 erg
#
#  PHASE = rest-frame days since explosion
#  LAMBDA_REST = rest-frame wavelength in angstroms
#  F_LAMBDA = rest-frame flux density in erg/s/cm^2/Angstrom
#             through a shell of radius 10pc
#
# Trimmed and reformatted by S.Rodney on 2014.08.15
#
# PHASE  LAMBDA_REST    F_LAMBDA"""%( outsedfile, modelname, progenitor, mass, energy )

    phaselist = np.array( [ float(filename.split('_')[1][:7]) for filename in inputfilelist ] )
    fout = open(outsedfile, 'w')
    print >> fout, header
    for filename,phase in zip( inputfilelist, phaselist ) :
        wave,flux = trim1( filename, Nwaveskip=Nwaveskip, wavemax=20000.)
        fluxdensity = flux / A_10pc
        phasearray = np.zeros( len(wave) ) + phase
        for p,w,f in zip( phasearray, wave, fluxdensity):
            print >> fout, '%7.3f  %10.3f  %12.8f'%(p,w,f)
    fout.close()
    return( outsedfile )



