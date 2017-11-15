import os
import numbers
import sys
import numpy as np
import scipy.interpolate

import transphere_python.natconst as nc


def ratranRun(r=0.0, rho=0.0, temp=0.0, db=0.0, abund=0.0, vr=0.0, tdust=0.0, dustonly=0, 
              file='transphere.mdl', dpc=0.0, imsize=129, pixel=0.5, channels=50, 
              channel_width=0.1, trans='220.0e9', writeonly=0, skyonly=1, molfile='', 
              outputfile="ratranResult", unit='yada'):
    """

    Parameters
    ----------
    r : array of floats
        Radii (in cm) for computation
    rho : array of floats
        Mass density (in g cm-3) at each radius
    temp : array of floats
        Temperature (in K) at each radius
    db : array of floats or float
        Doppler b-parameter linewidth (km s-1). Specifically, the 1/e half-width of the 
        line profile.
    abund : float or array of floats
        Number abundance (molecules per H2 molecule) of species in question
    vr : array of floats or float
        Radial velocity (km s-1)
    tdust : float
        Dust temperature (K)
    dustonly : float (0 or nonzero)
        Flag to set dust temperature `tdust` at `temp` if `dustonly` equals 0.0.
    file : str
        Name of input transphere model `.mdl` file.
    dpc : float 
        Distance (pc) to simulated source
    imsize : int
        Number of pixels in the image. Always square.
    pixel : float
        Scale in arcseconds of each pixel.
    trans : str
        Transition numbers. See online documentation for further explanation.
        If the input file has level populations defined, trans contains the transition numbers
        to be calculated. These are the numbers at the start of lines (10+NLEV) to (10+NLEV+NLIN) 
        in the molecular data file. For linear molecules with singlet sigma ground states without 
        hyperfine structure such as CO, CS and HCO+, trans is equal to Jup. 
        If there are no populations defined in the input model, then a pure continuum calculation 
        is done, and 'trans' contains the frequencies (in Hz).
    writeonly : int (1 or not 1)
        Flag: whether to output .inp files, etc. If not 1, then it won't.
    skyonly : int (0 or nonzero)
        Flag: whether to run AMC. If not 0, then AMC won't run.
    molfile : str
        Name of the molecular data file (usually '`molname`.dat').
        Often from the LAMDA Leiden database.
        Sometimes an '@' symbol (e.g. as in 'h13cn@xpol.dat') will crash Sky. 
        Rename the .dat file in such cases.
    outputfile : str
        Prefix of the output file of the ratran code.
    unit : str, 'K' or 'Jypx' or 'Wm2Hzsr'
        Units of the output image.

    Information
    -----------
    See this page: https://personal.sron.nl/~vdtak/ratran/running.html#amc-model
    for more information on the model parameters.

    """
    if tdust == 0.0:
        tdust = temp

    # Dust density * gas2dust
    nh = rho*100.0/nc.muh2/nc.mp

    print('Writing model in '+file)

    r1 = np.insert(r[0:-1], 0, 0)
    r2 = np.array(r)

    rr = np.zeros((len(r)), float)
    rr[1:] = 10**((np.log10(r1[1:])+np.log10(r2[1:]))/2.0)
    rr[0] = rr[1]

    nhf = scipy.interpolate.interp1d(np.log10(r2), np.log10(nh))
    tkf = scipy.interpolate.interp1d(np.log10(r2), np.log10(temp))
    tdf = scipy.interpolate.interp1d(np.log10(r2), np.log10(tdust))
    abf = scipy.interpolate.interp1d(np.log10(r2), np.log10(abund))

    nhint = 10**nhf(np.log10(rr))
    tkint = 10**tkf(np.log10(rr))
    tdint = 10**tdf(np.log10(rr))
    abint = 10**abf(np.log10(rr))

    nhint[0] = 0.0

    if isinstance(db, numbers.Number):
        db_arr = np.ones_like(rr) * db
    else:
        db_arr = db

    if isinstance(vr, numbers.Number):
        vr_arr = np.ones_like(rr) * vr
    else:
        vr_arr = vr

    # Check consistency with input grid here!!!

    f = open(file, 'w')
    f.write('# Ratran input file based on Transphere results'+'\n')
    if skyonly == 1:
        f.write('# ... intended for (SKY) continuum calculations only.'+'\n')
    f.write('rmax=%9.3E' % (max(r2)/1.0e2)+'\n')
    f.write('ncell=%i' % (len(r2))+'\n')
    f.write('tcmb=2.728\n')
    f.write('columns=id,ra,rb,nh,nm,tk,td,db,vr\n')
    f.write('gas:dust=100\n')
    if skyonly == 1:
        f.write('kappa=jena,thin,e6\n')
    f.write('@\n')
    for ii in range(0, len(r1)):
        # remember ra/r1 and rb/r2 should be in  m (NOT cm)
        f.write("%4i %12.5E %12.5E %12.5E %12.5E %12.5E %12.5E %12.5E %12.5E" % (
            ii+1, r1[ii]/100.0, r2[ii]/100.0, nhint[ii], nhint[ii]*abint[ii], tkint[ii], tdint[ii], db_arr[ii], vr_arr[ii])+'\n')
    f.close()

    if writeonly != 1:
        if skyonly == 0:
            if molfile == '':
                sys.exit(
                    'Error: for AMC calculations the molecular datafile (molfile) needs to be set.')
            f = open("amc.inp", 'w')
            f.write("source="+file+'\n')
            f.write("outfile=populations.pop\n")
            f.write("molfile="+molfile+'\n')
            f.write("snr=20\n")
            f.write("nphot=100\n")
            f.write("kappa=jena,thin,e6\n")
            f.write("seed=1971\n")
            f.write("go\n")
            f.write("q\n")
            f.write(" \n")
            print("Starting AMC calculation...")
            f.close()
            os.system('amc amc.inp')

        f = open("sky.inp", 'w')
        if skyonly == 0:
            f.write("source=populations.pop\n")
        else:
            f.write("source="+file+"\n")
        f.write("format=fits\n")
        f.write("outfile="+outputfile+"\n")
        f.write("trans="+trans+"\n")
        f.write("pix="+str(imsize)+","+str(pixel)+",32,2\n")
        if skyonly == 0:
            f.write("chan="+str(channels)+","+str(channel_width)+"\n")
        else:
            f.write("chan=1,1.0\n")
        f.write("distance="+str(dpc)+"\n")
        if skyonly == 0 and unit != 'jypx':
            f.write("units=K\n")
        else:
            f.write("units=jypx\n")
        f.write("go\n")
        f.write("q\n")
        f.close()

        print("Starting SKY calculation...")
        os.system("sky sky.inp")
