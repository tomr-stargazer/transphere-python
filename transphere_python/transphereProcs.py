#
# Module to read-in dust opacity files
# This one adapted from readopac.pro
#  --- Jes Jorgensen (jeskj@nbi); Dec 2011
#

import transphere_python.natconst as nc


def readopac(nr=-1):
    import numpy as np
    import sys

    if nr == -1:
        nr = 1
    filename = 'dustopac_'+str(nr)+'.inp'
    print("Reading "+filename)
    f = open(filename, 'r')
    nf, ns = f.readline().strip().split()
    nf, ns = int(nf), int(ns)
    f.readline()
    if nf >= 1:
        cabs = np.zeros((nf, ns), float)
        csca = np.zeros((nf, ns), float)
        dum = 0.e0
        print(ns)
        for kk in range(0, nf):
            for iss in range(0, ns):
                dum = float(f.readline().strip())
                cabs[kk, iss] = dum

        for kk in range(0, nf):
            for iss in range(0, ns):
                dum = float(f.readline().strip())
                csca[kk, iss] = dum

        nrtrange = 1
        trange = [0.e0, 0.e0]

    else:
        nf = ns
        ismooth = 0
        nrtrange = 0
        line = f.readline().strip()
        line = f.readline().strip()
        ns, ismooth, nrtrange = line.split()
        ns, ismooth, nrtrange = int(ns), int(ismooth), int(nrtrange)
        if ismooth != 0:
            sys.exit('Error: Smoothing not yet allowed.')
        cabs = np.zeros((nf, ns), float)
        csca = np.zeros((nf, ns, nrtrange), float)
        dum = 0.e0
        trange = np.zeros(nrtrange+1, float)

        for ir in range(0, nrtrange):
            a, b = f.readline().strip().split()
            a, b = int(a), int(b)
            trange[ir] = b
            for kk in range(0, nf):
                for iss in range(0, ns):
                    dum = float(f.readline().split())
                    cabs[kk, iss, ir] = dum

            for kk in range(0, nf):
                for iss in range(0, ns):
                    dum = float(f.readline().split())
                    csca[kk, iss, ir] = dum

    f.close

# dd=findfile(file,count=count)
# if(count le 0) then begin
#    print,("Could not find frequency.inp. Taking frequency.dat")
#    file='frequency.dat'
#    dd=findfile(file,count=count)
#    if(count le 0) then begin
#       print,("Could not find frequency.dat either")
#       stop
#    endif
# endif

    file = 'frequency.inp'
    f = open(file, 'r')
    nf = int(f.readline().strip())
 #   if nnf != nf: sys.exit("ERROR: frequency file has different nr of points as dustopac file")
    dum = f.readline()
    freq = np.zeros(nf, float)
    wave = np.zeros(nf, float)
    for kk in range(0, nf):
        dum = float(f.readline().strip())
        freq[kk] = dum
        wave[kk] = 2.9979e14 / dum
    f.close

    opacity = {'ns': ns, 'nf': nf, 'freq': freq, 'wave': wave,
               'cabs': cabs, 'csca': csca, 'nrt': nrtrange, 'trange': trange}

    return opacity

#
# Module to write actual opacity as function of R
# Adapted from problem_session.pro
#  --- Jes Jorgensen (jeskj@nbi); Dec 2011
#


def writeopac(opacity, localdust, nr):
    import os

    if localdust == 1:
        os.system('rm -f dustopac_1.inp')
        f = open('dustopac.inp', 'w')
        f.write('1               Format number of this file'+'\n')
        f.write('1               Nr of dust species'+'\n')
        f.write(
            '============================================================================'+'\n')
        f.write('-1              Way in which this dust species is read (-1=file)'+'\n')
        f.write('0               Extension of name of dustopac_***.inp file'+'\n')
        f.write(
            '----------------------------------------------------------------------------'+'\n')
        f.close
        f = open('dustopac_0.inp', 'w')
        f.write(nr)
        f.write(str(opacity['nf'])+' 1\n')
        f.write(' ')
        redux = 1.e0
        for ir in range(0, nr):
            for inu in range(0, nr):
                f.write(opacity['cabs'][inu]*redux)
            for inu in range(0, opacity['nf']):
                f.write(opacity['csca'][inu]*redux)
            f.write(' ')
        f.close
    else:
        os.system('rm -f dustopac_0.inp')
        f = open('dustopac.inp', 'w')
        f.write('1               Format number of this file'+'\n')
        f.write('1               Nr of dust species'+'\n')
        f.write(
            '============================================================================'+'\n')
        f.write('-1              Way in which this dust species is read (-1=file)'+'\n')
        f.write('1               Extension of name of dustopac_***.inp file'+'\n')
        f.write(
            '----------------------------------------------------------------------------'+'\n')
        f.close()

        f = open('dustopac_1.inp', 'w')
        f.write(str(opacity['nf'])+' 1\n')
        f.write(' '+'\n')
        for inu in range(0, opacity['nf']):
            f.write(str(opacity['cabs'][inu][0])+'\n')
        for inu in range(0, opacity['nf']):
            f.write(str(opacity['csca'][inu][0])+'\n')
        f.close

#
# Find kappa at 0.55 micron
# Adapted from problem_session.pro
#  --- Jes Jorgensen (jeskj@nbi); Dec 2011
#


def findKappa(localdust, opacity):
    import numpy as np
    if localdust == 1:
        f = open('dustopac_0.inp', 'r')
        nr = int(f.readline().strip())
        nf, idum = readlines().strip().split()
        nf, idum = int(nf), int(idum)
#        if nnr != nr: sys.exit('ERROR when reading opacity file in findKappa (nnr != nr)')
#        if nnf != nf: sys.exit('ERROR when reading opacity file in findKappa (nnf != nf)')
        kap = np.zeros((nf, 2, nr), float)
        for kk in range(0, nr):
            for jj in range(0, 2):
                for ii in range(0, nf):
                    kap[ii, jj, kk] = float(f.readline().strip())
        f.close()
        ivis = 0
        fvis = 1e4 * nc.cc / 0.55e0
        for inu in range(1, nf):
            if (o.freq[inu]-fvis)*(o.freq[inu-1]-fvis) < 0.e0:
                ivis = inu
        eps = (fvis-o.freq[ivis-1]) / (o.freq[ivis]-o.freq[ivis-1])
        kappa = (1.e0-eps)*kap[ivis-1, 0, :]+eps * \
            kap[ivis, 0, :]  # Removed transpose
    else:
        f = open('dustopac_1.inp', 'r')
        nf, idum = f.readline().strip().split()
        nf, idum = int(nf), int(idum)
        #       if nnf != nf: sys.exit('ERROR when reading opacity file in findKappa (nnf != nf)')
        f.readline()
        kap = np.zeros((nf, 2), float)
        for jj in range(0, 2):
            for ii in range(0, nf):
                kap[ii, jj] = float(f.readline().strip())
        f.close()
        ivis = 0
        fvis = 1e4 * nc.cc / 0.55e0
        for inu in range(1, nf-1):
            if (opacity['freq'][inu]-fvis)*(opacity['freq'][inu-1]-fvis) < 0.e0:
                ivis = inu
        eps = (fvis-opacity['freq'][ivis-1]) / \
            (opacity['freq'][ivis]-opacity['freq'][ivis-1])
        kappa = (1.e0-eps)*kap[ivis-1, 0]+eps*kap[ivis, 0]

    return kappa

#
# Make radial grid
# Adapted from problem_session.pro
#  --- Jes Jorgensen (jeskj@nbi); Dec 2011
#


def makeRadialGrid(nref, rin, rout, rref, nr):
    import math
    import numpy as np

    if nref == 0:
        r = rin * (rout/rin)**(np.arange(nr)/(nr-1.e0))  # Simple Log grid
    else:
        #
        # Log grid with refinement at inner edge
        #
        lgr0 = math.log10(rin)
        lgr1 = math.log10(rout)
        lgrr = math.log10(rref)
        n1 = nref
        n2 = nr-nref
        n = np.arange(nr)-n1
        c = (lgr1-lgrr)/(1.e0*n2)
        b = c/(lgrr-lgr0)
        r = np.zeros((nr), float)
        r[0:n1] = (lgrr-lgr0)*math.e**(b*n[0:n1])
        r[n1:n1+n2] = (lgrr-lgr0)+c*n[n1:n1+n2]
        r = r-(r[0]+(r[0]-r[1]))
        r = r*(lgr1-lgr0)/(max(r))+lgr0
        r = 10.0**r
    return r

#
# Write Transphere input files
# Adapted from problem_session.pro
#  --- Jes Jorgensen (jeskj@nbi); Dec 2011
#


def writeTransphereInput(model):
    import natconst as nc
    import math
    import astroProcs
    import numpy as np
# Transphere input file
    f = open('transphere.inp', 'w')
    f.write(str(2)+'\n')
    f.write(str(model['nriter'])+'\n')
    f.write(str(model['convcrit'])+'\n')
    f.write(str(model['ncst'])+'\n')
    f.write(str(model['ncex'])+'\n')
    f.write(str(model['ncnr'])+'\n')
    f.write(str(model['itypemw'])+'\n')
    f.write(str(model['idump'])+'\n')
    f.close()
#
# Make the stellar information file
# (mstar and tstar are irrelevant; they are there for historical reasons)
#
    f = open('starinfo.inp', 'w')
    f.write(str(1)+'\n')
    f.write(str(model['rstar'])+'\n')
    f.write(str(model['mstar'])+'\n')
    f.write(str(model['tstar'])+'\n')
    f.close()
#
# The stellar spectrum
#
    f = open('starspectrum.inp', 'w')
    f.write(str(len(model['freq']))+'\n')
    sspec = (model['rstar']/nc.pc)**2*math.pi * \
        astroProcs.bplanck(model['freq'], model['tstar'])
    for inu in range(0, len(model['freq'])):
        f.write(str(model['freq'][inu])+' '+str(sspec[inu])+'\n')
    f.close()
    #
    # The exterior spectrum
    #
    if model['tbg'] == 0.0 or model['isrf'] == 0:
        bgspec = np.zeros((len(model['freq'])), float)
    elif model['tbg'] == -1:
        f = open('isrf.inp', 'r')
        nf = int(f.readline().strip())
        bgspec = np.zeros((len(model['freq'])), float)
        for ii in range(0, nf):
            bgspec[ii] = float(f.readline().strip())*model['isrf']
    else:
        if tbg > 0:
            bgspec = astroProcs.bplanck(
                model['freq'], model['tbg'])*model['isrf']

    f = open('external_meanint.inp', 'w')
    f.write(str(len(model['freq']))+'\n')
    for inu in range(0, len(model['freq'])):
        f.write(str(model['freq'][inu])+' '+str(bgspec[inu])+'\n')
    f.close()
#
# Write the envelope structure
#
    f = open('envstruct.inp', 'w')
    f.write(str(len(model['r']))+'\n')
    f.write(' '+'\n')
    for ir in range(0, len(model['r'])):
        f.write("%13.6E %13.6E %13.6E" % (model['r'][ir], model[
                'rho'][ir], 0.e0)+'\n')  # ,format='(3(E13.6,1X))'
    f.close()


def readEnvstruct(ext=0):
    import numpy as np

    if n_elements(ext) == 0:
        ext = ''
    file = 'envstruct'+ext+'.dat'
    f = open(file, 'r')
    nr = int(f.readline().strip())
    dat = np.zeros((3, nr), float)
    for ii in range(0, 3):
        for jj in range(0, nr):
            dum = float(f.readline().strip())
            dat[ii, jj] = dum
    f.close()

    r = dat[0, :]
    rho = dat[1, :]
    temp = dat[2, :]

    class envstruct:
        r = r
        rho = rho
        temp = temp

    return envstruct


def readConvhist():
    import numpy as np
    f = open('convhist.info', 'r')
    nn = int(f.readline().strip().split()[0])
    f.close()

    f = open('convhist.dat', 'r')
    nr = int(f.readline().strip())
    dat = np.zeros((9, nn, nr), float)
    for jj in range(0, nn):
        for kk in range(0, nr):
            dum = f.readline().strip().split()
            if dum == []:
                dum = f.readline().strip().split()
            dat[0:9, jj, kk] = np.array(dum, dtype=float)
    f.close()

# if nn gt 1 then idx=[1,2,0] else idx=[1,0]. Note transpose commands not
# executed...
    temp = dat[0, :, :]
    jjme = dat[1, :, :]
    hhme = dat[2, :, :]
    jj = dat[3, :, :]
    hh = dat[4, :, :]
    kapt = dat[5, :, :]
    kapj = dat[6, :, :]
    kaph = dat[7, :, :]
    fj = dat[8, :, :]

    f = open('envstruct.inp')
    nr = int(f.readline().strip())
    dat = np.zeros((3, nr), float)
    for ii in range(0, nr):
        dum = f.readline().strip().split()
        if dum == []:
            dum = f.readline().strip().split()
        dat[0:3, ii] = np.array(dum, dtype=float)
    r = dat[0, :]
    f.close()

    convhist = {'r': r, 'temp': temp, 'jjme': jjme, 'hhme': hhme, 'jj': jj,
                'hh': hh, 'kapt': kapt, 'kapj': kapj, 'kaph': kaph, 'fj': fj}

    return convhist


def readSpectrum(ext=0, file=0):
    import numpy as np
    import os.path
    fr_units = 0
    #
    # Read possible dust information
    #
    if os.path.isfile('dustopac.inp'):
        print("Found dustopac.inp, so assuming dust spectrum")
        fr_units = -1

    #
    # Read the total spectrum
    #
    if file == 0:
        if ext == 0:
            filename = 'spectrum.dat'
        else:
            filename = 'spectrum_'+str(ext)+'.dat'
    else:
        filename = file

    f = open(filename, 'r')

    nrfr = int(f.readline().strip())
    f.readline()
    dat = np.zeros((2, nrfr), float)
    for ii in range(0, nrfr):
        dum = f.readline().strip()
        dat[:, ii] = np.array(dum.split(), float)

    freq = dat[0, :]
    spectrum = dat[1, :]
    f.close()

    readspectrum = {'spectrum': spectrum, 'nfr': nrfr,
                    'freq': freq, 'fr_units': fr_units, 'obs': 0}

    return readspectrum
