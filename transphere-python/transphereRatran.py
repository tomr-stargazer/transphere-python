def ratranRun(r=0.0, rho=0.0, temp=0.0, db=0.0, abund=0.0, vr=0.0, tdust=0.0, dustonly=0, file='transphere.mdl', dpc=0.0, imsize=129, pixel=0.5, trans='220.0e9', writeonly=0, skyonly=1, molfile='', outputfile="ratranResult", unit='yada'):
    import natconst as nc
    import numpy as np
    import scipy.interpolate
    import sys
    import os

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
            ii+1, r1[ii]/100.0, r2[ii]/100.0, nhint[ii], nhint[ii]*abint[ii], tkint[ii], tdint[ii], db, vr)+'\n')
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
            f.write("chan=50,0.1\n")
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
