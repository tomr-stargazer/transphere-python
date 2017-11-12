def readSed(filename, sep=None):
    f = open(filename, 'r')
    inpd=f.readlines()

    sed={'wave': [], 'flux': [], 'instrument': []}

    for ii in range(0, len(inpd)):
        aa, bb, cc=inpd[ii].strip().split(sep=sep)
        sed['wave'].append(float(aa))
        sed['flux'].append(float(bb))
        sed['instrument'].append(cc)

    return sed