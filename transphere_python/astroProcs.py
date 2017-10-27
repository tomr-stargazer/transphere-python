import transphere_python.natconst as nc

def bplanck(nu, T):
    import numpy as np
    import math
    n = len(nu)
    bpl = np.zeros((n), float)
    if T != 0.0:
        x = nc.hh*nu/(nc.kk*T)
        bpl = (2.e0*nc.hh*nu**3/nc.cc**2)/(math.e**x-1.e0)

    return bpl
