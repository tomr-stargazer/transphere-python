def plotSpectrum(a, dpc=0, jy=0, pstyle='', xlog=1, ylog=1):
    from matplotlib import rc
    rc('text', usetex=True)
    rc('font', family='serif')
    import natconst as nc
    import sys
    from matplotlib import pyplot as plt
    xcoord = 1.0e4*nc.cc / a['freq']
    
    if dpc == 0: sys.exit('Error: distance needs to be set when plotting flux')

    distfact = 1.e0/ (dpc**2)

    if jy != 0:
        lumfact=1e+23
    else:
        lumfact=a['freq']

    plt.plot(xcoord,distfact*lumfact*a['spectrum'],pstyle)
    plt.xlabel(r'$\lambda\, [\mu \textrm{m}]$',fontsize=14)

    if jy != 0:
        plt.ylabel(r'$F_\nu$\, [Jy]',fontsize=14)
    else:
        plt.ylabel(r'$\nu F_\nu \, [\textrm{erg cm}^{-2}\, \textrm{s}^{-1}]$',fontsize=14)
        
    if xlog == 1: plt.xscale('log')
    if ylog == 1: plt.yscale('log')



def plotTemperature(a, model=-1, pstyle='', xlog=1, ylog=1):
    from matplotlib import rc
    import numpy as np
    rc('text', usetex=True)
    rc('font', family='serif')
    import natconst as nc
    from matplotlib import pyplot as plt

    plt.plot(a['r']/nc.AU,a['temp'][model,:],pstyle)
    plt.xlabel(r'$R\, [\textrm{AU}]$',fontsize=14)
    plt.ylabel(r'$T\, [\textrm{K}]$',fontsize=14)
    if xlog == 1: plt.xscale('log')
    if ylog == 1: plt.yscale('log')

    xmin=np.floor(np.log10(min(a['r']/nc.AU)))
    xmax=np.ceil(np.log10(1.2*max(a['r']/nc.AU))*10.0)
    plt.xlim((10**xmin,10**(xmax/10.0)))
    plt.ylim((5,150))
