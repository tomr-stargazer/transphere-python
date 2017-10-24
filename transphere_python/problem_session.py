"""
This is a script which does the following:

- Imports some python modules, including local ones that wrap the `transphere` and `ratran` executables
- Defines physical and code-control constants for the modeling
- Calls `transphere` 
- Plots the output of `transphere`
- Optionally

"""

import os
import sys
import math
import pdb

import numpy as np
from matplotlib import pyplot as plt
from scipy.integrate import trapz as integrate
import scipy.interpolate

import astroProcs
import natconst as nc
import transphereProcs as tP
import plotTransphere

### Parameters of physical model

lstar    = 3                 # Stellar luminosity in units of solar luminosities
tstar    = 1000.             # Stellar temperature
rstar    = lstar**0.5*(tstar/5785.0)**(-2.0)*nc.RS
rin      = 1.0 * nc.AU       # Inner radius of shell
rout     = 1.0e4 * nc.AU     # Outer radius of shell
r0       = 10.0 * nc.AU      # Reference radius
#rho0     = 9.0e-16          # Density at reference radius (gas mass). Can be used instead of Menv (need to change code further down)
menv     = 0.5               # Envelope mass in M_sun
plrho    = -1.5              # Powerlaw for rho
isrf     = 0.0               # Scaling of ISRF
dist     = 235.0             # Distance in pc
tbg      = -1                # Spectral shape of ISRF (Blackbody equivalent temperature). With -1 the spectrum is read from the file isrf.inp and scaled by ISRF.

### Parameters related to control of code

nr       = 200               # Nr of radial shells for grid
nref     = 100               # Nr of refinement points
rref     = 10. * nc.AU       # Refinement radius
nriter   = 30                # Maximum nr of iterations
convcrit = 0.00001           # Convergence criterion
ncst     = 10                # Nr of rays for star
ncex     = 30                # Nr of rays between star and Rin
ncnr     = 1                 # Nr of rays per radial grid point
itypemw  = 1                 # Type of mu weighting
idump    = 1                 # Dump convergence history
localdust= 0                 # Dust opacity local?

o = tP.readopac(nr='oh5')           ### Read-in dust opacity file
tP.writeopac(o, localdust, nr)      ### Write-out opacity as function of radius
kappa = tP.findKappa(localdust, o)  ### Find Kappa at 0.55 micron

r = tP.makeRadialGrid(nref, rin, rout, rref, nr)  # Make radial grid

# print 4.0*math.pi/(3.0+plrho)*rho0*r0**(-plrho)*(rout**(3.0+plrho)-rin**(3.0+plrho))/1.989e33
# If rho0 is given above then uncomment here.
rho0 = menv/(4.0*math.pi/(3.0+plrho)*r0**(-plrho) *
             (rout**(3.0+plrho)-rin**(3.0+plrho))/1.989e33)

rho = 1e-2*rho0 * (r/r0)**(plrho)
#tau    = integrate(rho*kappa,r)
# print 'Tau = ',tau

model = {"rstar": rstar, "mstar": 1.0, "tstar": tstar, "r": r, "rho": rho, "isrf": isrf, "tbg": tbg, 
    "freq": o['freq'], "nriter": nriter, "convcrit": convcrit, "ncst": ncst, "ncex": ncex, 
    "ncnr": ncnr, "itypemw": itypemw, "idump": idump}
tP.writeTransphereInput(model)

os.system('transphere')  # Change this to a popen call or something

s = tP.readSpectrum()
a = tP.readConvhist()

# Plot SED

import readSed

sed = readSed.readSed('sed.dat')

plt.figure(1)
plotTransphere.plotSpectrum(s, dpc=dist, jy=1, pstyle='b-')
plt.xscale('log')
plt.yscale('log')
plt.xlim((1, 3.0e3))
plt.ylim((1.0e-9, 1.0e3))
z = {"freq": o['freq'], "spectrum": math.pi *
     astroProcs.bplanck(o['freq'], tstar)*rstar**2/nc.pc**2}
plotTransphere.plotSpectrum(z, dpc=dist, jy=1, pstyle='g--')

plt.plot(sed['wave'], sed['flux'], 'r*')


# Plot temperature

plt.figure(3)
plotTransphere.plotTemperature(a, -1, pstyle='b-')
## for i in range(0,len(a['temp'][:,0])): plotTransphere.plotTemperature(a,i,pstyle='g--')

plt.show()
# sys.exit()

# Script usually stops here. Following lines can be used to create and run
# Ratran models.

import transphereRatran

rx = tP.makeRadialGrid(0, rin, rout, 0, 60.0)
rhof = scipy.interpolate.interp1d(np.log10(r), np.log10(rho))
tempf = scipy.interpolate.interp1d(np.log10(r), np.log10(a['temp'][-1, :]))

rhox_short = 10**rhof(np.log10(rx[1:]))
rhox = np.insert(rhox_short, 0, rho[0])
tempx_short = 10**tempf(np.log10(rx[1:]))
tempx = np.insert(tempx_short, 0, a['temp'][-1, 0])
abund = np.zeros(len(rx))
abund[(tempx > 90.0).nonzero()] = 1.0e-7
abund[(tempx < 90.0).nonzero()] = 1.0e-9

pdb.set_trace()

transphereRatran.ratranRun(r=rx, rho=rhox, temp=tempx, db=0.5, abund=abund, dpc=dist,
                           trans='3', pixel=0.1, molfile='hco+.dat', writeonly=0, skyonly=0, 
                           unit='jypx')

os.system('rm -rf image.sky image.conv image.mom')
os.system('fits in=ratranResult_003.fits out=image.sky op=xyin')
os.system('puthd in=image.sky/restfreq value=257.558')
os.system('restor model=image.sky fwhm=0.46 mode=convol out=image.conv')
os.system('fits in=image.conv out=conimage.fits op=xyout')
os.system('moment in=image.conv out=image.mom')
os.system('imfit in=image.mom object=gaussian | grep "Major axis"')
os.system('imfit in=image.mom object=gaussian | grep "Peak"')
