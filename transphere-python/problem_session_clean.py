import natconst as nc
import transphereProcs as tP
#~ from scipy.integrate import trapz as integrate # needed for tau calc
import plotTransphere
import astroProcs
import math
import os
from matplotlib import pyplot as plt

### Parameters of physical model

rstar    = 4.0e0 * nc.RS     # Stellar radius
mstar    = 1.0 * nc.MS       # Stellar mass ( *** forget about this - it does not matter *** )
tstar    = 5780.             # Stellar temperature
rin      = 23.4 * nc.AU       # Inner radius of shell
rout     = 2.1e4 * nc.AU     # Outer radius of shell
r0       = 23.4 * nc.AU     # Reference radius
n0       = 1.3e9               # H2 number density at reference radius
rho0     = n0*2.3*1.667e-24  # Density at reference radius GAS
#~ rho0     = 3.0e-17           # Density at reference radius
plrho    = -1.8              # Powerlaw for rho
isrf     = 0.0               # Scaling of ISRF
tbg      = 0                # Spectral shape of ISRF (Blackbody equivalent temperature). With -1 the spectrum is read from the file isrf.inp and scaled by ISRF.
dpc      = 220.0             # Distance in pc

### Parameters related to control of code

nr       = 200              # Nr of radial shells for grid
nref     = 0                # Nr of refinement points
rref     = 0. * nc.AU       # Refinement radius

nriter   = 30                # Maximum nr of iterations
convcrit = 1E-5           # Convergence criterion
ncst     = 10                # Nr of rays for star
ncex     = 30                # Nr of rays between star and Rin
ncnr     = 1                 # Nr of rays per radial grid point
itypemw  = 1                 # Type of mu weighting
idump    = 1                 # Dump convergence history
localdust= 0                 # Dust opacity local?

o   = tP.readopac(nr='sil')       ### Read-in dust opacity file
tP.writeopac(o,localdust,nr)      ### Write-out opacity as function of radius
#~ kappa=tP.findKappa(localdust,o)   ### Find Kappa at 0.55 micron

r = tP.makeRadialGrid(nref, rin, rout, rref, nr) ### Make radial grid

rho    = 1e-2 * rho0 * (r/r0)**(plrho)  # DUST
#tau    = integrate(rho*kappa,r)  # need kappa at 0.55 micron
#print 'Tau = ',tau

model={"rstar": rstar, "mstar": mstar, "tstar": tstar, "r": r, "rho": rho, "isrf": isrf, "tbg": tbg, "freq": o['freq'], "nriter": nriter, "convcrit": convcrit, "ncst": ncst, "ncex": ncex, "ncnr": ncnr, "itypemw": itypemw, "idump": idump}

tP.writeTransphereInput(model)

os.system('../src/transphere')

s = tP.readSpectrum()
a = tP.readConvhist()

plt.figure(2)
plotTransphere.plotTemperature(a,-1,pstyle='b-')
for i in range(0,len(a['temp'][:,0])): plotTransphere.plotTemperature(a,i,pstyle='g--')


import sys; sys.exit()
## Everything here is commented out. These are examples of plots that can be done.

## plt.figure(1)
## plt.clf()
## plotTransphere.plotSpectrum(s,dpc=dpc,jy=1,pstyle='b-')
## plt.xscale('log')
## plt.yscale('log')
## plt.xlim((1,3.0e3))
## plt.ylim((1.0e-3,6.0e2))
## z={"freq": o['freq'], "spectrum": math.pi*astroProcs.bplanck(o['freq'],tstar)*rstar**2/nc.pc**2}
## plotTransphere.plotSpectrum(z,dpc=dpc,jy=1,pstyle='g--')

## from readSed import *
## sed1=readSed('bhr71_sed.dat')
## sed2=readSed('bhr71_sed_herschel.dat')

## plt.plot(sed1['wave'],sed1['flux'],'r.')
## plt.plot(sed2['wave'],sed2['flux'],'g.')
## ## plt.figure(1)
## ## plotTransphere.plotSpectrum(s,dpc=dpc,pstyle='b-')
## ## plt.xscale('log')
## ## plt.yscale('log')
## ## plt.xlim((1,3.0e3))
## ## plt.ylim((1.0e-19,1.0e-7))
## ## z={"freq": o['freq'], "spectrum": math.pi*astroProcs.bplanck(o['freq'],tstar)*rstar**2/nc.pc**2}
## ## plotTransphere.plotSpectrum(z,dpc=dpc,pstyle='g--')

## plt.figure(2)
## plotTransphere.plotTemperature(a,-1,pstyle='b-')
## for i in range(0,len(a['temp'][:,0])): plotTransphere.plotTemperature(a,i,pstyle='g--')

## plt.show()

## sys.exit('Halting...')

from scipy import ones, where

i = where(a['temp'][-1]>100)[0]
x_in = 1e-6/560
x_out = 1e-10/560
mol_abundance = ones(len(rho))*x_out
mol_abundance[i] = x_in



import transphereRatran; reload(transphereRatran)

pixelsize=0.2
imsize=251

transphereRatran.ratranRun(r=r, rho=rho, molfile='ph2o-h2.dat', db=1.0, abund=mol_abundance, temp=a['temp'][-1], trans='7', dpc=dpc, imsize=imsize, pixel=pixelsize, skyonly=0, snr=10, fixset=1e-3, ncell=10, writeonly=1)

## The following lines convolves the resulting image with a beam in
## Miriad. If the output is in K, you need to introduce a scale factor
## in the convolution (check Christian or Floris van der Tak's
## homepages)

os.system('rm -rf image.sky image.conv')
os.system('fits in=ratranResult_3.530E+11.fits out=image.sky op=xyin')
os.system('restor model=image.sky fwhm=1.0 mode=convol out=image.conv')
os.system('fits in=image.conv out=conimage.fits op=xyout')

## More plotting. Probably not necessary.

## import daoproc

## class options: pass
## options.filename='conimage.fits'

## map=daoproc.readMap(options)
## plt.figure(3)
## daoproc.plotMap(map,options)

## plt.figure(4)
## plt.plot(map.xaxis,map.data[math.floor(imsize/2.0),:])
## plt.yscale('log')
## plt.xlim(0,40)
## plt.ylim(0.00001,0.01)
## plt.xlabel('b ["]')
## plt.ylabel('I [Jy beam$^{-1}]$')
## plt.grid(True,which="both")
