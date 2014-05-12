#### Serpens South 2D Gaussian Fit ####
## Nicholas S. Kern
## nkern@umich.edu

### Import Modules
import numpy as np
import matplotlib.pyplot as mp
import scipy.optimize
import astropy.io.fits as fits
from mpl_toolkits.mplot3d import Axes3D
from GaussFit import gaussian, moments, fitgaussian

### Flags



### Program

## Load Data
# Load FITS
HDU1 = fits.open('SERPS_2.upper.mask.deeper2.briggs.fits')
SS_4cm = HDU1[0].data[0][0]

## Load Approximate Source Locations
# Load SS_Targets_pxl.csv
Y,X,ID = np.loadtxt('SS_Targets_pxl.csv',delimiter=',',unpack=True)
ID = np.array(ID,int)
sort = np.argsort(ID)
X,Y,ID = X[sort],Y[sort],ID[sort]

## Cut regions
X_sub = X-5
X_add = X+5
Y_sub = Y-5
Y_add = Y+5
SS_4cm_sources = []
for i in range(len(ID)):
	SS_4cm_sources.append(SS_4cm[X_sub[i]:X_add[i],Y_sub[i]:Y_add[i]])

SS_4cm_sources = np.array(SS_4cm_sources)

## Find Maximum Value per Source
max_val = np.array(map(np.max,SS_4cm_sources))

## Fit 2D Gaussian
PARAMS = []
for i in range(len(ID)):
	PARAMS.append(fitgaussian(SS_4cm_sources[i]))

PARAMS=np.array(PARAMS)

## Integrate FWHM
# 2D Gaussian Functionalform
# f = A*np.exp(-( (x-mu_x)**2/(2*sig_x**2) + (y-mu_y)**2/(2*sig_y**2) ))
def f(x,y,A,mu_x,mu_y,sig_x,sig_y):
	return A*np.exp(-( (x-mu_x)**2/(2*sig_x**2) + (y-mu_y)**2/(2*sig_y**2) ))




## Plot
i = 5
mp.matshow(SS_4cm_sources[5])
X,Y = np.meshgrid(range(10),range(10))
CS = mp.contour(f(X,Y,*PARAMS[5]))
mp.show()


