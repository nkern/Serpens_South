"""
This code takes the gaussian fitted data and derives the Spectral Index of the sources with fits

May, 2014
"""

### Import Modules ###
import numpy as np
import matplotlib.pyplot as mp
import scipy.stats as stats

### Flags, Constants ###
lower_freq	= 4750.00		# MHz
upper_freq	= 7250.00		# MHz
lower_passband	= [4238.00,5262.00]	# MHz
upper_passband	= [6738.00,7762.00]	# MHz
Source_Num	= 21
Image_Num	= 2
c		= 2.9979e8		# m/s


### Program ###

# Load Data, arrays with 0.0000 are nulls
data = np.loadtxt('Gauss_Fitting/Final_Fit_Results.txt',usecols=(1,2,3,4,5,6,7,8,9,10))
sources = np.loadtxt('Gauss_Fitting/Final_Fit_Results.txt',usecols=(0,),dtype=str)

lower_data = data[0:21]
lower_sources = sources[0:21]
upper_data = data[21:]
upper_sources = sources[21:]

# data arrays have 2nd axis metadata of:
# flux,flux_e,ra,dec,ra_e,dec_e,bmaj,bmin,pa,rms

# Find sources with fit in both upper and lower
fitted = np.array([0,1,2,3,4,5,6,7,10,11,12,14,16,18,19,20])	# Indicies

# Load Gutermuth Cross Match data
RA_VLA,Dec_VLA,jmag,hmag,kmag,i1mag,i2mag,i3mag,i4mag,m24mag,p70mag = np.loadtxt('vla_serps_merger_gutermuth_try1.csv',delimiter=',',unpack=True)
# -100 are nulls
# jhk mags are 2mass near infrared 
# i1,i2,i3,i4 are IRAC on Spitzer
# m24 is MIPS 24 micron on Spitzer 
# p70 is PACS 70 micron on Herschel

# Write ds9 region file for RA and DEC coordinates
output = open('Gutermuth_CrossMatch.reg','w')
output.write('# Region file format: DS9 version 4.1\n')
output.write('global color=magenta dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
output.write('fk5\n')
for i in range(len(RA_VLA)):
	output.write('point('+str(RA_VLA[i])+','+str(Dec_VLA[i])+') # point=x text={'+str(i)+'}\n')

output.close()
raise NameError


## Convert magnitude to flux
# Flux Zero Points
jmag_lam	= 1.235e-6	# meter
jmag_lam_err	= 0.006e-6	# meter
jmag_zero	= 1594		# Jy
jmag_zero_err	= 28		# Jy

hmag_lam	= 1.662e-6	
hmag_lam_err	= 0.009e-6
hmag_zero	= 1024
hmag_zero_err	= 20

kmag_lam	= 2.159e-6
kmag_lam_err	= 0.011e-6
kmag_zero	= 666.7
kmag_zero_err	= 12.6

i1mag_lam	= 3.6e-6
i1mag_zero	= 280.9
i1mag_zero_err	= 4.1

i2mag_lam	= 4.5e-6
i2mag_zero	= 179.7
i2mag_zero_err	= 2.6

i3mag_lam	= 5.8e-6
i3mag_zero	= 115.0
i3mag_zero_err	= 1.7

i4mag_lam	= 8.0e-6
i4mag_zero	= 64.9
i4mag_zero_err	= 0.9

m24mag_lam	= 23.68e-6
m24mag_zero	= 7.17
m24mag_zero_err	= 0.11

p70mag_lam	= 71.42e-6
p70mag_zero	= 0.778
p70mag_zero_err	= 0.012

# Convert Mag to Flux
def mag2flux(mag,zero):
	# Converts magnitude to flux, given a zero point flux
	return 10**(-mag/2.5)*zero

jflux		= mag2flux(jmag,jmag_zero)

hflux		= mag2flux(hmag,hmag_zero)

kflux		= mag2flux(kmag,kmag_zero)

i1flux		= mag2flux(i1mag,i1mag_zero)

i2flux		= mag2flux(i2mag,i2mag_zero)

i3flux		= mag2flux(i3mag,i3mag_zero)

i4flux		= mag2flux(i4mag,i4mag_zero)

m24flux		= mag2flux(m24mag,m24mag_zero)

p70flux		= mag2flux(p70mag,p70mag_zero)


guter_fluxes = np.vstack([jflux,hflux,kflux,i1flux,i2flux,i3flux,i4flux,m24flux,p70flux]).T
# Outrageously high fluxes are nulls (above 1e10 Jansky)

# Plot SEDs

raise NameError
## Calculate Spectral Indicies ##
def radio_spix(flux1,flux2,freq1,freq2):
	return np.log(flux1/flux2)/np.log(freq2/freq1)

def infr_spix(flux1,flux2,lam1,lam2):
	# lam1 should be larger than lam2
	freq1 = c/lam1
	freq2 = c/lam2
	return (np.log(lam1*flux1)-np.log(lam2*flux2))/(np.log(lam1)-np.log(lam2))

# Infrared Spectral Index






alpha_ideal = []
for i in range(Source_Num):
	if (i in fitted) == True:
		alpha_ideal.append(SI_calc2(lower_data[i][0],upper_data[i][0],lower_freq,upper_freq))
	else:
		alpha_ideal.append(None)

alpha_ideal = np.array(alpha_ideal)











