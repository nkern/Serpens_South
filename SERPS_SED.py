"""
This code takes the gaussian fitted data and derives the Spectral Index of the sources with fits
Also does other stuffs...

December, 2014
nkern@umich.edu
"""

### Import Modules ###
print '...importing modules'
import numpy as np
import matplotlib.pyplot as mp
import scipy.stats as stats
import numpy.ma as ma
from mpl_toolkits.axes_grid1 import AxesGrid

print '...done importing'

### Flags, Constants ###
lower_freq	= 4750.00		# MHz
upper_freq	= 7250.00		# MHz
lower_passband	= [4238.00,5262.00]	# MHz
upper_passband	= [6738.00,7762.00]	# MHz
Source_Num	= 18
Image_Num	= 2
c		= 2.9979e8		# m/s

### Program ###


###############################################################################
########################### Load Data #########################################
###############################################################################

print '...loading Gutermuth Cross Match Data'
## Load Gutermuth Cross Match data
RA_VLA,Dec_VLA,jmag,hmag,kmag,i1mag,i2mag,i3mag,i4mag,m24mag,p70mag = np.loadtxt('vla_serps_merger_gutermuth_try1.csv',delimiter=',',unpack=True)
# -100 are nulls
# jhk mags are 2mass near infrared 
# i1,i2,i3,i4 are IRAC on Spitzer
# m24 is MIPS 24 micron on Spitzer 
# p70 is PACS 70 micron on Herschel

# Write ds9 region file for RA and DEC coordinates
#output = open('Gutermuth_CrossMatch.reg','w')
#output.write('# Region file format: DS9 version 4.1\n')
#output.write('global color=magenta dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
#output.write('fk5\n')
#for i in range(len(RA_VLA)):
#	output.write('point('+str(RA_VLA[i])+','+str(Dec_VLA[i])+') # point=x text={'+str(i+1)+'}\n')

#output.close()


## Convert magnitude to flux ##
# Flux Zero Points
jmag_lam	= 1.235e-6	# meter
jmag_zero	= 1594		# Jy
jmag_zero_err	= 28		# Jy

hmag_lam	= 1.662e-6	
hmag_zero	= 1024
hmag_zero_err	= 20

kmag_lam	= 2.159e-6
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

# Stack Fluxes into single array
guter_fluxes = np.vstack([jflux,hflux,kflux,i1flux,i2flux,i3flux,i4flux,m24flux,p70flux]).T # large numbers are nulls
wavelens = np.array([jmag_lam,hmag_lam,kmag_lam,i1mag_lam,i2mag_lam,i3mag_lam,i4mag_lam,m24mag_lam,p70mag_lam])*1e6 # microns

# Make nulls as zero
guter_fluxes[np.where(guter_fluxes>1e5)] = 0.0

# Magnitude Errors, for now we assume magnitude error is 0.15, which corresponds to +/- 15% in flux
mag_rel_errs = np.array([0.15]*9)

# Because flux is directly prop. to zero_offset, error in flux offset is just flux * relative error of zero_offset
zeroflux_rel_errs = np.array([jmag_zero_err,hmag_zero_err,kmag_zero_err,i1mag_zero_err,i2mag_zero_err,i3mag_zero_err,i4mag_zero_err,m24mag_zero_err,p70mag_zero_err])/np.array([jmag_zero,hmag_zero,kmag_zero,i1mag_zero,i2mag_zero,i3mag_zero,i4mag_zero,m24mag_zero,p70mag_zero])

# Total Error On Flux (Relative) add in quadrature
guter_rel_errs = np.sqrt( mag_rel_errs**2 + zeroflux_rel_errs**2)


########################################################################################
################################### Radio Spectral Indices #############################
########################################################################################

print '...Radio Spix Calculations'
lower_lam = 6.311	# cm
upper_lam = 4.135	# cm

# Load Data, arrays with 0.0000 are nulls
data = np.loadtxt('Gauss_Fitting/Final_Fit_Results_try2.txt',usecols=(1,2,3,4,5,6,7,8,9,10))
sources = np.loadtxt('Gauss_Fitting/Final_Fit_Results_try2.txt',usecols=(0,),dtype=str)

# Fix data's ra to positive degree value
data.T[2] = data.T[2]+360

# Separate Data
source_num = 18
lower_data = data[0:source_num]
lower_sources = sources[0:source_num]
upper_data = data[source_num:]
upper_sources = sources[source_num:]
# data arrays have 2nd axis metadata of:
# flux,flux_e,ra,dec,ra_e,dec_e,bmaj,bmin,pa,rms

## Radio Spectral Index, Shirley et al. 2007
def radio_spix(flux1,flux2,lam1,lam2):
	''' lam2 is larger than lam1'''
	return np.log(flux1/flux2)/np.log(lam2/lam1)

def radio_spix_err(flux1,flux1_err,flux2,flux2_err):
	# Add 5% errors to flux errors for calibration error
	return np.sqrt(  (flux1_err/flux1)**2 + (flux2_err/flux2)**2 )


# Load cm Peak Flux Data
lower_peak_flux = np.loadtxt('Lower_Peak_Fluxes_try2.csv',delimiter=',',unpack=True)
upper_peak_flux = np.loadtxt('Upper_Peak_Fluxes_try2.csv',delimiter=',',unpack=True)
cm_int_spix = []
cm_int_spix_err = []
cm_peak_spix = []
cm_peak_spix_err = []
upper_rms = 8.5e-6		# Jy
lower_rms = 11.1e-6		# Jy
for i in range(Source_Num):
		cm_int_spix.append(	radio_spix(upper_data[i][0],lower_data[i][0],upper_lam,lower_lam) )
		upper_data[i][1] = upper_data[i][0] * np.sqrt( (upper_data[i][1]/upper_data[i][0])**2 + 0.05**2 )
		lower_data[i][1] = upper_data[i][0] * np.sqrt( (lower_data[i][1]/lower_data[i][0])**2 + 0.05**2 )
		cm_int_spix_err.append(	radio_spix_err(upper_data[i][0],upper_data[i][1],lower_data[i][0],lower_data[i][1]) )
		cm_peak_spix.append(	radio_spix(upper_peak_flux[1][i],lower_peak_flux[1][i],upper_lam,lower_lam) )
		upper_err = upper_peak_flux[1][i] * np.sqrt( (upper_rms/upper_peak_flux[1][i])**2 + 0.05**2 )
		lower_err = lower_peak_flux[1][i] * np.sqrt( (lower_rms/lower_peak_flux[1][i])**2 + 0.05**2 )
		cm_peak_spix_err.append(radio_spix_err(upper_peak_flux[1][i],upper_err,lower_peak_flux[1][i],lower_err) )

cm_int_spix = np.array(cm_int_spix)
cm_int_spix_err = np.array(cm_int_spix_err)
cm_peak_spix = np.array(cm_peak_spix)
cm_peak_spix_err = np.array(cm_peak_spix_err)
# Note that error in spix is just error in top logarithm, which is: delta log(x) = delta x / x

# Which spectral index to trust based on observational properties
cm_spix_trust = ['int','int','int','int','int','peak','int','peak','peak','peak','peak','peak','peak','peak','int','peak','int','peak']

cm_spix = []
cm_spix_err = []
for i in range(source_num):
	if cm_spix_trust[i] == 'int':
		cm_spix.append(cm_int_spix[i])
		cm_spix_err.append(cm_int_spix_err[i])
	elif cm_spix_trust[i] == 'peak':
		cm_spix.append(cm_peak_spix[i])
		cm_spix_err.append(cm_peak_spix_err[i])

cm_spix = np.array(cm_spix)
cm_spix_err = np.array(cm_spix_err)	

cm_lower_snr = lower_peak_flux[1] / lower_rms
cm_upper_snr = upper_peak_flux[1] / upper_rms


##########################################################################################
################################ Infrared Spectral Index #################################
##########################################################################################

print '...Infrared Spix Calculations'
# Use 8 micron to 24 micron regime (i4flux to m24flux; guter_fluxes[6] to guter_fluxes[7])
# 7 sources with 8 micron and 24 micron infrared data

# match gutermuth coordinates with our coordinates with 3 arcsec tolerance
epsilon = 3./3600.

ir_source_id = []
upper_dec = upper_data.T[3]
upper_ra = upper_data.T[2]

for i in range(source_num):
	try:
		ir_source_id.append( np.where((Dec_VLA<upper_dec[i]+epsilon) & (Dec_VLA>upper_dec[i]-epsilon) & (RA_VLA<upper_ra[i]+epsilon) & (RA_VLA>upper_ra[i]-epsilon))[0][0] )
		#print np.sqrt( (Dec_VLA[ir_source_id[-1]]-upper_dec[i])**2 + (RA_VLA[ir_source_id[-1]]-upper_ra[i])**2) 
		if guter_fluxes[ir_source_id[-1]].sum() == 0.0: ir_source_id[-1] = -99
	except IndexError:
		ir_source_id.append(-99)
	
ir_source_id = np.array(ir_source_id)

guter_fluxes_temp = np.copy(guter_fluxes)
guter_fluxes = []
for i in range(source_num):
	if ir_source_id[i] == -99:
		guter_fluxes.append(np.array([np.nan]*9))
	else:
		guter_fluxes.append(guter_fluxes_temp[ir_source_id[i]])

guter_fluxes = np.array(guter_fluxes)

def infr_spix(flux1,flux2,lam1,lam2,flux_rel_err1,flux_rel_err2):
	# lam1 should be larger than lam2
	freq1 = c/lam1
	freq2 = c/lam2
	spix = (np.log(freq1*flux1/freq2/flux2))/(np.log(lam1/lam2))
	err = np.sqrt( ((freq1*flux_rel_err1*flux1)/(freq1*flux1))**2 + ((freq2*flux_rel_err2*flux2)/(freq2*flux2))**2) / np.log(lam1/lam2)
	return spix,err

ir_lower_lam = 8e-6	# m
ir_upper_lam = 23.68e-6	# m

ir_spix = []
ir_spix_err = []

# Sources with 8 and 24 micron data
ir_fitted = np.array([6,10,12,15,16])

# infrared pointings
for i in range(source_num):
	if (i in ir_fitted) == True:
		ir_spix.append( infr_spix(guter_fluxes[i][7],guter_fluxes[i][6],ir_upper_lam,ir_lower_lam,guter_rel_errs[7],guter_rel_errs[6])[0] )
		ir_spix_err.append( infr_spix(guter_fluxes[i][7],guter_fluxes[i][6],ir_upper_lam,ir_lower_lam,guter_rel_errs[7],guter_rel_errs[6])[1] )
	else:
		ir_spix.append(np.nan)
		ir_spix_err.append(np.nan)	

ir_spix = np.array(ir_spix)
ir_spix_err = np.array(ir_spix_err)



##########################################################################################
################################ Combine IR and CM Data ##################################
##########################################################################################

def cross_match(array1,array2):
	'''Takes array1 and matches values from array2, w/ -99 being np.nan'''
	array1_temp = []
	for i in range(len(array2)):
		if array2[i] == -99: array1_temp.append(np.nan)
		else: array1_temp.append(array1[array2[i]])
	return np.array(array1_temp)

def sigfigs(array,cut_str):
	''' iterates a float manipulation command
	try using "'%.3f'" as cut_str '''
	return np.array(map(lambda x:cut_str % x, array))
	
def replace_nan(array,sub,round_num=2,flt=True):
	''' replaces numpy array's np.nan with sub '''
	length = len(array)
	new_array = []
	try: nans = np.isnan(array)
	except: nans = array == 'nan'
	for i in range(length):
		if nans[i] == True:
			new_array.append(sub)
		else:
			if flt == True:
				new_array.append(('%.'+str(round_num)+'f')%array[i])
			else:
				new_array.append(array[i])
	return np.array(new_array)

## Create data dictionary for arrays with ordering identical to SS_Final_Sources
data = {}

# Create Source Names Array
source_names = []
for i in range(source_num):
	source_names.append('VLA_'+str(i+1))
source_names = np.array(source_names)

# Order by descending Dec, and match cm sources with infra sources
# 18 cm sources
# 23 infrared sources

# Load Gutermuth Spitzer Protostellar Data
# index, RA, DEC, Class
spitz_id,spitz_ra,spitz_dec = np.loadtxt('gutermuth2008_serps_ysos1.txt',usecols=(0,1,2),unpack=True)
spitz_class = np.loadtxt('gutermuth2008_serps_ysos1.txt',usecols=(3,),unpack=True,dtype='str')

# Match Spitzer Data with 3 arcsec tolerance
spitz_source_id = []
for i in range(source_num):
	try:
		spitz_source_id.append( np.where((spitz_dec<upper_dec[i]+epsilon) & (spitz_dec>upper_dec[i]-epsilon) & (spitz_ra<upper_ra[i]+epsilon) & (spitz_ra>upper_ra[i]-epsilon))[0][0] )
	except IndexError:
		spitz_source_id.append(-99)

spitz_source_id = np.array(spitz_source_id)

# Match w/ SS_Final_Source Ordering
ir_spitz_class = cross_match(spitz_class,spitz_source_id)


### Dunham08 70micron vs. Linternal ###
Lint_dunham = (guter_fluxes.T[-1] * c / 70e-6 * 1e-23) ** (0.94) * 3.3e8
Lint_dunham_err = 3.3e8 * guter_fluxes.T[-1]*guter_rel_errs[-1]*c/70e-6 * 1e-23		# Lsun

### Kryukova12 Lmir vs. Lbol ###
dist1 = 429	# pc
dist2 = 260	# pc

Lmir1 = (19.79*guter_fluxes.T[0] + 16.96*guter_fluxes.T[1] + 10.49*guter_fluxes.T[2] + 5.50*guter_fluxes.T[3] + 4.68*guter_fluxes.T[4] + 4.01*guter_fluxes.T[5] + 4.31*guter_fluxes.T[6] + 0.81*guter_fluxes.T[7]) * 1e-6 * dist1**2

Lmir1_err = (19.79*guter_fluxes.T[0]*guter_rel_errs[0] + 16.96*guter_fluxes.T[1]*guter_rel_errs[1] + 10.49*guter_fluxes.T[2]*guter_rel_errs[2] + 5.50*guter_fluxes.T[3]**guter_rel_errs[3] + 4.68*guter_fluxes.T[4]*guter_rel_errs[4] + 4.01*guter_fluxes.T[5]*guter_rel_errs[5] + 4.31*guter_fluxes.T[6]*guter_rel_errs[6] + 0.81*guter_fluxes.T[7]*guter_rel_errs[7]) * 1e-6 * dist1**2

Lmir2 = (19.79*guter_fluxes.T[0] + 16.96*guter_fluxes.T[1] + 10.49*guter_fluxes.T[2] + 5.50*guter_fluxes.T[3] + 4.68*guter_fluxes.T[4] + 4.01*guter_fluxes.T[5] + 4.31*guter_fluxes.T[6] + 0.81*guter_fluxes.T[7]) * 1e-6 * dist2**2

Lmir2_err = (19.79*guter_fluxes.T[0]*guter_rel_errs[0] + 16.96*guter_fluxes.T[1]*guter_rel_errs[1] + 10.49*guter_fluxes.T[2]*guter_rel_errs[2] + 5.50*guter_fluxes.T[3]**guter_rel_errs[3] + 4.68*guter_fluxes.T[4]*guter_rel_errs[4] + 4.01*guter_fluxes.T[5]*guter_rel_errs[5] + 4.31*guter_fluxes.T[6]*guter_rel_errs[6] + 0.81*guter_fluxes.T[7]*guter_rel_errs[7]) * 1e-6 * dist2**2

## Flat and Rising IR Sources
flat = np.array([6,16])
rising = np.array([12])#10,12,15])
Lbol_sources = np.concatenate([flat,rising])

Lbol_kryukova1 = np.concatenate([ Lmir1[flat] * (-0.466 * np.log10(0.3) + 0.337)**2, Lmir1[rising] / (-0.466 * np.log10(ir_spix[rising]) + 0.337)**2 ])
Lbol_kryukova2 = np.concatenate([ Lmir2[flat] * (-0.466 * np.log10(0.3) + 0.337)**2, Lmir2[rising] / (-0.466 * np.log10(ir_spix[rising]) + 0.337)**2 ])

Delta = (-0.466 * np.log10(cm_spix[Lbol_sources]) + 0.337)
Delta_err = np.sqrt( (0.014 * cm_spix_err[Lbol_sources]/cm_spix[Lbol_sources])**2 + 0.033**2)
Beta = Delta**2
Beta_err = Beta * np.sqrt((Delta_err/Delta)**2 * 2)

Lbol_kryukova1_err = Lbol_kryukova1 * np.sqrt( (Lmir1_err[Lbol_sources]/Lmir1[Lbol_sources])**2 + (Beta_err/Beta)**2 )
Lbol_kryukova2_err = Lbol_kryukova2 * np.sqrt( (Lmir2_err[Lbol_sources]/Lmir2[Lbol_sources])**2 + (Beta_err/Beta)**2 )

### Dunham13 Lmir vs. Lbol ###
Lbol_dunham1 = np.concatenate( [Lmir1[flat] / (-0.298 * np.log10(0.3) + 0.270)**2, Lmir1[rising] / (-0.298 * np.log10(ir_spix[rising]) + 0.270)**2] )
Lbol_dunham2 = np.concatenate( [Lmir2[flat] / (-0.298 * np.log10(0.3) + 0.270)**2, Lmir2[rising] / (-0.298 * np.log10(ir_spix[rising]) + 0.270)**2] )

# Classification
classes = ['Extragal.?','Extragal.','Extragal.','Extragal.','Extragal.?','Extragal.','Class I','Extragal.?','Class 0?','Unknown','Class 0','Class 0','Class 0','Class 0','Extragal.','Class II','Class I','Extragal.?']

# Convert RA and DEC to J2000
ra = upper_ra
ra_frac = ra/360. * 24
ra_h = np.array(np.floor(ra_frac),int)
ra_m = np.array(np.floor((ra_frac - ra_h)*60),int)
ra_s = np.around(((ra_frac - ra_h)*60 - ra_m)*60,2)
ra_s = np.array(map(lambda x:"%05.2f" % x, ra_s))
ra_sexig = np.array([a+' '+b+' '+c for a,b,c in zip(np.array(ra_h,str),np.array(ra_m,str),np.array(ra_s,str))])

dec = upper_dec
dec_deg = np.ceil(dec)
dec_m = np.array(np.floor((dec_deg-dec)*60),int)
dec_s = np.around( ((dec_deg-dec)*60 - dec_m)*60,1)
dec_deg = np.array(map(lambda x: "%03d" % (x), dec_deg))
dec_m = np.array(map(lambda x:"%02d" % (x), dec_m))
dec_s = np.array(map(lambda x:"%04.1f" % (x), dec_s))
dec_sexig = np.array([a+' '+b+' '+c for a,b,c in zip(np.array(dec_deg,str),np.array(dec_m,str),np.array(dec_s,str))])

## Put Data into "Data" Dictionary
# Order arrays into SS_Final_Sources ordering
data['source_names'] = np.array(map(lambda x: x[0:3]+'\\'+x[3:],source_names))
data['ra'] = np.array(ra,str)
data['dec'] = np.array(dec,str)
data['ra_sexig'] = ra_sexig
data['dec_sexig'] = dec_sexig
data['flux_upper'] = np.array(np.array(upper_data.T[0]*1e6,int),str)	# Put into  microJy unit and make integer
data['flux_upper_err'] = np.array(np.array(upper_data.T[1]*1e6,int),str)
data['flux_upper_peak'] = np.array(np.array(upper_peak_flux[1]*1e6,int),str)
data['flux_upper_peak_err'] = np.str(np.around(upper_rms*1e6,1))
data['flux_lower'] = np.array(np.array(lower_data.T[0]*1e6,int),str)	# Put into  microJy unit and make integer
data['flux_lower_err'] = np.array(np.array(lower_data.T[1]*1e6,int),str)
data['flux_lower_peak'] = np.array(np.array(lower_peak_flux[1]*1e6,int),str)
data['flux_lower_peak_err'] = np.str(np.around(lower_rms*1e6,1))
data['cm_int_spix'] = np.array(np.around(cm_int_spix,2),str)
data['cm_int_spix_err'] = np.array(np.around(cm_int_spix_err,2),str)
data['cm_peak_spix'] = np.array(np.around(cm_peak_spix,2),str)
data['cm_peak_spix_err'] = np.array(np.around(cm_peak_spix_err,2),str)
data['ir_spix'] = replace_nan(ir_spix,'--',1)
data['ir_spix_err'] = replace_nan(ir_spix_err,'--',1)
data['ir_spitz_class'] = replace_nan(ir_spitz_class,'--',flt=False)
data['classes'] = classes
data['cm_spix_trust'] = np.array(cm_spix_trust)
data['cm_spix'] = np.array(cm_spix)
data['cm_spix_err'] = np.array(cm_spix_err)

data['Lbol_kryukova1'] = Lbol_kryukova1
data['Lbol_kryukova2'] = Lbol_kryukova2
data['Lbol_kryukova1_err'] = Lbol_kryukova1_err
data['Lbol_kryukova2_err'] = Lbol_kryukova2_err
data['Lbol_dunham1'] = Lbol_dunham1
data['Lbol_dunham2'] = Lbol_dunham2
data['Lbol_sources'] = np.array(Lbol_sources)

##################################################
########## PRINT STRINGS FOR FLUX TABLE ##########
##################################################

print_strings = False
if print_strings == True:

	tc = ['source_names','ra_sexig','dec_sexig','flux_lower','flux_lower_err',
	'flux_lower_peak','flux_upper','flux_upper_err','flux_upper_peak','cm_int_spix',
	'cm_int_spix_err','cm_peak_spix','cm_peak_spix_err','classes'
	]

	length = len(data['source_names'])

	for i in range(length):
		print ""
		print data[tc[0]][i]+"\t&\t"+data[tc[1]][i]+"\t&\t"+data[tc[2]][i]+"\t&\t"+data[tc[3]][i]+" $\pm$ "+data[tc[4]][i]+"\t&\t"+data[tc[5]][i]+"\t&\t"+data[tc[6]][i]+" $\pm$ "+data[tc[7]][i]+"\t&\t"+data[tc[8]][i]+"\t&\t"+data[tc[9]][i]+" $\pm$ "+data[tc[10]][i]+"\t&\t"+data[tc[11]][i]+" $\pm$ "+data[tc[12]][i]+"\t&\t"+data[tc[13]][i]+"\\\\[1ex]"



#############################################
######### W R I T E   D A T A ###############
#############################################

write_pkl = True
if write_pkl == True:
	import pickle as pkl
	file = open('Radio_Data.pkl','wb')
	output = pkl.Pickler(file)
	output.dump(data)
	file.close()

write_csv = False
if write_csv == True:		# Write Radio Properties Data Files 
	file = open('VLA_Sources_Radio_Properties.csv','w')
	file.write('# EVLA C band C configuration Observation of Serpens South Central Filament in July 2013.\n')
 	file.write('# This file contains radio properties for 18 VLA radio sources in SS.\n')
	file.write('#\n#\n# Nicholas Kern, nkern@umich.edu\n#\n# DESCRIPTION OF DATA:\n')
	file.write('# 1  Sourcename, string\n')
	file.write('# 2  Right Ascension, degrees-float\n')
	file.write('# 3  Declination, degrees-float\n')
	file.write('# 4  6.3cm GaussFit bmaj, arcsec\n')
	file.write('# 5  6.3cm GaussFit bmin, arcsec\n')
	file.write('# 6  6.3cm GaussFit PA, degrees\n')
	file.write('# 4  4.1cm GaussFit bmaj, arcsec\n')
	file.write('# 5  4.1cm GaussFit bmin, arcsec\n')
	file.write('# 6  4.1cm GaussFit PA, degrees\n')
	file.write('# 4  Integrated 6.3 cm Flux, microJy\n')
	file.write('# 5  Int. 6.3 cm Flux Error, microJy\n')
	file.write('# 6  Peak 6.3 cm Flux, microJy/beam\n')
	file.write('# 7  Integrated 4.1 cm Flux, microJy\n')
	file.write('# 8  Int. 4.1 cm Flux Error, microJy\n')
	file.write('# 9  Peak 4.1 cm Flux, microJy/beam\n')
	file.write('# 10 Integrated Flux Spix, unitless\n')
	file.write('# 11 Int. Flux Spix Error, unitless\n')
	file.write('# 12 Peak Flux Spix, unitless\n')
	file.write('# 13 Peak Flux Spix, Error, unitless\n')
	file.write('# 14 Image RMS, microJy/beam\n')
	file.write('# 15 SNR Ratio Peak Flux to RMS, unitless\n') 







###################################################################################
################################ P L O T T I N G ##################################
###################################################################################

## Plot SEDs ##
plot = False
oneplot = True
savefig = False
if plot == True:	
	print '...Making Infrared SED Plots'
	# Allow Tex Rendering, its slow though...
	#mp.rc('text', usetex=True)
	#mp.rc('font',**{'family':'sans-serif'})
	infra_seds = np.array([6,9,10,11,12,15,16,17])
	c = 2.99e8	# m/s

	if oneplot == True:
		fig = mp.figure(figsize=(13,9))
		grid = AxesGrid(fig, 111,share_all=False,
			nrows_ncols=(2,4),axes_pad=0,label_mode='L',aspect=False)
		#grid[-2].axis('off')
		#grid[-1].axis('off')
		j = 0
		for i in range(source_num):
			if i in infra_seds:
				# Add Axes Labels
				if j == 4:
					grid[j].set_xlabel('log10(wavelength / 1 micron)',fontsize=16)
					grid[j].set_ylabel('log10(flux / Jy)',fontsize=16)

				# Add Band Plots
				grid[j].annotate(s="",xy=(.25,0),xycoords='data',xytext=(.25,0.3),textcoords='data',
					arrowprops=dict(arrowstyle='-[',connectionstyle='arc3'))
				grid[j].annotate(s="H, K",xy=(.09,0.4),size=12)

				grid[j].annotate(s="",xy=(.7,1.5),xycoords='data',xytext=(.7,1.8),textcoords='data',
					arrowprops=dict(arrowstyle='-[,widthB=1.75',connectionstyle='arc3'))
				grid[j].annotate(s="IRAC",xy=(.5,1.9),size=12)

				grid[j].annotate(s="",xy=(1.35,1.8),xycoords='data',xytext=(1.35,2.1),textcoords='data',
					arrowprops=dict(arrowstyle='-[,widthB=0.6',connectionstyle='arc3'))
				grid[j].annotate(s="MIPS",xy=(1.2,2.2),size=12)

				grid[j].annotate(s="",xy=(1.85,1.8),xycoords='data',xytext=(1.85,2.1),textcoords='data',
					arrowprops=dict(arrowstyle='-[,widthB=0.6',connectionstyle='arc3'))
				grid[j].annotate(s="PACS",xy=(1.7,2.2),size=12)

				# Plot Data
				mask = ~ma.masked_array(guter_fluxes[i],mask=guter_fluxes[i]>1e5).mask
				grid[j].errorbar(np.log10(wavelens[mask]),np.log10(guter_fluxes[i][mask]),
					yerr=.434*guter_rel_errs[mask],fmt='s',color='darkblue')
				#grid[j].set_xscale('log')
				#grid[j].set_yscale('log')
				#grid[j].set_xlabel(r'$\lambda$',fontsize=16)
				#grid[j].set_ylabel(r'$\lambda\cdot S$',fontsize=16)
				grid[j].set_xlim(0,2.2)
				grid[j].set_ylim(-5,5)
			
				#grid[j].grid()	
				grid[j].get_xaxis().set_ticks([0.2,0.6,1.0,1.4,1.8])
				grid[j].get_yaxis().set_ticks([-4,-3,-2,-1,0,1,2,3,4])
				#grid[j].tick_params(which='minor',length=3)	
				#grid[j].tick_params(which='major',length=5)	

				props = dict(boxstyle='round', facecolor='white', alpha=0.8)
				grid[j].axes.text(0.1,0.9,'VLA '+str(i+1),transform=grid[j].axes.transAxes,fontsize=14,bbox=props)

				j += 1

	#	grid[6].axes.set_xticklabels(grid[5].axes.get_xticklabels())
	#	grid[7].axes.set_xticklabels(grid[5].axes.get_xticklabels())
	
		if savefig == True:
			mp.savefig('Paper/figures/Infra_SEDs_OnePlot.eps')
		else:
			mp.show()
	
	else:
		for i in range(source_num):
			if ir_source_id[i] == -99: continue
			print '...working on source',i+1
			fig = mp.figure()
			ax = fig.add_subplot(111)
			mask = ~ma.masked_array(guter_fluxes[i],mask=guter_fluxes[i]>1e5).mask
			ax.errorbar(wavelens[mask],guter_fluxes[i][mask],
				yerr=guter_rel_errs[mask]*guter_fluxes[i][mask],fmt='s',color='darkblue')
			ax.set_xscale('log')
			ax.set_yscale('log')
			ax.set_xlabel(r'Wavelength ($\mu$m)',fontsize=16)
			ax.set_ylabel('Flux Density (Jy)',fontsize=16)
			ax.set_xlim(1,200)
			ax.set_ylim(1e-5,3000)
			ax.set_title('Infrared SED of VLA_'+str(i+1))
			if savefig == True:
				mp.savefig('Infra_SED_Source'+str(i)+'.png')
			else:
				mp.show()
		#	mp.close()	








