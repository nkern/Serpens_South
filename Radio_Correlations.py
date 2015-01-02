## This file loads in Scaling Relationship Data, like L_radio vs. L_bolometric
## It does Chi Square fitting and makes plots

## Import Modules ##
import numpy as np
import matplotlib.pyplot as mp
import astropy.io.fits as fits
from scipy.optimize import leastsq
from AttrDict import AttrDict
from scipy import stats
import pickle as pkl
from matplotlib import rc

## Flags ##


## Constants ##
dist1 = 0.429	# kpc
dist2 = 0.260	# kpc
c = 2.99e8	# m/s

## Load In Data ##

# Load Shirley07 Data
shirley = AttrDict({})
shirley.Lbol_36, shirley.flux_36 = np.loadtxt('Shirley07Data/3.6cmcont.dat',unpack=True)
shirley.Lbol_60, shirley.flux_60 = np.loadtxt('Shirley07Data/6cmcont.dat',unpack=True)

# Load Scaife11a Data
scaifea = AttrDict({})
data = np.loadtxt('Scaife11a_Data.tab',dtype='str',delimiter='\t',usecols=(0,1,2,3,4,5,6,7,8,9))
match = np.where(data.T[4]=='y')[0]
match = np.delete(match,[1,11,15,27])
scaifea.classes = data.T[2][match]
scaifea.flux_18 = np.array(data.T[5][match],float)
scaifea.Menv = data.T[6][match]
scaifea.Tbol = data.T[7][match]
scaifea.Lbol = np.array(data.T[8][match],float)
scaifea.Fout = data.T[9][match]

data2 = np.loadtxt('Scaife11a_Data2.tab',dtype='str',delimiter='\t',usecols=(0,1,2,3,4,5,6,7))
change_dist = np.array(['093','106','104','107','105','063','071',
			'065','064','088','084','060','068','080','090','092','109'])
dist_fix = np.array(map(lambda x:np.where(data2.T[0]==x)[0][0],change_dist))
dist_fix = dist_fix[np.argsort(dist_fix)]
new_dist = np.array(data2.T[5][dist_fix],float)*1e3
scaifea.dist = 0.25	#kpc

# Load Scaife11b Data
scaifeb = AttrDict({})
data = np.loadtxt('Scaife11b_Data.tab',dtype='str',delimiter='\t',usecols=(0,1,2,3,4,5,6,7,8,9))
match = np.where(data.T[3]!='-')[0]
scaifeb.classes = data.T[9][match]
scaifeb.distances = np.array(map(lambda x: x[0:3],data.T[1]),float)[match]/1e3
scaifeb.flux_18 = np.array(data.T[8][match],float)
scaifeb.Lbol = np.array(data.T[3][match],float)

# Merge Scaife Data
scaife = AttrDict({})
scaife.Lbol = np.concatenate([scaifea.Lbol,scaifeb.Lbol])
scaife.flux_18 = np.concatenate([np.log10(scaifea.flux_18/1e3*scaifea.dist**2),scaifeb.flux_18])
scaife.distances = np.concatenate([[scaifea.dist]*scaifea.Lbol.shape[0], scaifeb.distances])

# Load Centimeter Data
data = pkl.Unpickler(open('Radio_Data.pkl','rb')).load()
for i in data:
	try:
		data[i] = np.array(data[i],float)
	except ValueError:
		pass
kern = AttrDict(data)
kern.Lbol_sources = np.array(kern.Lbol_sources,int)

# Load Mallick Data
mallick_data = np.loadtxt("Mallick13_Data.tab",dtype='str',unpack=True,delimiter='\t')
mallick = AttrDict({})
mallick.ra = np.array(map(lambda x: x[0:10], mallick_data))
mallick.dec = np.array(map(lambda x: x[11:20], mallick_data))
mallick.jmag = np.array(map(lambda x: x[21:27], mallick_data))
mallick.jmag_err = np.array(map(lambda x: x[28:33], mallick_data)) 
mallick.hmag = np.array(map(lambda x: x[34:40], mallick_data)) 
mallick.hmag_err = np.array(map(lambda x: x[41:46], mallick_data)) 
mallick.kmag = np.array(map(lambda x: x[47:53], mallick_data)) 
mallick.kmag_err = np.array(map(lambda x: x[54:59], mallick_data)) 
mallick.i1 = np.array(map(lambda x: x[60:66], mallick_data))
mallick.i1_err = np.array(map(lambda x: x[67:72], mallick_data))
mallick.i2 = np.array(map(lambda x: x[73:79], mallick_data))
mallick.i2_err = np.array(map(lambda x: x[80:85], mallick_data))
mallick.i3 = np.array(map(lambda x: x[86:92], mallick_data))
mallick.i3_err = np.array(map(lambda x: x[93:98], mallick_data))
mallick.i4 = np.array(map(lambda x: x[99:105], mallick_data))
mallick.i4_err = np.array(map(lambda x: x[106:111], mallick_data))
mallick.classes = np.array(map(lambda x: x[112:118], mallick_data))


# Load Rodriguez Data
rodriguez_data = np.loadtxt('Rodriguez10_Data.tab',delimiter='\t',unpack=True,dtype='str')
rodrig = AttrDict({})
rodrig.ra = np.array(map(lambda x: float(x[0:2]), rodriguez_data[1]))*360./24. + \
np.array(map(lambda x: float(x[3:5]), rodriguez_data[1]))*360./(24.*60) + \
np.array(map(lambda x: float(x[6:11]), rodriguez_data[1]))*360./(24.*60*60)

rodrig.dec = np.array(map(lambda x: float(x[0:3]), rodriguez_data[2])) - \
np.array(map(lambda x: float(x[4:6]), rodriguez_data[2]))/60. - \
np.array(map(lambda x: float(x[7:12]), rodriguez_data[2]))/(60.*60)

rodrig.flux_36 = np.array(rodriguez_data[3],float)	# mJy

# Can't use Mallick13 or Rodriguez10 data b/c we don't have MIPS 24 micron data from W40, which we need to
# calculate IR Spix to be used in Kryukova12 Lmir to Lbol relationship..


################################
##### Least Square Fitting #####
################################

## Model is log-log
# log(y) = a*log(x) + b	

def residual(vars, x, data):
    """ model for data, subtract data """
    a = vars[0]
    b = vars[1]
    model = a*x + b
    return (data-model)**2


def standard_err(x_data,y_data,a,b):
	n = len(x_data)
	err_a = np.sqrt( (1/np.float(n-2)) * np.sum((y_data - x_data*a - b)**2) / np.sum((x_data - np.median(x_data))**2) )
	err_b = err_a * np.sqrt(np.sum(x_data**2)/n)
	return err_a,err_b



## Scaife 1.8 cm Fit
vars = [0.7,-2.0]	# Guess
scaife_x = np.log10(scaife.Lbol)
scaife_y = scaife.flux_18

x_data = np.copy(scaife_x)
y_data = np.copy(scaife_y)

out = leastsq(residual, vars, args=(x_data, y_data))

scaife_a = out[0][0]
scaife_b = out[0][1]
scaife_a_err, scaife_b_err = standard_err(x_data,y_data,scaife_a,scaife_b)
scaife_xrange = np.arange(scaife_x.min()-.5,scaife_x.max()+.5,.1)
scaife_model = scaife_a*scaife_xrange + scaife_b

#p1, = mp.plot(x_data,y_data,'bo',alpha=.6)
#mp.plot(scaife_x,scaife_model,'b')

## Shirley 3.6cm Fit
vars = [0.75,-2.2]
shirley36_x = np.delete(shirley.Lbol_36,14)
shirley36_y = np.delete(shirley.flux_36,14)

x_data = np.copy(shirley36_x)
y_data = np.copy(shirley36_y)

out = leastsq(residual, vars, args=(x_data, y_data) )

shirley36_a = out[0][0]
shirley36_b = out[0][1]
shirley36_a_err, shirley36_b_err = standard_err(x_data,y_data,shirley36_a,shirley36_b)
shirley36_xrange = np.arange(-2,4,.1)
shirley36_model = shirley36_a*shirley36_xrange + shirley36_b

#p2, = mp.plot(x_data,y_data,'go',alpha=.6)
#mp.plot(shirley36_x,shirley36_model,'g')

## Kern 4.1cm Fit
vars = [0.8,-2.3]
kern_x_kryukova1 = np.log10(kern.Lbol_kryukova1)
kern_x_kryukova2 = np.log10(kern.Lbol_kryukova2)
kern_x_dunham1 = np.log10(kern.Lbol_dunham1)
kern_x_dunham2 = np.log10(kern.Lbol_dunham2)
kern41_y1 = np.log10(kern.flux_upper/1e3*dist1**2)[kern.Lbol_sources]			# mJy * dist^2
kern41_y1_err = (kern.flux_upper_err/kern.flux_upper)[kern.Lbol_sources]
kern41_y2 = np.log10(kern.flux_upper/1e3*dist2**2)[kern.Lbol_sources]
kern41_y2_err = (kern.flux_upper_err/kern.flux_upper)[kern.Lbol_sources]

x_data = np.copy(kern_x_kryukova1)
y_data = np.copy(kern41_y2)

out = leastsq(residual, vars, args=(x_data,y_data) )

kern41_a = out[0][0]
kern41_b = out[0][1]
kern41_a_err, kern41_b_err = standard_err(x_data,y_data,kern41_a,kern41_b)
kern41_xrange = np.arange(-1,1,.1)
kern41_model = kern41_a*kern41_xrange + kern41_b

#p3, = mp.plot(x_data,y_data,'ko',alpha=.6)
#mp.plot(kern41_x,kern41_model,'k')

## Shirley 6.0cm Fit
vars = [0.9, -2.5]
shirley60_x = shirley.Lbol_60
shirley60_y = shirley.flux_60

x_data = np.copy(shirley60_x)
y_data = np.copy(shirley60_y)

out = leastsq(residual,vars,args=(x_data,y_data))

shirley60_a = out[0][0]
shirley60_b = out[0][1]
shirley60_a_err, shirley60_b_err = standard_err(x_data,y_data,shirley60_a,shirley60_b)
shirley60_xrange = np.arange(-2,4,.1)
shirley60_model = shirley60_a*shirley60_xrange + shirley60_b

#p4, = mp.plot(x_data,y_data,'mo',alpha=.6)
#mp.plot(shirley60_x,shirley60_model,'m')

## Kern 6.3cm Fit
vars = [1.0,-2.6]
kern63_y1 = np.log10(kern.flux_lower/1e3*dist1**2)[kern.Lbol_sources]
kern63_y1_err = (kern.flux_lower_err/kern.flux_lower)[kern.Lbol_sources]
kern63_y2 = np.log10(kern.flux_lower/1e3*dist2**2)[kern.Lbol_sources]
kern63_y2_err = (kern.flux_lower_err/kern.flux_lower)[kern.Lbol_sources]

x_data = np.copy(kern_x_kryukova1)
y_data = np.copy(kern63_y1)

out = leastsq(residual, vars, args=(x_data,y_data))

kern63_a = out[0][0]
kern63_b = out[0][1]
kern63_a_err, kern63_b_err = standard_err(x_data,y_data,kern63_a,kern63_b)
kern63_xrange = np.arange(-0.5,0.75,.1)
kern63_model = kern63_a*kern63_xrange + kern63_b

#p5, = mp.plot(x_data,y_data,'co',alpha=.6)
#mp.plot(kern63_x,kern63_model,'c')

#mp.legend([p1,p2,p3,p4,p5],['scaife','shirley36','kern41','shirley60','kern63'])





##################################
######## P L O T T I N G #########
##################################


## Flags
plot1 = True
plot2 = False

## Plot 4.1 cm over Shirley 3.6 cm data and 6.3 cm over Shirley 6.0 cm data
if plot1 == True:

	figname1 = 'Paper/figures/Shirley36_Correlation.eps'
	figname2 = 'Paper/figures/Shirley60_Correlation.eps'
	figname1 = 'AAS_2015/Shirley36_Correlation.eps'
	figname2 = 'AAS_2015/Shirley60_Correlation.eps'

	savefig = True

	# Extrapolate 4.1 cm to 3.6 cm and 6.3 cm to 6.0 cm w/ spectral indices	
	trust_int_spix = np.where(kern.cm_spix_trust == 'int')[0]
	trust_peak_spix = np.where(kern.cm_spix_trust == 'peak')[0]

	kern36_flux1 = []
	kern36_flux_err1 = []
	kern36_flux2 = []
	kern36_flux_err2 = []

	kern60_flux1 = []
	kern60_flux_err1 = []
	kern60_flux2 = []
	kern60_flux_err2 = []

	for i in range(len(kern.Lbol_sources)):
		if kern.Lbol_sources[i] in trust_int_spix:
			kern36_flux1.append(kern41_y1[i] + np.log(8.46e9/7.25e9) * kern.cm_int_spix[i])
			kern36_flux_err1.append(np.sqrt(kern41_y1_err[i]**2 + (np.log(8.46e9/7.25e9) * kern.cm_int_spix_err[i])**2))
			kern36_flux2.append(kern41_y2[i] + np.log(8.46e9/7.25e9) * kern.cm_int_spix[i])
			kern36_flux_err2.append(np.sqrt(kern41_y2_err[i]**2 + (np.log(8.46e9/7.25e9) * kern.cm_int_spix_err[i])**2))

			kern60_flux1.append(kern63_y1[i] + np.log(4.86e9/4.75e9) * kern.cm_int_spix[i])
			kern60_flux_err1.append(np.sqrt(kern63_y1_err[i]**2 + (np.log(4.86e9/4.75e9) * kern.cm_int_spix_err[i])**2))
			kern60_flux2.append(kern63_y2[i] + np.log(4.86e9/4.75e9) * kern.cm_int_spix[i])
			kern60_flux_err2.append(np.sqrt(kern63_y2_err[i]**2 + (np.log(4.86e9/4.75e9) * kern.cm_int_spix_err[i])**2))

		if kern.Lbol_sources[i] in trust_peak_spix:
			kern36_flux1.append(kern41_y1[i] + np.log(8.46e9/7.25e9) * kern.cm_peak_spix[i])
			kern36_flux_err1.append(np.sqrt(kern41_y1_err[i]**2 + (np.log(8.46e9/7.25e9) * kern.cm_peak_spix_err[i])**2))
			kern36_flux2.append(kern41_y2[i] + np.log(8.46e9/7.25e9) * kern.cm_peak_spix[i])
			kern36_flux_err2.append(np.sqrt(kern41_y2_err[i]**2 + (np.log(8.46e9/7.25e9) * kern.cm_peak_spix_err[i])**2))

			kern60_flux1.append(kern63_y1[i] + np.log(4.86e9/4.75e9) * kern.cm_peak_spix[i])
			kern60_flux_err1.append(np.sqrt(kern63_y1_err[i]**2 + (np.log(4.86e9/4.75e9) * kern.cm_peak_spix_err[i])**2))
			kern60_flux2.append(kern63_y2[i] + np.log(4.86e9/4.75e9) * kern.cm_peak_spix[i])
			kern60_flux_err2.append(np.sqrt(kern63_y2_err[i]**2 + (np.log(4.86e9/4.75e9) * kern.cm_peak_spix_err[i])**2))

	# Plot 3.6 correlation
	fig1,ax1 = mp.subplots()

	p1, = ax1.plot(shirley36_x,shirley36_y,'ko',alpha=.8)
	p2 = ax1.errorbar(kern_x_kryukova1,kern36_flux1,xerr=kern.Lbol_kryukova1_err/kern.Lbol_kryukova1,yerr=kern36_flux_err1,fmt='o',color='b',alpha=.6)
#	p3 = ax1.errorbar(kern_x_kryukova2,kern36_flux2,xerr=0.5,yerr=kern36_flux_err2,fmt='o',color='r',alpha=.9)
#	p4 = ax1.errorbar(kern_x_dunham1,kern36_flux1,xerr=0.5,yerr=kern36_flux_err1,fmt='o',color='r',alpha=.6)
#	p5 = ax1.errorbar(kern_x_dunham2,kern36_flux2,xerr=0.5,yerr=kern36_flux_err2,fmt='o',color='m',alpha=.7)
#	mp.legend([p1,p2,p4],["Shirley+07","This work and Kryukova+ 2012","This work and Dunham+ 2013"],loc=2)

	mp.rc('text', usetex=True)
	mp.rc('font', family='serif')
	ax1.set_xlabel(r'log10( L$_{bol}$ / 1 L$_{\odot}$ )',fontsize=20)
	ax1.set_ylabel(r'log10( $S\cdot$D$^{2}$ / 1mJy kpc$^{2}$ )',fontsize=20)

	if savefig == True:
		fig1.savefig(figname1,adjust_bbox='tight',transparent=True)

	# Plot 6.0 corelation
	fig2,ax2 = mp.subplots()

	p1, = ax2.plot(shirley60_x,shirley60_y,'ko',alpha=.8)
	p2 = ax2.errorbar(kern_x_kryukova1,kern60_flux1,xerr=kern.Lbol_kryukova1_err/kern.Lbol_kryukova1,yerr=kern60_flux_err1,fmt='o',color='b',alpha=.6)
#	p3 = ax2.errorbar(kern_x_kryukova2,kern60_flux2,xerr=0.5,yerr=kern60_flux_err2,fmt='o',color='r',alpha=.9)
#	p4 = ax2.errorbar(kern_x_dunham1,kern60_flux1,xerr=0.5,yerr=kern60_flux_err1,fmt='o',color='r',alpha=.6)
#	p5 = ax2.errorbar(kern_x_dunham2,kern60_flux2,xerr=0.5,yerr=kern60_flux_err2,fmt='o',color='m',alpha=.7)
#	mp.legend([p1,p2,p4],["Shirley+07","This work and Kryukova+ 2012","This work and Dunham+ 2013"],loc=2)

	mp.rc('text', usetex=True)
	mp.rc('font', family='serif')
	ax2.set_xlabel(r'log10( L$_{bol}$ / 1 L$_{\odot}$ )',fontsize=20)
	ax2.set_ylabel(r'log10( $S\cdot$D$^{2}$ / 1 mJy kpc$^{2}$ )',fontsize=20)

	if savefig == True:
		fig2.savefig(figname2,adjust_bbox='tight',transparent=True)
	


## Normalized Flux vs. Lbol diagram to 5GHz
if plot2 == True:
	# Iterate over different global spectral indices, normalized to 5 GHz = 5.98 cm
	global_spix = [.6]
	freqs = [16.07e9, 8.46e9, 7.25e9, 4.86e9, 4.75e9]

	# Normalize to 5 GHz, outlier at global_data = 63
	global_spix = 0.6
	scaife_flux = scaife.flux_18 + np.log(5e9/freqs[0]) * global_spix
	shirley36_flux = shirley36_y + np.log(5e9/freqs[1]) * global_spix
	kern41_flux = kern41_y1 + np.log(5e9/freqs[2]) * global_spix
	shirley60_flux = shirley60_y + np.log(5e9/freqs[3]) * global_spix
	kern63_flux = kern63_y1 + np.log(5e9/freqs[4]) * global_spix

	global_ydata = np.concatenate([scaife_flux,shirley36_flux,kern41_flux,shirley60_flux,kern63_flux])
	global_xdata = np.concatenate([scaife_x,shirley36_x,kern_x_kryukova1,shirley60_x,kern_x_kryukova1])

	# Chi Square
	global_x = np.copy(global_xdata)
	global_y = np.copy(global_ydata)
	vars = [1.0,-2.3]
	out = leastsq(residual, vars, args=(global_xdata,global_ydata))

	global_a, global_b = out[0][0], out[0][1]
	global_a_err, global_b_err = standard_err(global_x,global_y,global_a,global_b)
	global_xrange = np.arange(-1.25,3.5,.1)
	global_model = global_a*global_xrange + global_b


	fig,ax = mp.subplots()

	p1, = ax.plot(scaife_x,scaife_flux,'ro',alpha=.6)
	p2, = ax.plot(kern_x_kryukova1,kern41_flux,'bo',alpha=.6)
	p3, = ax.plot(shirley60_x,shirley60_flux,'go',alpha=.6)
	
	fig.legend([p1,p2,p3],["Scaifea+11;Scaifeb+11","This Work","Shirley+07"])




