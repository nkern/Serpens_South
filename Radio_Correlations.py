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
scaifeb.distances = np.array(map(lambda x: x[0:3],data.T[1]),float)/1e3
scaifeb.flux_18 = np.array(data.T[8][match],float)
scaifeb.Lbol = np.array(data.T[3][match],float)

# Merge Scaife Data
scaife = AttrDict({})
scaife.Lbol = np.concatenate([scaifea.Lbol,scaifeb.Lbol])
scaife.flux_18 = np.concatenate([np.log10(scaifea.flux_18/1e3*scaifea.dist**2),scaifeb.flux_18])

# Load Centimeter Data
data = pkl.Unpickler(open('Radio_Data.pkl','rb')).load()
for i in data:
	try:
		data[i] = np.array(data[i],float)
	except ValueError:
		pass
kern = AttrDict(data)













##############################
##### Chi Square Fitting #####
##############################

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
x_data = np.log10(scaife.Lbol)
y_data = scaife.flux_18
out = leastsq(residual, vars, args=(x_data, y_data))

scaife_a = out[0][0]
scaife_b = out[0][1]
scaife_a_err, scaife_b_err = standard_err(x_data,y_data,scaife_a,scaife_b)
scaife_x = np.arange(-1,1.3,.1)
sciafe_model = scaife_a*scaife_x + scaife_b

## Shirley 3.6cm Fit
vars = [0.7,-2.2]
#x_data = 
#y_data = 
#out = leastsq( )

shirley36_a = out[0][0]
shirley36_b = out[0][1]
shirley36_a_err, shirley36_b_err = standard_err(x_data,y_data,shirley36_a,shirley36_b)
shirley36_x = np.arange(-2,4,.1)
shirley36_model = shirley36_a*shirley36_x + shirley36_b

## Kern 4.1cm Fit
vars = [0.8,-2.3]
x_data1 = np.log10(kern.flux_upper/1e3*dist1**2)	# mJy * dist^2
x_data2 = np.log10(kern.flux_upper/1e3*dist2**2)
y_data_kryukova_rising1 = np.log10(kern.Lbol_kryukova_rising1)
y_data_kryukova_rising2 = np.log10(kern.Lbol_kryukova_rising2)
y_data_kryukova_flat1 = np.log10(kern.Lbol_kryukova_flat1)
y_data_kryukova_flat2 = np.log10(kern.Lbol_kryukova_flat2)
y_data_dunham_rising1 = np.log10(kern.Lbol_dunham_rising1)
y_data_dunham_rising2 = np.log10(kern.Lbol_dunham_rising2)
y_data_dunham_flat1 = np.log10(kern.Lbol_dunham_flat1)
y_data_dunham_flat2 = np.log10(kern.Lbol_dunham_flat2)

raise NameError
out = leastsq()

kern41_a = out[0][0]
kern41_b = out[0][1]
kern41_a_err, kern41_b_err = standard_err(x_data,y_data,kern41_a,kern41_b)
kern41_x = np.arange()
kern41_model = kern41_a*kern41_x + kern41_b

## Shirley 6.0cm Fit
vars = [0.9,-2.5]
#x_data = 
#y_data = 
out = leastsq()

shirley60_a = out[0][0]
shirley60_b = out[0][1]
shirley60_a_err, shirley60_b_err = standard_err(x_data,y_data,shirley60_a,shirley60_b)
shirley60_x = np.arange(-2,4,.1)
shirley60_model = shirley60_a*shirley60_x + shirley60_b

## Kern 6.3cm Fit
vars = [1.0,-2.6]
#x_data = 
#y_data = 
#out = leastsq()

kern63_a = out[0][0]
kern63_b = out[0][1]
kern63_a_err, kern63_b_err = standard_err(x_data,y_data,kern63_a,kern63_b)
kern63_x = np.arange()
kern63_model = kern63_a*kern63_x + kern63_b


## Combined Fit
# Iterate over different global spectral indices, normalized to 5 GHz = 5.98 cm
global_spix = np.arange(-0.1,1.1,0.2)






