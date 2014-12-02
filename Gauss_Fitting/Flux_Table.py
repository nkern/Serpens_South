'''
This script takes data and prints out Latex Table strings to be copy-pasted into tablular float in master paper.tex file
It is dependent on Table.py and sigfig.py
'''
import numpy as np
import Table
from wcsaxes import WCS
import astropy.io.fits as fits

SN = 18		# Source Number

#################
### Load Data ###
#################
print '...Loading Data'

## Load flux data
sourcenames = np.array(map(lambda x:"VLA "+x,np.array(np.arange(1,SN+1),str)))
sources = np.loadtxt('Final_Fit_Results_try2.txt',usecols=(0,),dtype='str')
flux,flux_err,ra,dec,ra_err,dec_err,bmaj,bmin,pos_angle,rms=np.loadtxt('Final_Fit_Results_try2.txt',usecols=(1,2,3,4,5,6,7,8,9,10),unpack=True)
peak_flux = []
j = 0
for i in np.concatenate([np.arange(1,SN+1),np.arange(1,SN+1)]):
	if j < SN+1:
		peak_flux.append(np.loadtxt('Fit_Files/SERPS_lower_try2.S'+str(i)+'-gaussest_new.txt',delimiter=',',usecols=(0,)))
	elif j > SN:
		peak_flux.append(np.loadtxt('Fit_Files/SERPS_upper_try2.S'+str(i)+'-gaussest_new.txt',delimiter=',',usecols=(0,)))
peak_flux = np.array(peak_flux)
print "Using Peak Fluxes from IMFIT's 2D Gaussian Estimates"

# Put flux data into micro jansky
flux *= 1e6
flux_err *= 1e6
rms *= 1e6

# Configure Flux Data
lower_flux,upper_flux = flux[0:SN],flux[SN:]
lower_flux_err,upper_flux_err = flux_err[0:SN],flux_err[SN:]
lower_peak_flux,upper_peak_flux = peak_flux[0:SN],peak_flux[SN:]
lower_peak_flux_err, upper_peak_flux_err = 11.1e-5, 8.5e-6
lower_ra,upper_ra = ra[0:SN],ra[SN:]
lower_dec,upper_dec = dec[0:SN],dec[SN:]
lower_ra_err,upper_ra_err = ra_err[0:SN],ra_err[SN:]
lower_dec_err,upper_dec_err = dec_err[0:SN],dec_err[SN:]
lower_bmaj,upper_bmaj = bmaj[0:SN],bmaj[SN:]
lower_bmin,upper_bmin = bmin[0:SN],bmin[SN:]
lower_posang,upper_posang = pos_angle[0:SN],pos_angle[SN:]
lower_rms,upper_rms = rms[0:SN],rms[SN:]

# Convert RA and DEC to J2000
ra = 360 + ra 
ra_frac = ra/360. * 24
ra_h = np.array(np.floor(ra_frac),int)
ra_m = np.array(np.floor((ra_frac - ra_h)*60),int)
ra_s = np.around(((ra_frac - ra_h)*60 - ra_m)*60,1)
ra_s = np.array(map(lambda x:"%02g" % (x,), ra_s))
ra_sexig = np.array([a+':'+b+':'+c for a,b,c in zip(np.array(ra_h,str),np.array(ra_m,str),np.array(ra_s,str))])

dec_deg = np.ceil(dec)
dec_m = np.array(np.floor((dec_deg-dec)*60),int)
dec_s = np.around( ((dec_deg-dec)*60 - dec_m)*60,1)
dec_deg = np.array(map(lambda x: "%03d" % (x), dec_deg))
dec_m = np.array(map(lambda x:"%02d" % (x), dec_m))
dec_s = np.array(map(lambda x:"%02g" % (x), dec_s))
for i,k in zip(dec_s,np.arange(len(dec_s))): 
	if np.float(i)/10 < 1: dec_s[k]='0'+dec_s[k]
	if np.float(i)%1.0 == 0: dec_s[k]=dec_s[k]+'.0'
dec_sexig = np.array([a+':'+b+':'+c for a,b,c in zip(np.array(dec_deg,str),np.array(dec_m,str),np.array(dec_s,str))])


#######################################
### Print Out Latex Tabular Strings ###
#######################################


#	Order,		Array,			Units :
#	1 Source,	sourcenames		string
#	2 RA,		ra_sexig		sexigismal
#	3 Dec.,		dec_sexig		sexigismal
#	4 S_6.3,	lower_flux		microJansky float w/ error
#	5 S^peak_6.3,	lower_peak_flux		microJansky/beam float w/ error
#	6 S_4.1,	upper_flux		microJansky float w/ error
#	7 S^peak_4.1,	upper_peak_flux		microJanksy/beam float w/ error
#	8 spix,		cm_int/peak_spix	float w/ error







# Set Table Constants
#col_num	= 6						# Number of Columns in Table
#justs	= 'l'+'c'*(col_num-1)				# Justifications (left, center, right)
#title	= 'Radio Properties of VLA Sources'		# Caption to Table
#label	= 'flux_table'					# Label
#header	= ['Source','$S_{6.3 cm}$','$S_{4.1 cm}$','RA','Dec.']#,'Bmaj','Bmin','Pos. Angle','RMS']
#header2	= ['','(mJy)','(mJy)','(mJy)','(J2000)','(J2000)']#,'(arcsec)','(arcsec)','(degrees)','(mJy)']
#keys	= ['sourcenames','flux','flux_err','peak_flux','ra','dec']#,'bmaj','bmin','pos_angle','rms']
# Create Table
#make_table = False
#if make_table == True:
#	print '...Writing Latex Table'
#	fout = open('Fit_Latex_Table2.tex','w')
#	fout.write("\\begin{landscape}\n\\tabletypesize{\\footnotesize}\n")
#	t = Table.Table(col_num, justs=justs, caption=title, label=label)
#	t.add_header_row(header,header2)
#	cols = []
#	for i in range(col_num):
#		cols.append(eval(keys[i]))
#	t.add_data(cols, sigfigs=4)
#	t.print_table(fout)
#	fout.write("\end{landscape}")
#	fout.close()

