'''
This script takes data and writes a Latex Table.
It is dependent on Table.py and sigfig.py
'''
import numpy as np
import Table

Source_Num = 21

# Turn into Latex Table (need Table.py and sigfig.py in working directory)
print '...Loading Data'

# Load and Configure Necessary Data for Table
sourcename = np.loadtxt('Final_Fit_Results.txt',dtype=str,usecols=(0,),unpack=True)
sourcenames = np.array(map(lambda x:x[6:],sourcename))
flux,flux_err,ra,dec,ra_err,dec_err,bmaj,bmin,pos_angle,rms=np.loadtxt('Final_Fit_Results.txt',usecols=(1,2,3,4,5,6,7,8,9,10),unpack=True)
peak_flux = []
j = 0
for i in np.concatenate([np.arange(1,Source_Num+1),np.arange(1,Source_Num+1)]):
	j += 1
	if j < 22:
		peak_flux.append(np.loadtxt('SERPS_lower.S'+str(i)+'-gaussest.txt',delimiter=',',usecols=(0,)))
	elif j > 21:
		peak_flux.append(np.loadtxt('SERPS_upper.S'+str(i)+'-gaussest.txt',delimiter=',',usecols=(0,)))
peak_flux = np.array(peak_flux)

# Convert flux values to mJy
flux *= 1e3
flux_err *= 1e3
peak_flux *= 1e3
rms *= 1e3

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
dec_sexig = np.array([a+':'+b+':'+c for a,b,c in zip(np.array(dec_deg,str),np.array(dec_m,str),np.array(dec_s,str))])

# Set Table Constants
col_num	= 12			# Number of Columns in Table
justs	= 'l'+'c'*(col_num-1)	# Justifications (left, center, right)
caption = 'Caption...'	# Caption to Table
label	= 'Properties of 2D Gaussian Fits'	# Table Title
header	= ['Source Name','Int. Flux','Int. Flux Err','Peak Flux','RA','RA Err','DEC','DEC Err','Bmaj','Bmin','Pos. Angle','RMS']
header2	= ['','(mJy)','(mJy)','(mJy)','(J2000)','(arcsec)','(J2000)','(arcsec)','(arcsec)','(arcsec)','(degrees)','(mJy)']
keys	= ['sourcenames','flux','flux_err','peak_flux','ra','ra_err','dec','dec_err','bmaj','bmin','pos_angle','rms']

# Create Table
make_table = False
if make_table == True:
	print '...Writing Latex Table'
	fout = open('Fit_Latex_Table.tex','w')
	fout.write("\\begin{landscape}\n")
	t = Table.Table(col_num, justs=justs, caption=caption, label=label)
	t.add_header_row(header,header2)
	cols = []
	for i in range(col_num):
		cols.append(eval(keys[i]))

	t.add_data(cols, sigfigs=4)
	t.print_table(fout)
	fout.write("\end{landscape}")
	fout.close()

