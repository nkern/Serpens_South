"""
 - This code takes the two SERPS.**.fits files and uses imfit-sources.py to fit 
Gaussians to sources in the images.

 - To be run within CASA: execfile('SERPS_Gauss_Fit.py')

 - Dependent on readcol.py
	http://code.google.com/p/agpy/source/browse/trunk/agpy/readcol.py

May, 2014
"""

### Import Modules ###
print '...Loading Modules'
print '-'*30
import numpy as np
import matplotlib.pyplot as mp
from readcol import *
import os, sys


### Define Flags ###


### Functions ###
def extractcomponents(sourcename=None,componentlist=None,comp=None,rms=None):
	print comp
	#outputfilefits = open(sourcename+'-fitresults.txt','w')
	#outputfileseps = open(sourcename+'-separations.txt','w')
	cl.open(componentlist)
	separations=np.zeros(comp, dtype=float)
	e_separations=np.zeros(comp, dtype=float)
	rakeep=np.zeros(comp, dtype=float)
	deckeep=np.zeros(comp, dtype=float)
	e_rakeep=np.zeros(comp, dtype=float)
	e_deckeep=np.zeros(comp, dtype=float)

	for j in range(comp):
		fit=cl.getcomponent(j)
		flux = fit['flux']['value']
		e_flux = fit['flux']['error']
		print flux[0]
		ra = fit['shape']['direction']['m0']['value']*180.0/3.14159265
		dec = fit['shape']['direction']['m1']['value']*180.0/3.14159265
		e_ra = fit['shape']['direction']['error']['latitude']['value']
		e_dec = fit['shape']['direction']['error']['longitude']['value']
		bmaj = fit['shape']['majoraxis']['value']*60.0     
		bmin = fit['shape']['minoraxis']['value']*60.0       
		pa =   fit['shape']['positionangle']['value']
		print ra,dec,e_ra,e_dec,bmaj,bmin
		line = str(sourcename[:])+" "+str(flux[0])+" "+str(e_flux[0])+ " " +str(ra)+" "+str(dec)+" "+str(e_ra)+" "+str(e_dec)+" " +str(bmaj)+" "+str(bmin)+" "+str(pa)+" "+str(rms)+"\n"
		rakeep[j]=ra
		deckeep[j]=dec
		e_rakeep[j]=e_ra
		e_deckeep[j]=e_dec
      
		outputfilefits.writelines(line)
	print range(comp)
	iter=0
	for j in range(comp):
		for k in range(comp):
			print e_rakeep[j], e_deckeep[j],e_rakeep[k], e_deckeep[k]
			if j != k and k > j:
				print iter
				dra=rakeep[j]* 3600.0-rakeep[k]* 3600.0
				ddec=deckeep[j]* 3600.0-deckeep[k]* 3600.0

				e_dra=(e_rakeep[j]**2 + e_rakeep[k]**2)**0.5
				e_ddec=(e_deckeep[j]**2 + e_deckeep[k]**2)**0.5

				separations[iter]=sqrt( (dra)**2 + (ddec)**2)  
				e_rasq=(dra*e_dra)**2
				e_decsq=(ddec*e_ddec)**2

				print e_rasq, e_decsq
				e_separations[iter] =  sqrt( e_rasq + e_decsq )/separations[iter]
				iter += 1

	#for j in range(iter):
		#line = str(separations[j])+ " "+str(e_separations[j])+ " "+str(separations[j]*dist) + " "+str(e_separations[j]*dist) +"\n"
		# outputfileseps.writelines(line)
	cl.close()
	#outputfilefits.close()
	#outputfileseps.close()
	print 'end of extract'



#######################
####### Program #######
#######################

## Constants ##
Source_Num = 18		# Number of Sources in Field of View, see SS_Targets.reg
Image_Num = 2		# Number of Images with same sources

run_program = True
if run_program == True:
	print ''
	print '...Beginning Program'
	print '-'*30

	# Define Constants
	# Initialize imfit function and image & source list
	inp(taskname='imfit')

	cl.close()

	# Get Image List
	imagelist=np.concatenate([['../FITS/SERPS_2.lower.briggs05.deeper.image.pbcorr.fits']*Source_Num,['../FITS/SERPS_2.upper.briggs05.deeper.image.pbcorr.fits']*Source_Num])

	# Get Source List
	sources = []
	source_id=np.array(np.arange(1,Source_Num+1),str)
	for p in range(Source_Num*Image_Num):
		if p < Source_Num:
			sources.append('SERPS_lower_try2.S'+source_id[p])	
		else:
			sources.append('SERPS_upper_try2.S'+source_id[p-Source_Num])
	sources = np.array(sources)

	# Get Region List
	regions=[]
	with open('SS_Final_Sources.crtf') as fp:
		for line,i in zip(fp,range(Source_Num+1)):
			if i == 0:
				continue
			regions.append(line[0:-1])
	regions = np.array(regions*2)


	# Open File to Write Data
	outputfilefits = open('Final_Fit_Results_try2.txt','w')
	outputfilefits.write('#sourcename, flux, flux_err, ra, dec, ra_err, dec_err, bmaj, bmin, pos_angle, rms\n')
	# Start Loop Over Sources
	for i in range(len(sources)):
		print ''
		print '...Working on Source:',sources[i]
		print '-'*25

		# Define Filenames	
		inp(taskname='imfit')
		region		= regions[i]
		estimates	= 'Fit_Files/'+sources[i]+'-gaussest.txt'
		imagename	= imagelist[i]
		complist	= 'Fit_Files/'+sources[i]+'-fitcomps.cl'
		overwrite	= True
		newestimates	= 'Fit_Files/'+sources[i]+'-gaussest_new.txt'

		while True:
			inp(taskname='imfit')
			result = imfit()
			if result == None:
				region_size = int(region[39:40])
				new_size = region_size - 1
				region = region.replace(str(region_size)+'arcsec',str(new_size)+'arcsec')

				if new_size == 0:	
					result = {'converged':[False]}
					outputfilefits.write(sources[i]+' '+'0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0\n')
					break
			else:
				break
		print ''
		print '... did fit converge:',result['converged'][0]
		ia.open(imagename)
		stats=ia.statistics()
		rms=stats["rms"][0]
		ia.close()

		print '...running extractcomponents'
		extractcomponents(sources[i],complist,1,rms)



	outputfilefits.close()

### End main program ###


### Turn Final_Fit_Results.txt into ds9 and CASA region file ###
write_files = True
if write_files == True:
	print '...Beginning Region File Write'
	# Load Fit Data File
	sourcename = np.loadtxt('Final_Fit_Results_try2.txt',dtype=str,usecols=(0,),unpack=True)
	flux,flux_err,ra,dec,ra_err,dec_err,bmaj,bmin,pos_angle,rms=np.loadtxt('Final_Fit_Results_try2.txt',usecols=(1,2,3,4,5,6,7,8,9,10),unpack=True)

	pos_angle = pos_angle-90

	# Write out files for lower
	casa_file = open('GaussFits_lower.crtf','w')
	casa_file.write('#CRTFv0 CASA Region Text Format version 0\n')
	ds9_file = open('GaussFits_lower.reg','w')
	ds9_file.write('# Region file format: DS9 version 4.1\n')
	ds9_file.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
	ds9_file.write('fk5\n')

	for i in np.arange(0,Source_Num):
		if flux[i] == 0.0:	# Skip null values marked as 0.0
			continue
		print '...working on source:',sourcename[i]

		casa_file.write('ellipse [['+str(np.around(ra[i],5))+', '+str(np.around(dec[i],5))+'], ['+str(np.around(bmaj[i],5))+'arcsec, '+str(np.around(bmin[i],5))+'arcsec], '+str(np.around(pos_angle[i],5))+'deg] coord=J2000, corr=[I], linewidth=1, linestyle=-, symsize=1, symthick=1, color=green, font=Helvetica, fontsize=10, fontstyle=bold, usetex=false, label="'+str(sourcename[i][13:])+'", labelcolor=green, labelpos=top\n')
	 	
		ds9_file.write('ellipse('+str(np.around(ra[i],5))+','+str(np.around(dec[i],5))+','+str(np.around(bmaj[i],5))+'",'+str(np.around(bmin[i],5))+'",'+str(np.around(pos_angle[i],5))+') # text={'+str(sourcename[i][17:])+'}\n')
	
	casa_file.close()
	ds9_file.close()

	# Write out files for upper
	casa_file = open('GaussFits_lower.crtf','w')
	casa_file.write('#CRTFv0 CASA Region Text Format version 0\n')
	ds9_file = open('GaussFits_upper.reg','w')
	ds9_file.write('# Region file format: DS9 version 4.1\n')
	ds9_file.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
	ds9_file.write('fk5\n')

	for i in np.arange(Source_Num,2*Source_Num):
		if flux[i] == 0.0:	# Skip null values marked as 0.0
			continue
		print '...working on source:',sourcename[i]

		casa_file.write('ellipse [['+str(np.around(ra[i],5))+', '+str(np.around(dec[i],5))+'], ['+str(np.around(bmaj[i],5))+'arcsec, '+str(np.around(bmin[i],5))+'arcsec], '+str(np.around(pos_angle[i],5))+'deg] coord=J2000, corr=[I], linewidth=1, linestyle=-, symsize=1, symthick=1, color=green, font=Helvetica, fontsize=10, fontstyle=bold, usetex=false, label="'+str(sourcename[i][13:])+'", labelcolor=green, labelpos=top\n')
	 	
		ds9_file.write('ellipse('+str(np.around(ra[i],5))+','+str(np.around(dec[i],5))+','+str(np.around(bmaj[i],5))+'",'+str(np.around(bmin[i],5))+'",'+str(np.around(pos_angle[i],5))+') # text={'+str(sourcename[i][17:])+'}\n')
	
	casa_file.close()
	ds9_file.close()

	




###################
### End Program ###
###################










