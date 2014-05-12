"""
This code interactively takes in data to write out estimates files needed by SERPS_Gauss_Fit.py
"""
import numpy as np

# Define Constants
Source_Num = 21		# Number of Sources in Field of View, see SS_Targets.reg
Image_Num = 2		# Number of Images with same sources

sources = []
source_id=np.array(np.arange(1,Source_Num+1),str)
for p in range(Source_Num*Image_Num):
	if p < 21:
		sources.append('SERPS_lower.S'+source_id[p])	
	else:
		sources.append('SERPS_upper.S'+source_id[p-21])

sources = np.array(sources)

for i in np.arange(0,Source_Num*Image_Num):
	file = open(sources[i]+'-gaussest.txt','r')
	cont = False
	while cont == False:
		print ''
		print 'Working on Source:',sources[i]
		print '-'*25
		flux = raw_input('Flux density of center pixel? (Jy) : ')
		x = raw_input('X pixel coordinate of center pixel? : ')
		y = raw_input('Y pixel coordinate of center pixel? : ')
		major_axis = '3.15arcsec' #raw_input('Major axis of beam? (#.##arcsec) : ')
		minor_axis = '2.52arcsec' #raw_input('Minor axis of beam? (#.##arcsec) : ')
		position_angle = '12.7deg' #raw_input('Position angle of beam? (###.#deg) : ')
		print 'Given Values:'
		print 'flux =',flux
		print 'x =',x
		print 'y =',y
		print 'major_axis =',major_axis
		print 'minor_axis =',minor_axis
		print 'position_angle =',position_angle
		read = raw_input('Is this correct? (y/n) : ')
		if read == 'y':
			cont = True
		else:
			print 'Try again...'	

	print '...Writing Data for Source:',sources[i]
	file.write(str(flux)+','+str(x)+','+str(y)+','+str(major_axis)+','+str(minor_axis)+','+str(position_angle)+'\n')
	file.close()







