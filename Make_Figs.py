## Make Paper Figures
import numpy as np
import matplotlib.pyplot as mp
import astropy.io.fits as fits
from astropy.wcs import WCS
from wcsaxes import WCSAxes
from astropy import units as u
import DictEZ as dez
import aplpy
from PIL import Image
import matplotlib
from astropy import coordinates as coord
from astropy import units as u
from matplotlib.widgets import Slider
from AttrDict import AttrDict
import re


## Constants
c = 2.99e8		# m / s
au_per_pc = 206264.806	# au


## Flags
savefig = False

plot1 = False
plot2 = False
plot3 = True
plot4 = False

def sexig_to_deg(ra,dec):
	length = len(ra)
	ra_deg = []
	dec_deg = []
	for i in range(length):
		_ra_deg = np.float(ra[i][0:2])*360./24. + np.float(ra[i][3:5])*360./24./60. + np.float(ra[i][6:10])*360./24./60./60.
		_dec_deg = (np.float(dec[i][1:3]) + np.float(dec[i][4:6])/60. + np.float(dec[i][7:11])/3600.)*np.float(dec[i][0]+'1')
		ra_deg.append(_ra_deg)
		dec_deg.append(_dec_deg)
	ra_deg = np.array(ra_deg)
	dec_deg = np.array(dec_deg)
	return ra_deg,dec_deg

def get_zorder(layers):
	''' Gets zorder of a aplpy figure given dictionary *layers* (fig._layers) '''
	zorder_values = {}
	zorder_paths = {}
	for a in layers:
		try:
			zorder = layers[a].get_zorder()
			path = layers[a]
		except AttributeError:
			try:
				zorder = layers[a].collections[0].get_zorder()
				path = layers[a].collections[0]
			except AttributeError:
				try:
					zorder = layers[a].artistlist[0].get_zorder()
					path = layers[a].artistlist[0]	
				except:
					print "Couldn't get zorder for layer:",a
					zorder = np.nan
					paths = np.nan
		zorder_values[a] = zorder
		zorder_paths[a] = path
	return zorder_values, zorder_paths


### Herschel map w/ Spitzer sources and VLA primary beam
if plot1 == True:

	figname = 'AAS_2015/SerpSouth_FOV.eps'

	savefig = False

	# Choose Image
	color_im = 'FITS/SPIRE_350_micron.fits'
	contour_im = 'FITS/SPIRE_350_micron.fits'

	# Load Header and FITS Data for Contour
	contour_head = fits.open(contour_im)[0].header
	contour_data = fits.open(contour_im)[1].data
	
	# Get RMS
	rms = np.std(contour_data[~np.isnan(contour_data)])

	# Initialize Upper Subband Plot
	fig = aplpy.FITSFigure(color_im, subplot=(1,1,1), figsize=(11,9))

	fig.show_colorscale(cmap='Greys',stretch='power',exponent=.25,pmin=0,pmax=100.001)
	fig.show_contour(contour_im,alpha=1,cmap='Greys_r',levels=[365])

	## Plot Regions
	# Add Spitzer Infrared Sources
	add_spitz = True
	if add_spitz == True:
		spitzer = np.loadtxt('gutermuth2008_serps_ysos1.txt',unpack=True,dtype='str')
		spitz_radec = np.array([np.array(spitzer[1],float),np.array(spitzer[2],float)])
		classI = spitzer[3]=='I'
		fig.show_markers(spitz_radec[0][classI],spitz_radec[1][classI],s=75,marker='^',facecolor='blue',edgecolor='blue',lw=1,zorder=1,alpha=0.9)
		fig.show_markers(spitz_radec[0][~classI],spitz_radec[1][~classI],s=75,marker='^',facecolor='red',edgecolor='red',lw=1,zorder=1,alpha=0.9)

	# Box Regions for Figures
	make_boxes = True
	if make_boxes == True:
		fig.show_rectangles(277.515,-2.040,0.065,0.065,color='Lime',linewidth=5)
		#fig.add_label(277.535,-2.062,'(a)',clip_on=True)

	# Miscellaneous
	fig.ticks.set_color('#000000')
	fig.ticks.set_linewidth(0.75)
	fig.ticks.set_length(8)
	fig.ticks.set_minor_frequency(0)
	fig.ticks.show_x()
	fig.ticks.show_y()
	#fig.axis_labels.set_xtext('Right Ascension (J2000)')
	#fig.axis_labels.set_ytext('Declination (J2000)')
	fig.axis_labels.set_xtext('')
	fig.axis_labels.set_ytext('')
	fig._ax1.set_xticklabels([''])
	fig._ax1.set_yticklabels([''])
	fig.axis_labels.set_font(size=0)
	fig.axis_labels.set_xpad(0)
	fig.axis_labels.set_ypad(0)
	fig._figure.axes[0].tick_params(labelsize=0)

	# Add Scalebar
	fig.add_scalebar(5./60)
	fig.scalebar.set_corner('bottom right')
	fig.scalebar.set_frame(True)
	fig.scalebar.set(linewidth=2,color='black')
	fig.scalebar.set_label('0.6 pc')
	fig.scalebar.set_font(size=18)
	fig.refresh()

	fig.recenter(277.520,-2.060,width=0.225,height=0.325)

	if savefig == True:
		fig.save(figname,adjust_bbox='tight',transparent=False)


### Radio Contours
if plot2 == True:

	zoom = False
	vla4 = True
	vla10 = False
	vla12 = False

	figname = 'Paper/figures/SerpSouth_VLA4_contour.eps'
	#figname = 'Paper/figures/SerpSouth_zoom_contour.eps'

	savefig = True

	# Choose Image
	color_im = 'FITS/PACS_100_micron_cropped.fits'
	contour_im = 'FITS/SERPS_2.upper.briggs05.deeper.fits'

	# Load Header and FITS Data for Contour
	contour_head = fits.open(contour_im)[0].header
	contour_data = fits.open(contour_im)[0].data
	
	contour_data = contour_data[0][0]

	# Get RMS
	rms = np.std(contour_data[~np.isnan(contour_data)])

	# Initialize Upper Subband Plot
	fig = aplpy.FITSFigure(color_im, subplot=(1,1,1), figsize=(11,9))

	fig.show_colorscale(cmap='CMRmap',stretch='power',exponent=.30,pmin=0,pmax=100.005)
	fig.hide_colorscale()
	fig.show_contour(contour_im,alpha=1,cmap='Greys_r',levels=[rms*4,rms*5,rms*6,rms*8,rms*10,rms*15,rms*20,rms*30,rms*50])

	if vla4 == True or vla10 == True or vla12 == True:
		fig.hide_layer('contour_set_1')
		fig.show_contour(contour_im,alpha=0.6,cmap='Blues_r',levels=[rms*3,rms*4,rms*5,rms*6,rms*8,rms*10,rms*15,rms*20,rms*30,rms*50])
		contour_im = 'FITS/SERPS_2.lower.briggs05.deeper.fits'
		contour_data = fits.open(contour_im)[0].data
		rms = np.std(contour_data[~np.isnan(contour_data)])
		fig.show_contour(contour_im,alpha=0.6,cmap='Reds_r',levels=[rms*3,rms*4,rms*5,rms*6,rms*8,rms*10,rms*15,rms*20,rms*30,rms*50],linestyles='dashed')

	# Draw Beam
	beam_dict = {'bmaj':contour_head['bmaj'],'bmin':contour_head['bmin'],'bpa':contour_head['bpa']}
	fig._header.update(beam_dict)
	fig.show_beam(facecolor='black',frame=True)

	# Miscellaneous
	fig.ticks.set_color('#000000')
	fig.ticks.set_linewidth(0.75)
	fig.ticks.set_length(8)
	fig.ticks.set_minor_frequency(5)
	fig.ticks.show_x()
	fig.ticks.show_y()
	fig.axis_labels.set_xtext('Right Ascension (J2000)')
	fig.axis_labels.set_ytext('Declination (J2000)')
	fig.axis_labels.set_xtext('')
	fig.axis_labels.set_ytext('')
	fig.axis_labels.set_font(size=18)
	fig.axis_labels.set_xpad(15)
	fig.axis_labels.set_ypad(10)
	fig._figure.axes[0].tick_params(labelsize=16)
	fig.refresh()

	fig.recenter(277.520,-2.048,width=0.055,height=0.055)

	## Draw Regions ##

	# Box Regions for Figures
	make_boxes = False
	if make_boxes == True:
		fig.show_rectangles(277.518,-2.043,0.04,0.045,color='red',linewidth=1)
		#fig.add_label(277.535,-2.062,'(a)',clip_on=True)
		
	# Add Gaussian Fits
	add_gaussfits = False
	if add_gaussfits == True:
		fig.show_regions('Gauss_Fitting/GaussFits_lower.reg',zorder=3)
		fig.hide_layer('region_set_1_txt')
		#fig.hide_layer('region_set_1')

	# Add IRAM Millimeter sources
	add_iram = True
	if add_iram == True:
		iram_source_name, iram_ra, iram_dec, iram_Lbol, iram_class, iram_S, iram_S_err, iram_FWHM = np.loadtxt('IRAM_sources.txt',unpack=True,delimiter=',',dtype='str')
		iram_ra_deg, iram_dec_deg = sexig_to_deg(iram_ra, iram_dec)
		fig.show_markers(iram_ra_deg,iram_dec_deg,s=550, marker='x', c='green',zorder=1,alpha=.7,lw=1)

	# Add Spitzer Infrared Sources
	add_spitz = True
	if add_spitz == True:
		spitzer = np.loadtxt('gutermuth2008_serps_ysos1.txt',unpack=True,dtype='str')
		spitz_radec = np.array([np.array(spitzer[1],float),np.array(spitzer[2],float)])
		classI = spitzer[3]=='I'
		fig.show_markers(spitz_radec[0][classI],spitz_radec[1][classI],s=550,marker='x',facecolor='None',edgecolor='blue',lw=1,zorder=1,alpha=1)
		fig.show_markers(spitz_radec[0][~classI],spitz_radec[1][~classI],s=550,marker='x',facecolor='None',edgecolor='red',lw=1,zorder=1,alpha=1)

	# Add Teixeira Sources
	add_teixeira = False
	if add_teixeira == True:
		fig.show_regions('Teixeira12_Sources.reg')	

	# Add VLA Labels
	add_labels = True
	if add_labels == True:
		label_ra,label_dec,label_num = np.loadtxt('SS_Labels.txt',delimiter=',',unpack=True)
		label_num = np.array(label_num,int)
		for i in range(len(label_num)):
			fig.add_label(label_ra[i],label_dec[i],'VLA '+str(label_num[i]),clip_on=True,zorder=5,size=18)

	# Add arrows 
	if vla4 == True:
		pass
	if vla10 == True:
		fig.show_arrows(277.5116,-2.0391,0,-0.00045,width=.15,head_width=.5,head_length=.4,facecolor='k')
		fig.show_arrows(277.5093,-2.0421,0.0015,0,width=.15,head_width=.5,head_length=.4,facecolor='k')
		fig.show_arrows(277.5125,-2.0468,0.00075,0.0005,width=.15,head_width=.5,head_length=.4,facecolor='k')

	elif vla12 == True:
		fig.hide_layer('label_11')
		fig.show_arrows(277.5195,-2.0500,-0.0015,-0.0005,width=.15,head_width=.5,head_length=.4,facecolor='k')
		fig.show_arrows(277.5122,-2.0525,0.0011,0,width=.15,head_width=.5,head_length=.4,facecolor='k')
		fig.show_arrows(277.5149,-2.0544,0.0005,0.0005,width=.15,head_width=.5,head_length=.4,facecolor='k')
		fig.show_arrows(277.5205,-2.0550,-0.00075,0.00075,width=.15,head_width=.5,head_length=.4,facecolor='k')
		fig.show_arrows(277.5115,-2.0585,0.00125,0.0005,width=.15,head_width=.5,head_length=.4,facecolor='k')

	# Add SPIRE 350 micron contour from figure 1
	add_spire = False
	if add_spire == True:
		fig.show_contour('FITS/SPIRE_350_micron.fits', alpha=.8, cmap='Greys_r', levels=[380], hdu=1)
	if zoom == True:
		fig.recenter(277.518,-2.043,width=0.045,height=0.045)
		fig.hide_layer('rectangle_set_1')
	elif vla4 == True:
		fig.recenter(277.512,-2.012,width=0.02,height=0.015)
		fig.hide_layer('label_6')
		# Add Scalebar
		fig.add_scalebar(10./3600)
		fig.scalebar.set_corner('bottom')
		fig.scalebar.set_frame(False)
		fig.scalebar.set(linewidth=2,color='black')
		fig.scalebar.set_label('0.02 pc')
		fig.scalebar.set_font(size=18)
	elif vla10 == True:
		fig.recenter(277.512,-2.043,width=0.01,height=0.01)
		# Add Scalebar
		fig.add_scalebar(10./3600)
		fig.scalebar.set_corner('bottom right')
		fig.scalebar.set_frame(False)
		fig.scalebar.set(linewidth=2,color='black')
		fig.scalebar.set_label('0.02 pc')
		fig.scalebar.set_font(size=18)
	elif vla12 == True:
		# Re-center
		fig.recenter(277.516,-2.054,width=0.015,height=0.0145)
		# Decrease label frequency
		labels = fig._ax1.get_xticklabels()
		labels = np.array(map(lambda x: x.get_text(),labels))
		labels[1::2] = u''
		fig._ax1.set_xticklabels(labels)
		# Add Scalebar
		fig.add_scalebar(10./3600)
		fig.scalebar.set_corner('bottom')
		fig.scalebar.set_frame(False)
		fig.scalebar.set(linewidth=2,color='black')
		fig.scalebar.set_label('0.02 pc')
		fig.scalebar.set_font(size=18)
		fig.refresh()
	else: 
		fig.recenter(277.525,-2.038,width=0.082,height=0.082)

	if savefig == True:
		fig.save(figname,adjust_bbox='tight',transparent=True)


### Radio Contours over Herschel Maps
if plot3 == True:

	zoom = False
	zoom2 = True
	center = False
	
	figname = 'AAS_2015/Spitzer80_Contour.eps'

	savefig = True

	# Choose Image
	color_im = 'FITS/2MASS_Kband.fits'
	color_im = 'FITS/Spitzer_8.0_micron_crop.fits'
	#color_im = 'FITS/Spitzer_24_micron.fits'
	#color_im = 'FITS/PACS_70_micron_cropped.fits'
	#color_im = 'FITS/SPIRE_250_micron.fits'
	#contour_im = 'FITS/SERPS_2.lower.briggs05.deeper.fits'
	contour_im = 'FITS/SERPS_2.upper.briggs05.deeper.fits'

	# Load Header and FITS Data for Contour
	contour_head = fits.open(contour_im)[0].header
	contour_data = fits.open(contour_im)[0].data
	
	contour_data = contour_data[0][0]

	# Get RMS
	rms = np.std(contour_data[~np.isnan(contour_data)])

	# Initialize Upper Subband Plot
	fig = aplpy.FITSFigure(color_im, subplot=(1,1,1), figsize=(10.5,9))

	fig.show_colorscale(cmap='YlGnBu_r',stretch='power',exponent=.70,pmin=0,pmax=99.7)
	fig.show_contour('FITS/SERPS_2.upper.briggs05.deeper.fits',cmap='bwr_r',alpha=1,levels=[rms*4,rms*5,rms*6,rms*8,rms*10,rms*15,rms*20,rms*30,rms*50],zorder=10,linewidth=4)

	# Draw Beam
	draw_beam = False
	if draw_beam == True:
		beam_dict = {'bmaj':contour_head['bmaj'],'bmin':contour_head['bmin'],'bpa':contour_head['bpa']}
		fig._header.update(beam_dict)
		fig.show_beam(facecolor='black',frame=True)

	# Miscellaneous
	fig.ticks.set_color('#000000')
	fig.ticks.set_linewidth(0.75)
	fig.ticks.set_length(8)
	fig.ticks.set_minor_frequency(5)
	fig.ticks.show_x()
	fig.ticks.show_y()
	fig.axis_labels.set_xtext('Right Ascension (J2000)')
	fig.axis_labels.set_ytext('Declination (J2000)')
	fig.axis_labels.set_xtext('')
	fig.axis_labels.set_ytext('')
	fig.axis_labels.set_font(size=16)
	fig.axis_labels.set_xpad(15)
	fig.axis_labels.set_ypad(10)
	fig.axis_labels.set_xpad(0)
	fig.axis_labels.set_ypad(0)
	fig._ax1.set_xticklabels([''])
	fig._ax1.set_yticklabels([''])
	fig._figure.axes[0].tick_params(labelsize=16)
	fig.refresh()

	fig.recenter(277.520,-2.048,width=0.055,height=0.055)
	#fig.recenter(277.53,-2.037,width=0.07,height=0.072)

	## Draw Regions ##

	# Add Gaussian Fits
	add_gaussfits = False
	if add_gaussfits == True:
		fig.show_regions('Gauss_Fitting/GaussFits_upper.reg',zorder=1)
		fig.hide_layer('region_set_1_txt')
		fig.hide_layer('region_set_1')

	# Add IRAM Millimeter sources
	add_iram = False
	if add_iram == True:
		iram_source_name, iram_ra, iram_dec, iram_Lbol, iram_class, iram_S, iram_S_err, iram_FWHM = np.loadtxt('IRAM_sources.txt',unpack=True,delimiter=',',dtype='str')
		iram_ra_deg, iram_dec_deg = sexig_to_deg(iram_ra, iram_dec)
		fig.show_markers(iram_ra_deg,iram_dec_deg,s=200, marker='x', c='black',zorder=1,alpha=.8,lw=1)

	# Add Spitzer Infrared Sources
	add_spitz = False
	if add_spitz == True:
		spitzer = np.loadtxt('gutermuth2008_serps_ysos1.txt',unpack=True,dtype='str')
		spitz_radec = np.array([np.array(spitzer[1],float),np.array(spitzer[2],float)])
		classI = spitzer[3]=='I'
		fig.show_markers(spitz_radec[0][classI],spitz_radec[1][classI],s=250,marker='+',facecolor='None',edgecolor='blue',lw=2,zorder=1)
		fig.show_markers(spitz_radec[0][~classI],spitz_radec[1][~classI],s=250,marker='+',facecolor='None',edgecolor='red',lw=2,zorder=1)

	# Add VLA Labels
	add_labels = False
	if add_labels == True:
		label_ra,label_dec,label_num = np.loadtxt('SS_Labels.txt',delimiter=',',unpack=True)
		label_num = np.array(label_num,int)
		for i in range(len(label_num)):
			fig.add_label(label_ra[i],label_dec[i],'VLA '+str(label_num[i]),clip_on=True,zorder=5)

	# Add Teixeira Sources
	add_teixeira = False
	if add_teixeira == True:
		fig.show_regions('Teixeira12_Sources.reg')	

	# Add Caption
	add_cap = True
	if add_cap == True:
		#fig.add_label(0.1,0.9,'(f)',size=24,zorder=5,color='k',relative=True)
		props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
		fig.add_label(.23,.92,'8.0 micron', size=30, zorder=5, color='k', relative=True, bbox=props)

	if zoom == True:
		fig.recenter(277.518,-2.043,width=0.0425,height=0.045)
	elif zoom2 == True:
		fig.recenter(277.516,-2.049,width=0.018,height=0.0225)
	elif center == True:
		fig.recenter(277.518,-2.046,width=0.04,height=0.042)
	else: 
		fig.recenter(277.525,-2.038,width=0.082,height=0.082)

	if savefig == True:
		fig.save(figname,adjust_bbox='tight',transparent=False)



## Spix vs. Class plot
if plot4 == True:
	import pickle as pkl

	# Load Radio Data
	data = pkl.Unpickler(open('Radio_Data.pkl','rb')).load()
	for i in data:
		try:
			data[i] = np.array(data[i],float)
		except ValueError:
			pass
	kern = AttrDict(data)
	kern.Lbol_sources = np.array(kern.Lbol_sources,int)


	# Plot
	class0,class1,class2,extragal = [],[],[],[]
	for i in range(len(kern.classes)):
		if 'Class 0' in kern.classes[i]: class0.append(i)
		if 'Class I' in kern.classes[i] and 'Class II' not in kern.classes[i]: class1.append(i)
		if 'Class II' in kern.classes[i]: class2.append(i)
		if 'Extragal' in kern.classes[i] and 'Extragal.?' not in kern.classes[i]: extragal.append(i)

	class0,class1,class2,extragal =  np.array(class0),np.array(class1),np.array(class2),np.array(extragal)

	trust_int_spix = np.where(kern.cm_spix_trust == 'int')[0]
	trust_peak_spix = np.where(kern.cm_spix_trust == 'peak')[0]

	fig,ax = mp.subplots()

	for i in range(18):
		if i in class0:
			if i in trust_int_spix:
				ax.errorbar([0],kern.cm_int_spix[i],yerr=kern.cm_int_spix_err[i],fmt='s',color='green',alpha=.7)
			elif i in trust_peak_spix:
				ax.errorbar([0],kern.cm_peak_spix[i],yerr=kern.cm_peak_spix_err[i],fmt='s',color='green',alpha=.7)
		elif i in class1:
			if i in trust_int_spix:
				ax.errorbar([1],kern.cm_int_spix[i],yerr=kern.cm_int_spix_err[i],fmt='s',color='blue',alpha=.7)
			elif i in trust_peak_spix:
				ax.errorbar([1],kern.cm_peak_spix[i],yerr=kern.cm_peak_spix_err[i],fmt='s',color='blue',alpha=.7)
		elif i in class2:
			if i in trust_int_spix:
				ax.errorbar([2],kern.cm_int_spix[i],yerr=kern.cm_int_spix_err[i],fmt='s',color='red',alpha=.7)
			elif i in trust_peak_spix:
				ax.errorbar([2],kern.cm_peak_spix[i],yerr=kern.cm_peak_spix_err[i],fmt='s',color='red',alpha=.7)
		elif i in extragal:
			if i in trust_int_spix:
				ax.errorbar([3],kern.cm_int_spix[i],yerr=kern.cm_int_spix_err[i],fmt='s',color='black',alpha=.7)
			elif i in trust_peak_spix:
				ax.errorbar([3],kern.cm_peak_spix[i],yerr=kern.cm_peak_spix_err[i],fmt='s',color='black',alpha=.7)

	ax.set_xlim(-1,4)
	ax.set_ylim(-4,4)
	mp.xticks([0,1,2,3],['Class 0','Class I','Class II','Extragal.'],fontsize=16)
	ax.set_ylabel('Radio Spectral Index',fontsize=16)

	if savefig == True:
		fig.savefig('Paper/figures/Spix_Class.eps',adjust_bbox='tight',transparent=True)
