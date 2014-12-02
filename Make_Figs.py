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

## Flags
savefig = False

plot1 = False
plot2 = False
plot3 = True


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
	# Load and Configure Images
	file1 = fits.open("FITS/SPIRE_250_micron.fits")
	file2 = fits.open("FITS/SPIRE_350_micron.fits")
	file3 = fits.open("FITS/SPIRE_500_micron.fits")

	data1 = file1[1].data
	wcs1 = WCS(file1[1].header)

	data2 = file2[1].data
	wcs2 = WCS(file2[1].header)

	data3 = file3[1].data
	wcs3 = WCS(file3[1].header)

	# Load Regions
	spitzer = np.loadtxt('gutermuth2008_serps_ysos1.txt',unpack=True,dtype='str')
	spitz_radec = np.array([np.array(spitzer[1],float),np.array(spitzer[2],float)])

	# Plot Image
	fig = mp.figure()
	ax = WCSAxes(fig,[0.1,0.1,0.8,0.8],wcs=wcs2)
	fig.add_axes(ax)

	ax.imshow(np.log(data2), origin='lower', cmap='Greys', interpolation='nearest',
		alpha=1.0)

	ax.contour(np.log(data2), levels=[5.9], colors='black', alpha=.8)

	#	ax.imshow(data2, origin='lower', cmap='Greens', extent=extent2,
	#		vmax=np.max(data2[~np.isnan(data2)])/4, alpha=.7 )
	#	ax.imshow(data1, origin='lower', cmap='Greys', extent=extent1, alpha=1,
	#		vmin=np.min(data1[~np.isnan(data1)])*.001, vmax=np.max(data1[~np.isnan(data1)])/3)

	## Plot Regions
	# Spitzer Protostars
	classI = spitzer[3]=='I'
	ax.plot(spitz_radec[0][classI],spitz_radec[1][classI],'x',color='#FF3333',markersize=8,transform=ax.get_transform('fk5'))
	ax.plot(spitz_radec[0][~classI],spitz_radec[1][~classI],'x',color='#3333FF',markersize=8,transform=ax.get_transform('fk5'))

	# VLA Field of View
	r = mp.Rectangle((277.485,-2.077),0.067,0.067,edgecolor='Lime',facecolor='None',lw=3,transform=ax.get_transform('fk5'))
	ax.add_patch(r)

	# Figure Formatting	limits = wcs2.wcs_world2pix([(277.35,-2.25),(277.68,-1.88)],0).T
	ax.set_xlim(limits[0][1],limits[0][0])
	ax.set_ylim(limits[1][0],limits[1][1])

	xax = ax.coords[0]
	yax = ax.coords[1]

	xax.set_axislabel('Right Ascension (J2000)',fontsize=13)
	yax.set_axislabel('Declination (J2000)',fontsize=13)

	xax.set_major_formatter('hh:mm:ss')
	yax.set_major_formatter('dd:mm')

	xax.set_ticks(spacing=5*u.arcmin, exclude_overlapping=True, size=6)
	yax.set_ticks(spacing=5*u.arcmin, exclude_overlapping=True, size=6)

	if savefig == True:
		fig.savefig('Paper/figures/SerpSouth_FOV.eps',bbox_inches='tight')


### Radio Contours
if plot2 == True:
	# Choose Image
	#im_name = 'FITS/SERPS_2.lower.mask.deeper2.briggs.fits'
	#im_name = 'FITS/SERPS_2.lower.briggs.source_clean.image.fits'
	#im_name = 'FITS/SERPS_2.upper.mask.deeper2.briggs.fits'
	#im_name = 'FITS/SERPS_2.upper.briggs.source_clean.image.fits'
	im_name = 'FITS/SERPS_2.lower.briggs05.deeper.fits'
	im_name = 'FITS/SERPS_2.upper.briggs05.deeper.fits'

	# Load Header and FITS Data
	header = fits.open(im_name)[0].header
	data = fits.open(im_name)[0].data

	data = data[0][0]

	# Get RMS
	rms = np.std(data[~np.isnan(data)])

	# Initialize Upper Subband Plot
	fig = aplpy.FITSFigure(im_name, subplot=(1,1,1))

	# Draw Contour
	fig.show_colorscale(cmap='gist_earth_r',stretch='power',exponent=.25,pmin=0,pmax=99.992)
	fig.hide_colorscale()
	fig.show_contour(colors='black',alpha=0.8,levels=[rms*3,rms*4,rms*5,rms*6,rms*8,rms*10,rms*15,rms*20,rms*30,rms*50])

	# Draw Beam
	fig.show_beam(facecolor='black',frame=True)

	# Draw Regions
	#fig.show_regions('SS_Targets_pxl.reg')

	# Miscellaneous
	fig.ticks.set_color('#000000')
	fig.ticks.set_linewidth(1)
	fig.ticks.set_length(8)
	fig.ticks.set_minor_frequency(1)
	fig.ticks.show_x()
	fig.ticks.show_y()
	fig.axis_labels.set_xtext('Right Ascension (J2000)')
	fig.axis_labels.set_ytext('Declination (J2000)')
	fig.axis_labels.set_font(size=14)

	fig.recenter(277.528,-2.038,width=0.075,height=0.075)

	if savefig == True:
		fig.save('Paper/figures/SerpSouth_upper_contourplot.eps',adjust_bbox='tight',dpi=300)


### Radio Contours over Herschel Maps
if plot3 == True:
	# Choose Image
	color_im = 'FITS/SPIRE_250_micron.fits'
	contour_im = 'FITS/SERPS_2.upper.briggs05.deeper.fits'

	# Load Header and FITS Data for Contour
	contour_head = fits.open(contour_im)[0].header
	contour_data = fits.open(contour_im)[0].data
	
	contour_data = contour_data[0][0]

	# Get RMS
	rms = np.std(contour_data[~np.isnan(contour_data)])

	# Initialize Upper Subband Plot
	fig = aplpy.FITSFigure(color_im, subplot=(1,1,1), figsize=(10.5,9))


	fig.show_colorscale(cmap='CMRmap',stretch='power',exponent=.30,pmin=0,pmax=100.005)
	fig.show_contour('FITS/SERPS_2.upper.briggs05.deeper.fits',alpha=1,levels=[rms*4,rms*5,rms*6,rms*8,rms*10,rms*15,rms*20,rms*30,rms*50])

	# Draw Beam
	beam_dict = {'bmaj':contour_head['bmaj'],'bmin':contour_head['bmin'],'bpa':contour_head['bpa']}
	fig._header.update(beam_dict)
	fig.show_beam(facecolor='black',frame=True)

	# Miscellaneous
	fig.ticks.set_color('#000000')
	fig.ticks.set_linewidth(1)
	fig.ticks.set_length(8)
	fig.ticks.set_minor_frequency(1)
	fig.ticks.show_x()
	fig.ticks.show_y()
	fig.axis_labels.set_xtext('Right Ascension (J2000)')
	fig.axis_labels.set_ytext('Declination (J2000)')
	fig.axis_labels.set_font(size=16)
	fig.axis_labels.set_xpad(15)
	fig.axis_labels.set_ypad(10)
	fig._figure.axes[0].tick_params(labelsize=14)
	fig.refresh()

	fig.recenter(277.520,-2.048,width=0.055,height=0.055)

	## Draw Regions ##

	# Add Gaussian Fits
	add_gaussfits = False
	if add_gaussfits == True:
		fig.show_regions('Gauss_Fitting/GaussFits_lower.reg',zorder=1)
		fig.hide_layer('region_set_1_txt')

	# Add IRAM Millimeter sources
	add_iram = True
	if add_iram == True:
		iram_source_name, iram_ra, iram_dec, iram_Lbol, iram_class, iram_S, iram_S_err, iram_FWHM = np.loadtxt('IRAM_sources.txt',unpack=True,delimiter=',',dtype='str')
		iram_ra_deg, iram_dec_deg = sexig_to_deg(iram_ra, iram_dec)
		fig.show_markers(iram_ra_deg,iram_dec_deg,s=200, marker='x', c='black',zorder=1,alpha=.8,lw=1)

	# Add Spitzer Infrared Sources
	add_spitz = True
	if add_spitz == True:
		spitzer = np.loadtxt('gutermuth2008_serps_ysos1.txt',unpack=True,dtype='str')
		spitz_radec = np.array([np.array(spitzer[1],float),np.array(spitzer[2],float)])
		classI = spitzer[3]=='I'
		fig.show_markers(spitz_radec[0][classI],spitz_radec[1][classI],s=225,marker='^',facecolor='None',edgecolor='Turquoise',lw=2,zorder=1)
		fig.show_markers(spitz_radec[0][~classI],spitz_radec[1][~classI],s=225,marker='^',facecolor='None',edgecolor='YellowGreen',lw=2,zorder=1)

	# Add VLA Labels
	add_labels = True
	if add_labels == True:
		label_ra,label_dec,label_num = np.loadtxt('SS_Labels.txt',delimiter=',',unpack=True)
		label_num = np.array(label_num,int)
		for i in range(len(label_num)):
			fig.add_label(label_ra[i],label_dec[i],'VLA '+str(label_num[i]),clip_on=True,zorder=5)


	if savefig == True:
		fig.save('Paper/figures/SerpSouth_uppercont_herschelcolor.eps')#,adjust_bbox='tight')












	
