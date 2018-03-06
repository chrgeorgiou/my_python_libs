import numpy as np
from astropy.io import fits
import os
import subprocess
import shutil
import sys

class segmentation_map(object):
	
	def __init__(self, image, segmap, sextractor_conf, x, y, ID):
		'''
		Class object that uses Source Extractor to create a segmentation map of an image and uses it to measure the isophotal radius
		:param image: The image on which Source Extractor will be run
		:param segmap: The location of the output segmentation map
		:param sextractor_conf: The configuration file that will be used
		:param x, y: The positions of the objects whose isophotal radius will be measured
		:param ID: The ID of the objects. Can turn the Source Extractor ID's to the ID's provided here.
		'''
		self.image = image
		self.segmap = segmap
		self.sextractor_conf = sextractor_conf
		checkimage = os.path.basename(self.segmap)
		catalog = checkimage.replace('.fits', '.cat')
		self.sextractor_catalogue = os.path.join(os.path.dirname(self.segmap), 'catalogs', catalog)
		self.x = x
		self.y = y
		self.ID = ID
	
	def create_segmentation_map(self, thresh, quiet=False, save_map=True, save_cat=True, append=False, catalogue=None, correct_ids=False):
		# thresh[0] is analysis and thresh[1] is detection threshold
		try: len(thresh)
		except: thresh = (thresh, thresh)
		
		# Move to the folder where the map will be stored, to avoid copying delays
		cwd = os.getcwd()
		os.chdir(os.path.dirname(self.segmap))
		
		# Check if the default files are present and if not, copy them from the class parent folder
		realpath = os.path.realpath(__file__)
		sex_def_dir = os.path.join( os.path.dirname(realpath), 'config')
		if not os.path.isfile('default.conv'): shutil.copy( os.path.join(sex_def_dir, 'default.conv'), os.getcwd())
		if not os.path.isfile('default.nnw'): shutil.copy( os.path.join(sex_def_dir, 'default.nnw'), os.getcwd())
		if not os.path.isfile('default.param'):	shutil.copy( os.path.join(sex_def_dir, 'default.param'), os.getcwd())
		if not os.path.isfile(os.path.basename(self.sextractor_conf)): shutil.copy(self.sextractor_conf, os.getcwd())

		checkimage = os.path.basename(self.segmap)
		catalog = checkimage.replace('.fits', '.cat')
		
		# Run SExtractor
		bashCommand = 'sex %s -c %s  -ANALYSIS_THRESH %s -DETECT_THRESH %s -CHECKIMAGE_NAME %s -CHECKIMAGE_TYPE SEGMENTATION -CATALOG_NAME %s' % (self.image, self.sextractor_conf, thresh[0], thresh[1], checkimage, catalog)
		if quiet:
			bashCommand += ' -VERBOSE_TYPE QUIET'
			process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE, stdin=subprocess.PIPE,
			                           stderr=open(os.devnull, 'w'))
		else: process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE, stdin=subprocess.PIPE)
		process.wait()

		if not save_map: os.remove(checkimage)
		if save_cat:
			if not os.path.exists(os.path.join( os.path.dirname(self.segmap),'catalogs' )): os.makedirs(os.path.join( os.path.dirname(self.segmap),'catalogs' ))
			shutil.move(catalog, self.sextractor_catalogue)
		else: os.remove(catalog)
	
		# Set the IDs to the desired value and retrieve the radius of the isophotal area
		if save_cat: self.segrad = self.get_ISOAREA_radius(catalog, append, catalogue)
		if correct_ids: self.correct_segmap_ids(segmap=checkimage)
		os.chdir(cwd)
	
	def correct_segmap_ids(self, segmap='self'):
		# Swap the Source Extractor Seq Numbers with the ID numbers provided (only on the segmentation map, not on the Source Extractor catalogue).
		if segmap == 'self': segmap = self.segmap
		seghdu = fits.open(segmap, ignore_missing_end=True) # Fix for when the END header keyword is somehow not read properly
		segdata = seghdu[0].data
		ymax, xmax = segdata.shape
		xint = np.array(self.x + 0.5, dtype=int)
		yint = np.array(self.y + 0.5, dtype=int)
		for i in range(len(xint)):
			# Skip cases where object's centre is exactly on the edge of the image
			if (xint[i] == xmax) or (yint[i] == ymax):
				continue
			seqnr = segdata[yint[i], xint[i]]
			seg_ps_dim = 150 # This should be large enough to encompass every object, but not too large to avoid running for too long
			seg_ps = segdata[yint[i] - seg_ps_dim:yint[i] + seg_ps_dim, xint[i] - seg_ps_dim:xint[i] + seg_ps_dim]
			mask = seg_ps == seqnr
			seg_ps[mask] = np.full_like(seg_ps[mask], self.ID[i])
			segdata[yint[i] - seg_ps_dim:yint[i] + seg_ps_dim, xint[i] - seg_ps_dim:xint[i] + seg_ps_dim] = seg_ps
		
		seghdu[0].data = segdata
		seghdu.writeto(segmap, clobber=True)
		
	def get_ISOAREA_radius(self, sex_catalogue, append=False, catalogue=None):
		from catalogues import spatial_matching # to match self.x, self.y with SExtractor's detected X, Y
		#if not os.path.isfile('catalogs/'+sex_catalogue): exit('SExtractor catalogue %s missing' %'catalogs/'+sex_catalogue)
		sex_data = fits.getdata('catalogs/'+sex_catalogue, 2)
		id_sex, x_sex, y_sex = sex_data['NUMBER'], sex_data['X_IMAGE'], sex_data['Y_IMAGE']
		
		index_list = spatial_matching(x_sex, y_sex, self.x, self.y, 5) # 5 pixels ~ 1 arcsec tolerance between GAMA and SExtractor centroid
		matched_obj = index_list[:, 0] != len(x_sex)
		sex_ids = index_list[matched_obj][:,0] # Matched object id's for SExtractor
		cat_ids = index_list[matched_obj][:,1] # Matched object id's for input X, Y, ID
		#sex_ids, cat_ids, dist = spatial_matching(x_sex, y_sex, self.x, self.y, 5)
		segrad = np.full(self.x.shape, 0.0) # Fill segrad array with zeros for all objects. Non-matched will be left as a zero
		segrad[cat_ids] = np.sqrt( sex_data['ISOAREA_IMAGE'][sex_ids]/np.pi ) # Fill the matched objects with the SExtractor ISOAREA values
		if append: self.save_segrad(catalogue, segrad)
		return segrad

	def save_segrad(self, catalogue, segrad):
		from astropy.table import Table
		t = Table.read(catalogue)
		t['IsoRadius'] = segrad
		t['IsoRadius'].unit = 'pixel'
		t.write(catalogue, format='fits', overwrite=True)

	def noise_calculation(self, plot=False):
		'''
		A simple noise calculation routine that is based on iteratively calculating the median of the image (subtracting the segmentation map). Fits a Gaussian to the flux curve of the image around the median of the image.
		:param plot: If True, a plot of the image flux around the noise level and noise standard deviation will be saved in the cwd
		'''
		try: segdata = fits.getdata(self.segmap)
		except: segdata = 0.
		imagedata = fits.getdata(self.image)
		mask = (segdata == 0.) & (imagedata != 0.)
		maskedimage = imagedata[mask]
		
		median = np.median(maskedimage)  # median of the whole data
		
		mask_low = maskedimage.flatten() < median
		median_low = np.median(
			maskedimage.flatten()[mask_low])  # find the mediam for the left tail of the noise (where no sources exist)
		mask = (maskedimage.flatten() > median_low) & (maskedimage.flatten() < 2. * median - median_low)
		median_final = np.median(
			maskedimage.flatten()[mask])  # find the median of the 'noise area' of the flux histogram
		
		mask = (maskedimage > 3 * median_low) & (maskedimage < 2 * median - 3 * median_low)
		self.noise = np.std(maskedimage[mask])
		if plot:
			import matplotlib.pyplot as plt
			import matplotlib.mlab as mlab
			n, bins, patches = plt.hist(maskedimage[mask], bins=300, log=True)
			plt.xlabel('Image Flux')
			plt.ylabel('#')
			y = mlab.normpdf(bins, median_final, self.noise)
			plt.plot(bins, len(maskedimage[mask]) * y / np.sum(y), 'r--')
			plt.savefig('flux_hist.pdf')
