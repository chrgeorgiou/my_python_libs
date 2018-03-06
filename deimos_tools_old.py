import numpy as np
import os
import subprocess
from astropy.io import fits
from astropy.table import Table

def read_config(config_file):
	#def read_configuration_file(config_file='deimos_configuration.txt'):
	#config_file = 'deimos_configuration.txt'
	config = open(config_file)
	outparams = dict()
	for lines in config:
		if (lines.startswith('#'))or(lines=='\n'): continue
		line = lines.split()
		if line[0] == 'image': outparams[line[0]] = line[1]
		if line[0] == 'image_folder': outparams[line[0]] = line[1]
		if line[0] == 'band': outparams[line[0]] = line[1]
		if line[0] == 'segmap_folder': outparams[line[0]] = line[1]
		if line[0] == 'rband_catalogue': outparams[line[0]] = line[1]
		if line[0] == 'sextractor_config': outparams[line[0]] = line[1]
		if line[0] == 'GAMA_catalogue': outparams[line[0]] = line[1]
		if line[0] == 'matched_catalogue_directory': outparams[line[0]] = line[1]
		if line[0] == 'post_stamp_max': outparams[line[0]] = int(line[1])
		if line[0] == 'post_stamp_min': outparams[line[0]] = int(line[1])
		if line[0] == 'psfmap_folder': outparams[line[0]] = line[1]
		if line[0] == 'sizes': outparams[line[0]] = [float(i) for i in line[1:len(line)-1]]+[line[-1]]
		if line[0] == 'correction_order': outparams[line[0]] = [int(i) for i in line[1:len(line)]]
		if line[0] == 'output_directory': outparams[line[0]] = line[1]
		if line[0] == 'deimos_script': outparams[line[0]] = line[1]
		if line[0] == 'run_type': outparams[line[0]] = line[1]
		if line[0] == 'Ncores': outparams[line[0]] = int(line[1])
	return outparams

def get_rband_segrad(ID, rband_cat, matched_cat=None, append=False):
	from astropy.io import fits
	rband = fits.getdata(rband_cat, 1)
	
	segrad = []
	for i in ID:
		idx = np.where(i == rband['CATAID'])[0]
		if len(idx) == 0:
			segrad.append(0.0)
		else:
			segrad.append(rband['IsoRadius'][idx[0]])
	segrad = np.array(segrad)
	
	if (append) and (matched_cat != None):
		from astropy.table import Table
		t = Table.read(matched_cat)
		t['IsoRadius'] = segrad
		t['IsoRadius'].unit = 'pixel'
		t.write(matched_cat, format='fits', overwrite=True)
	return segrad


def calculate_postage_stamps(X, Y, sextractor_out, noise, min_postage_stamp, max_postage_stamp):
	from scipy import optimize
	from astropy.io import fits
	sextr_data = fits.getdata(sextractor_out, 2)
	id_sex, x_sex, y_sex = sextr_data['NUMBER'], sextr_data['X_IMAGE'], sextr_data['Y_IMAGE']
	from catalogues import spatial_matching
	index_list = spatial_matching(x_sex, y_sex, X, Y, 10)
	matched_obj = index_list[:, 0] != len(x_sex)
	sex_ids = index_list[matched_obj][:,0]
	cat_ids = index_list[matched_obj][:,1]
	q, R, I, iso = 1./sextr_data['ELONGATION'][sex_ids], sextr_data['FLUX_RADIUS'][sex_ids], sextr_data['FLUX_AUTO'][sex_ids], np.sqrt(sextr_data['ISOAREA_IMAGE'][sex_ids]/np.pi)
	pds = np.full(len(X), max_postage_stamp, dtype=np.float)
	for i in range(len(q)):
		xarray = np.linspace(10, 1000,100)
		maxarg = np.argmax(Q(xarray, 6, 6, q[i], R[i], I[i], iso[i], noise))
		try: x0 = xarray[maxarg+1]
		except: x0 = xarray[maxarg] # For when the maximum value happens at x=1000
		try: res = np.abs( optimize.newton(Q, x0, args=(6, 6, q[i], R[i], I[i], iso[i], noise)) )
		except RuntimeError: res = max_postage_stamp/2 # If the optimization fails just use max poststamp
		if res > max_postage_stamp/2: pds[cat_ids[i]] = max_postage_stamp
		elif res < min_postage_stamp/2: pds[cat_ids[i]] = min_postage_stamp
		else: pds[cat_ids[i]] = 2*res
	return pds

def Q(x, i, n, q, R, F, s, noise=0.):
	from scipy.special import gamma
	bn = 1.9992*n-0.3271
	Ie = F * q * bn**(2*n) * np.exp(-bn) / ( 2 * np.pi * n * R**2 * gamma(2*n) )
	W = np.exp( - 2 * q**2. * x**2. / ( s**2. * (1+q)**2. ) )
	I = Ie * np.exp( - bn * ( (np.abs(x) * np.sqrt(q) / R)**(1./n) -1 ) )
	return I*W*x**i-noise

def Qprime(x, i, n, q, R, F, s, noise=0.):
	bn = 1.9992*n-0.3271
	return Q(x, i, n, q, R, F, s) * ( - bn * (np.sqrt(q)/R)**(1./n) * x**(1./n-1.) / n - 4 * q**2. * x / ( s**2. * (1+q)**2. ) + 1.*i / x )

class deimos_tools(object):
	
	def __init__(self, params):
		self.image = params['image']
		imagename = os.path.basename(self.image)
		self.segmap = os.path.join( params['segmap_folder'], imagename.replace('.fits', '_segmap.fits') )
#		RA = image.split('/')[-1].split('_')[1].replace('.', 'p')
#		DEC = image.split('/')[-1].split('_')[2].replace('.', 'p').replace('-', 'm')
		self.psfmap = os.path.join( params['psfmap_folder'], imagename+'.psf.map' )
		self.psf_catalogue = os.path.join( params['output_directory'], 'input', imagename + '_psf_moments.cat' )
		self.input_catalogue = os.path.join( params['output_directory'], 'input', imagename + '_input.cat' )
		self.sextractor_catalogue = os.path.join( params['output_directory'], 'input', imagename + '_sextractor.cat' )
		self.matched_catalogue = os.path.join( params['matched_catalogue_directory'], imagename.replace('fits', 'GAMA_matched_wcs.cat') )
		self.noise_level = os.path.join( params['output_directory'], 'input', imagename.replace('fits', 'noise_level.txt') )
		if not os.path.exists(os.path.dirname(self.input_catalogue)): os.makedirs(os.path.dirname(self.input_catalogue))

		self.output_catalogue = os.path.join( params['output_directory'], imagename + 'DEIMOS_shapes.cat')
		self.sizes = params['sizes']
		self.deimos_script = params['deimos_script']
		self.post_stamp_max = params['post_stamp_max']
		self.post_stamp_min = params['post_stamp_min']
		self.correction_order = params['correction_order']
		self.run_type = params['run_type']
		
	def create_object_catalogue_wcs(self, image, ra, dec, ID):
		'''
		Assumes image[y,x], standard for fits files.
		'''
		from astropy.wcs import wcs
		hdulist = fits.open(image)
		w = wcs.WCS(hdulist[0].header, hdulist)  # Read wcs from header
		x, y = w.all_world2pix(ra, dec, 1)  # find x,y of ra,dec in the w table
		ymax, xmax = hdulist[0].data.shape  # get image max y and x
		mask = (x > 0.0) & (x < xmax) & (y > 0.0) & (y < ymax)  # get the x,y that reside within the image
		if len(x[mask]) == 0:  # if no objects found, set nomatch to TRUE
			return False
		
		t = Table(
			[ID[mask], x[mask], y[mask], ra[mask], dec[mask]],
			names = ('CATAID', 'Xpos_Awe', 'Ypos_Awe', 'RA_GAMA', 'DEC_GAMA'),
			dtype=('int32', 'float32', 'float32', 'float32', 'float32'),
			)
		t['Xpos_Awe'].unit, t['Ypos_Awe'].unit, t['RA_GAMA'].unit, t['DEC_GAMA'].unit = 'pixel', 'pixel', 'deg', 'deg'
		t.write(self.matched_catalogue, format='fits', overwrite=True)
		return True
		
	def read_object_catalogue(self, sort=False):
		if not os.path.isfile(self.matched_catalogue):
			print 'Object catalogue missing.\n%s\n' % self.matched_catalogue
			return False
		data = fits.getdata(self.matched_catalogue, 1)
		self.ID = data['CATAID']
		if sort: ind = np.lexsort((self.ID, self.ID))
		else: ind = np.full(self.ID.shape, True, dtype=bool)
		self.ID = self.ID[ind]
		self.x = data['Xpos_Awe'][ind]
		self.y = data['Ypos_Awe'][ind]
		try:
			self.isorad = data['IsoRadius'][ind]
			remove0 = self.isorad > 0.  # filter out problematic measurements
			self.ID, self.x, self.y, self.isorad = self.ID[remove0], self.x[remove0], self.y[remove0], self.isorad[remove0]
		except:
			pass
		return True
		
	def weight_function(self):
		Ngal = len(self.ID)
		wfsize = []
		if self.sizes[-1] == 'pixels':
			if len(self.sizes) == 2:
				Nsizes = 1
				wfsize = np.full((Ngal, Nsizes), self.sizes[0], dtype=np.float)
			elif len(self.sizes) == 4:
				Nsizes = self.sizes[2]
				for i in range(Ngal): wfsize.append(np.linspace(self.sizes[0], self.sizes[1], self.sizes[2]))
				wfsize = np.array(wfsize)
		elif self.sizes[-1] == 'isophote':
			for i in range(Ngal): wfsize.append([self.isorad[i] * s for s in self.sizes[:-1]])
			wfsize = np.array(wfsize)
		elif self.sizes[-1] == 'kpc':
			exit('WORK IN PROGRESS, DON\'T USE kpc yet!')
		self.wfsize = wfsize
		
	def create_input_catalogues(self):
		try: self.post_stamp
		except: self.post_stamp = self.post_stamp_max

		# Save a catalog with NUMBER, XWIN_IMAGE, YWIN_IMAGE, MAG_AUTO, FLUX_RADIUS, FLUX_AUTO, FLUXERR_AUTO, FLAGS, CLASS_STAR (sample file is 'sim_m6m3.cat')
		# Dummy variables that are not used
		mag = np.full_like(self.x, 19.0)# kidscat['MAG_AUTO'][mask]
		frad = np.full_like(self.x, 3.0)# kidscat['FLUX_RADIUS'][mask]
		f = np.full_like(self.x, 1000.0)# kidscat['FLUX_AUTO'][mask]
		ferr = np.full_like(self.x, 10.0)# kidscat['FLUXERR_AUTO'][mask]
		flags = np.full_like(self.ID, 0)# kidscat['Flag'][mask]
		class_star = np.full_like(self.x, 0.3)# kidscat['CLASS_STAR'][mask]

		sex_cat_array = np.column_stack((self.ID,self.x,self.y,mag,frad,f,ferr,flags,class_star))
		np.savetxt(self.sextractor_catalogue, sex_cat_array, fmt='%i %s %s %s %s %s %s %i %s')
		
		# Save a catalog with ID, x, y, x/y-min/max, and SEx flag (=0)
		in_cat = open(self.input_catalogue, 'w')
		in_cat.write('# 1 ID\n# 2 XCENTROID_IMAGE\n# 3 YCENTROID_IMAGE\n# 4 XMIN_IMAGE\n# 5 XMAX_IMAGE\n# 6 YMIN_IMAGE\n# 7 YMAX_IMAGE\n# 8 FLAGS\n')
		for i in range(len(self.ID)): in_cat.write('%i %s %s %i %i %i %i %i\n' % (self.ID[i], self.x[i], self.y[i], self.x[i] - self.post_stamp[i] / 2 + 0.5, self.x[i] + self.post_stamp[i] / 2 + 0.5, self.y[i] - self.post_stamp[i] / 2 + 0.5, self.y[i] + self.post_stamp[i] / 2 + 0.5, 0))
		in_cat.close()
		
	def run_deimos_script(self, N_sizes, C, output_cat, psf_cat, sextractor_cat, input_cat, current_size=None):
		# Set the desired output filename
		if N_sizes > 1: output_cat = output_cat.replace('.cat', '_C_%s_multisize.cat' % C)
		else: output_cat = output_cat.replace('.cat', '_C_%s_size_%s.cat' % (C, current_size))

		bashCommand = '%s -n %s -f %s -o %s -s %s -C %i -N 2 -p %s -x %s -c %s -m -e' % (self.deimos_script, self.noise, self.image, output_cat, N_sizes, C, psf_cat, sextractor_cat, input_cat)
		if os.path.isfile(self.segmap): bashCommand += ' -S %s' % self.segmap
		#print bashCommand
		process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE, stdin=subprocess.PIPE)
		stdout, stderr = process.communicate()
		print stdout

	def deimos_shape_measurement(self):
		'''
		Run the MomentsGalPAR script for the different correction orders and sizes provided.
		'''
		if not os.path.exists( os.path.dirname(self.output_catalogue) ): os.makedirs( os.path.dirname(self.output_catalogue) )

		for C in self.correction_order:
			if (self.run_type=='multisize') or (self.run_type=='mixed'): self.run_deimos_script(len(self.sizes[:-1]), C, self.output_catalogue, self.psf_catalogue, self.sextractor_catalogue, self.input_catalogue)
			if (self.run_type=='mixed') or (self.run_type=='single'):
				for idx, size in enumerate(self.sizes[:-1]):
					psf_cat_onesize = open(os.path.basename(self.image)+'.psf_catalogue.cat', 'w')
					line = np.loadtxt(self.psf_catalogue, unpack=True, dtype=str)
					# workaround for the case when "line" is just a single line
					if type(line[0])!=np.ndarray:
						temp = []
						temp += [[line[i]] for i in range(len(line))]
						line = temp
					for i in range(len(line[0])):
						psf_cat_onesize.write('%s %s %s %s %s %s %s %s %s %s\n' % (line[0][i], line[1][i], line[2][i], line[3+7*idx][i], line[4+7*idx][i], line[5+7*idx][i], line[6+7*idx][i], line[7+7*idx][i], line[8+7*idx][i], line[9+7*idx][i]) )
					psf_cat_onesize.close()
					self.run_deimos_script(1, C, self.output_catalogue, os.path.basename(self.image)+'.psf_catalogue.cat', self.sextractor_catalogue, self.input_catalogue, current_size=size)
					try: os.remove(os.path.basename(self.image)+'.psf_catalogue.cat')
					except OSError: pass
