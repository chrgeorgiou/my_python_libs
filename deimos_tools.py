import numpy as np
import os
import subprocess
from astropy.io import fits
from astropy.table import Table

def read_config(config_file):
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

def get_rband_segrad(ID, rband_cat, matched_cat=None, append=False, matchto='CATAID'):
	from astropy.io import fits
	rband = fits.getdata(rband_cat, 1)
	
	segrad = []
	for i in ID:
		idx = np.where(i == rband[matchto])[0]
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

def weight_function(Ngal, sizes, isorad=None):
	"""
	Construct a weight function array.
	:param Ngal: Total number of galaxies.
	:param sizes: Sizes list. E.g. [1,'pixels'] for one value in pixels, [2,5,10, 'isophote'] for three different values, corresponding to multiples of the isophotal radius.
	:param isorad: Array containing the isophotal radius.
	:return: Weight function in the shape of (#objects, #sizes). E.g. [ [1,2], [1,2], [1,2] ] for three galaxies with two weight functions, at 1 and 2 pixels.
	"""
	wfsize = []
	if sizes[-1] == 'pixels':
		if len(sizes) == 2:
			Nsizes = 1
			wfsize = np.full((Ngal, Nsizes), sizes[0], dtype=np.float)
		elif len(sizes) == 4:
			Nsizes = sizes[2]
			for i in range(Ngal): wfsize.append(np.linspace(sizes[0], sizes[1], sizes[2]))
			wfsize = np.array(wfsize)
	elif sizes[-1] == 'isophote':
		for i in range(Ngal): wfsize.append([isorad[i] * s for s in sizes[:-1]])
		wfsize = np.array(wfsize)
	elif sizes[-1] == 'kpc':
		exit('WORK IN PROGRESS, DON\'T USE kpc yet!')
	
	print 'Weight function used:\n', sizes[:-1], sizes[-1], '\n'
	return wfsize

def create_input_catalogues(ID, x, y, input_catalogue, sextractor_catalogue, post_stamp=None, post_stamp_max=300):
	if post_stamp is None: post_stamp = np.full_like(x, post_stamp_max)

	# Save a catalog with NUMBER, XWIN_IMAGE, YWIN_IMAGE, MAG_AUTO, FLUX_RADIUS, FLUX_AUTO, FLUXERR_AUTO, FLAGS, CLASS_STAR (sample file is 'sim_m6m3.cat')
	# Dummy variables that are not used
	mag = np.full_like(x, 19.0)  # kidscat['MAG_AUTO'][mask]
	frad = np.full_like(x, 3.0)  # kidscat['FLUX_RADIUS'][mask]
	f = np.full_like(x, 1000.0)  # kidscat['FLUX_AUTO'][mask]
	ferr = np.full_like(x, 10.0)  # kidscat['FLUXERR_AUTO'][mask]
	flags = np.full_like(ID, 0)  # kidscat['Flag'][mask]
	class_star = np.full_like(x, 0.3)  # kidscat['CLASS_STAR'][mask]

	sex_cat_array = np.column_stack((ID, x, y, mag, frad, f, ferr, flags, class_star))
	np.savetxt(sextractor_catalogue, sex_cat_array, fmt='%i %s %s %s %s %s %s %i %s')

	# Save a catalog with ID, x, y, x/y-min/max, and SEx flag (=0)
	in_cat = open(input_catalogue, 'w')
	in_cat.write(
		'# 1 ID\n# 2 XCENTROID_IMAGE\n# 3 YCENTROID_IMAGE\n# 4 XMIN_IMAGE\n# 5 XMAX_IMAGE\n# 6 YMIN_IMAGE\n# 7 YMAX_IMAGE\n# 8 FLAGS\n')
	for i in range(len(ID)): in_cat.write('%i %s %s %i %i %i %i %i\n' % (
	ID[i], x[i], y[i], x[i] - post_stamp[i] / 2 + 0.5,
	x[i] + post_stamp[i] / 2 + 0.5, y[i] - post_stamp[i] / 2 + 0.5,
	y[i] + post_stamp[i] / 2 + 0.5, 0))
	in_cat.close()

def create_object_catalogue_wcs(image, ra, dec, ID, output):
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

	t = Table([ID[mask], x[mask], y[mask], ra[mask], dec[mask]],
		names=('CATAID', 'Xpos_Awe', 'Ypos_Awe', 'RA_GAMA', 'DEC_GAMA'),
		dtype=('int32', 'float32', 'float32', 'float32', 'float32'), )
	t['Xpos_Awe'].unit, t['Ypos_Awe'].unit, t['RA_GAMA'].unit, t['DEC_GAMA'].unit = 'pixel', 'pixel', 'deg', 'deg'
	if not os.path.exists(os.path.dirname(output)): os.makedirs(os.path.dirname(output))
	t.write(output, format='fits', overwrite=True)
	return True

def read_object_catalogue(catalogue, sort=False):
	data = fits.getdata(catalogue, 1)
	ID = data['CATAID']
	if sort:
		ind = np.lexsort((ID, ID))
	else:
		ind = np.full(ID.shape, True, dtype=bool)
	ID = ID[ind]
	x = data['Xpos_Awe'][ind]
	y = data['Ypos_Awe'][ind]
	try:
		isorad = data['IsoRadius'][ind]
		remove0 = isorad > 0.  # filter out problematic measurements
		ID, x, y, isorad = ID[remove0], x[remove0], y[remove0], isorad[remove0]
	except:
		isorad = None
	return ID, x, y, isorad

class deimos_tools(object):
	
	def __init__(self, *args, **kwargs):
		"""

		:param args:
		:param kwargs:
		"""
		self.image = args[0]
		self.segmap = args[1]
		self.psf_catalogue = args[2]
		self.sextractor_catalogue = args[3]
		self.input_catalogue = args[4]
		self.noise = args[5]
		self.output_catalogue = args[6]
		self.sizes = args[7]
		self.correction_order = args[8]
		self.deimos_script = args[9]
		self.run_type='single'
		self.Ncores = args[10]

	def run_deimos_script(self, N_sizes, C, output_cat, psf_cat, sextractor_cat, input_cat, current_size=None, size_unit=None):
		# Set the desired output filename
		if current_size!=None:
			output_cat += '_C_%s_%s_%s.dat' %(C, current_size, size_unit)

		bashCommand = '%s -n %s -f %s -o %s -s %s -C %i -N 2 -p %s -x %s -c %s -m -e -T %s' % (self.deimos_script, self.noise, self.image, output_cat, N_sizes, C, psf_cat, sextractor_cat, input_cat, self.Ncores)
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

		try: len(self.correction_order)
		except: c = [self.correction_order]
		print 'Correction order used:\n', self.correction_order, '\n'
		for C in self.correction_order:
			if (self.run_type=='multisize') or (self.run_type=='mixed'): self.run_deimos_script(len(self.Nsizes), C, self.output_catalogue, self.psf_catalogue, self.sextractor_catalogue, self.input_catalogue)
			if (self.run_type=='mixed') or (self.run_type=='single'):
				for idx in range(len(self.sizes[:-1])):
					psf_cat_onesize = open(os.path.basename(self.psf_catalogue)+'.temp', 'w')
					line = np.loadtxt(self.psf_catalogue, unpack=True, dtype=str)
					# workaround for the case when "line" is just a single line
					if type(line[0])!=np.ndarray:
						temp = []
						temp += [[line[i]] for i in range(len(line))]
						line = temp
					for i in range(len(line[0])):
						psf_cat_onesize.write('%s %s %s %s %s %s %s %s %s %s\n' % (line[0][i], line[1][i], line[2][i], line[3+7*idx][i], line[4+7*idx][i], line[5+7*idx][i], line[6+7*idx][i], line[7+7*idx][i], line[8+7*idx][i], line[9+7*idx][i]) )
					psf_cat_onesize.close()
					self.run_deimos_script(1, C, self.output_catalogue, os.path.basename(self.psf_catalogue)+'.temp', self.sextractor_catalogue, self.input_catalogue, current_size=self.sizes[idx], size_unit=self.sizes[-1])
					try: os.remove(os.path.basename(self.psf_catalogue)+'.temp')
					except OSError: pass
