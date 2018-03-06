import os
import subprocess
import shutil
import numpy as np
import sys
from sys import stdout
import math
import matplotlib.pyplot as plt
pi = np.pi

class psf_pipeline(object):

	def __init__(self, image, wfsizes, ID, X, Y, star_ID=None, star_x=None, star_y=None, output_dir = './'):
		"""
		Pipeline to calculate shapelets-based PSF moments at galaxy positions.

		:param image: image file to be processed
		:param wfsize: Array of shape ( len(ID), N ) where N is the number of different weight function sizes to be computed for every object
		:param ID: ID number(s) of the points where the psf moments are calculated
		:param X: x-axis position of the above object(s)
		:param Y: y-axis position of the above object(s)
		:param star_ID: ID of the star to be used for the shapelet-based PSF model (array-like).
		:param star_X: X-position of the star to be used for the shapelet-based PSF model (array-like).
		:param star_Y: Y-position of the star to be used for the shapelet-based PSF model (array-like).
		:param output_dir:
		"""
		self.image = image
		self.wfsizes = wfsizes
		self.ID = ID
		self.X = X
		self.Y = Y
		self.star_ID = star_ID
		self.star_x = star_x
		self.star_y = star_y
		self.output_dir = output_dir
		self.psfmap = os.path.join(self.output_dir, os.path.basename(self.image)+'.psf.map')

	def mkpsfcat(self, output_dir=None, fcut=20, plot=False, quiet=False, clean=True):
		"""
		Find the stars in a fits image by comparing size of objects to flux. Stars are located in the plane where the objects are very bright but also very compact.

		:param fcut: the flux cut to be applied. It represents the number of brightest objects to be discarded (which are due to very saturated pixels and cannot model the psf properly). Default is 20

		USAGE:
		from astro import PSF_pipeline
		p = PSF_pipeline.psf_pipeline('/home/georgiou/Documents/PAUS/test_image/red_paucam.6988.0957.0149.FT_NB535_NB605.2345042.std.03.fits',1,1,1,1)
		p.mkpsfcat()
		"""
		import heapq
		if output_dir==None: output_dir = self.output_dir
		if not os.path.exists(output_dir): os.makedirs(output_dir)

		sexconfig = os.path.join(os.environ['GAaP'], 'sexconfig')
		for suffix in ['.conv', '.nnw', '.param', '.sex']: shutil.copy(os.path.join(sexconfig, 'default'+suffix), './')
		tmpcat = os.path.basename(self.image) +'_tmp.cat'
		bashCommand = "sex %s -CATALOG_NAME %s" % (self.image, tmpcat)
		if quiet:
			bashCommand += ' -VERBOSE_TYPE QUIET'
			process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=open(os.devnull, 'w'))
		else: process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE, stdin=subprocess.PIPE)
		stdout, stderr = process.communicate()
		print stdout

		x, y, f, ferr, r, a, b, theta, num, flag = np.loadtxt(tmpcat, unpack=True)

		f_at_cut = heapq.nlargest(fcut, f)[-1]

		npcat = np.array([
			(x[i], y[i], f[i], ferr[i], r[i], a[i], b[i], theta[i], num[i], flag[i])
			for i in range(len(x))
		])

		mask1 = (f<f_at_cut) & (f>f_at_cut/30.) & (r>0.5)
		median1 = np.median(npcat[mask1][:,4])

		mask2 =  (f<f_at_cut) & (f>f_at_cut/100.) & (r>0.5 ) & (r<1.5*median1)
		median2 = np.median(npcat[mask2][:,4])

		mask3 = (r>median2/2.) & (r<median2+1.) & (f>0.)
		fmax = np.max(npcat[mask3][:,2])
		mask4 = mask3 & (f>fmax/100.)
		starcat = np.array(sorted(npcat[mask4], key=lambda x: x[2], reverse=True))
		starcat = starcat[int(0.1*len(starcat))+1:len(starcat)]

		self.star_ID = np.array(starcat[:,-2], dtype=np.int)
		self.star_x = starcat[:,0]
		self.star_y = starcat[:,1]

		np.savetxt(os.path.join(output_dir, os.path.basename(self.image)+'.cat'), starcat, fmt='%s %s %s %s %s %s %s %s %i %i')
		print '\nSources: %s \t Stars: %s\n' % (len(npcat), len(starcat))

		if clean:
			for filename in ['default.conv', 'default.nnw', 'default.param', 'default.sex', tmpcat]: os.remove(filename)

		flux_sources, flux_err_sources, radius_sources = f, ferr, r
		flux_stars, flux_err_stars, radius_stars = starcat[:,2], starcat[:,3], starcat[:,4]
		rmin, rmax, fmin, fmax = np.min(radius_stars), np.max(radius_stars), np.min(flux_stars), np.max(flux_stars)

		if plot:
			fig, a = plt.subplots()
			a.errorbar(radius_sources, flux_sources, yerr=flux_err_sources, fmt='bo')
			a.errorbar(radius_stars, flux_stars, yerr=flux_err_stars, fmt='rx')
			a.legend(['Sources','Stars'])
			a.set_xlim([rmin-rmin/10,rmax+rmax/10])
			a.set_ylim([0,fmax+fmax/5])
			a.set_xlabel('FLUX_RADIUS')
			a.set_ylabel('FLUX_AUTO')
			plt.show()
			fig.savefig('size_plot.pdf')

	def gaussianize_psf(self, Hermite_order=10, porder=4, output_dir=None):
		"""
		Construct shapelets-based PSF model from an input star catalogue.

		:param Hermite_order: Maximum order of Hermite polynomials to be used for the shapelets-based PSF model. Default = 10.
		:param porder: Maximum order of polynomial to be used for PSF interpolation in galaxies position. Default = 2, which seems to work fine with PAUS images. Potentially, higher (and more accurate) order can be achieved in fields with more stars.
		:param output_dir: Output directory.
		"""
		if output_dir==None: output_dir=self.output_dir
		if not os.path.exists(output_dir): os.makedirs(output_dir)
		code = os.path.join(os.environ["GAaP"], 'shapelets', 'kk')
		bigim = os.path.join(code, 'bigim')

		stars_file = open('stars_file.cat', 'w')
		for i in range(len(self.star_ID)):
			stars_file.write('%s %s %s %s %s %s %s %s %i %i \n' %(self.star_x[i], self.star_y[i], 1, 0.001, 1, 1, 1, 0, self.star_ID[i], 0))
		stars_file.close()
		psf = 'stars_file.cat'

		print "===================="
		print "PSF-GAUSSIANIZING IMAGE %s" %self.image

		if not os.path.isfile(self.image):
			print '==== INPUT IMAGE %s NOT FOUND ====' %self.image
			exit()
		if not os.path.isfile(psf):
			print '==== PSF CATAOLGUE %s NOT FOUND ====' %psf
			exit()

		psfsh = os.path.join(output_dir, os.path.basename(self.image)+'.psf.sh')
		psfmap = os.path.join(output_dir, os.path.basename(self.image)+'.psf.map')

		orders = open('orders.par', 'w')
		orders.write('%i\n%i' %(Hermite_order, porder))
		orders.close()

		if (not os.path.isfile(psfmap)) or psfmap=="":
			print "BUILDING SHAPELET-BASED PSF MAP"
			if os.path.isfile('inimage.fits'): os.remove('inimage.fits')
			os.symlink(self.image, 'inimage.fits')
			os.system('%s/psfcat2sherr_2 < %s > %s' %(bigim, psf, psfsh))
			os.system('%s/fitpsfmap2 < %s > %s' %(bigim, psfsh, psfmap))
			# The line above spawns the following: bgnoise.dat (bgr and rms), fitpsfmap.ps (?), $im.psf.map, noise-est.ps (noise histogram), psfstarpos.ps (star positions), psfstars.ps (fits to stars), psfsticks.ps (psf whisker plot), starsused.txt (list of stars actually used after QC).
			shutil.move('psfsticks.ps', os.path.join(output_dir, '%s.dirtypsfsticks.ps'%os.path.basename(self.image)))
			shutil.move('psfstars.ps', os.path.join(output_dir, '%s.dirtypsfstars.ps'%os.path.basename(self.image)))
			os.system('%s/showpsfmap < %s' %(code, psfmap))
			shutil.move('psfmap.ps', os.path.join(output_dir, '%s.dirtypsfmap.ps'%os.path.basename(self.image)))
		else: print "==== USING EXISTING PSF MAP %s" %psfmap

		for filename in ['stars_file.cat', 'orders.par', 'bgnoise.dat', 'fitpsfmap.ps', 'noise-est.ps', 'starsused.txt', 'fitkermap.ps', 'gpsfsig.dat', 'kerpos.ps', 'kersticks.ps', 'inimage.fits']:
			#shutil.move(filename, os.path.join(output_dir, filename))
			try: os.remove(filename)
			except OSError: pass
		try: shutil.move('psfstarpos.ps', os.path.join(output_dir, 'psfstarpos.ps'))
		except IOError: pass

	def read_psfmap(self):
		"""
		This reads the .psf.map file and stores the various parameters found in there
		"""
		# Reading only the first two lines here
		f = open(self.psfmap)
		lines = f.readlines()
		self.N, self.beta = map(float, lines[0].split())
		self.klmax, self.xmax, self.ymax = map(float, lines[1].split())
		f.close()
		self.N = int(self.N)
		self.klmax = int(self.klmax)
		# Image centre
		self.xi0 = self.xmax/2.
		self.eta0 = self.ymax/2.

		# Read the rest of the file: the coefficients sab;kl
		self.sabkl = np.loadtxt(self.psfmap, skiprows=3)

	def create_klab_lists(self, Nmax):
		"""
		The following two lists contain the numbering of k and l for the sum k,l from 0,0 up to k+l<=Nmax
		This will result in giving the following numbers: k,l=[00,10,01,20,11,02,30,21,12,03,...]
		"""
		list1 = []
		for k in range(Nmax+1):
			output = k
			while(output>=0):
				list1.append(output)
				output -= 1
		list1=np.array(list1)

		list2 = []
		for l in range(Nmax+1):
			output = 0
			while(output<=l):
				list2.append(output)
				output += 1
		list2=np.array(list2)
		return list1, list2

	def sab(self, xi, eta):
		"""
		Normalised coefficients for the PSF calculation
		(xi, eta)= position on the image. (xi0, eta0)= image's centre. sabkl are the coefficients of the psf map
		klist is the list of k values (so [0,1,0,2...]) and llist the l values (so [0,0,1,0,...])
		"""
		klist, llist = self.create_klab_lists(self.klmax)
		return np.sum(self.sabkl*((2.*xi/(2.*self.xi0)-1.)**klist)*((2.*eta/(2.*self.eta0)-1.)**llist), axis=1)

	def I(self, s):
		"""
		Returns the value of the integral I(i,a)=\int_{-inf}^{inf} e^(-c u^2) H_a(u) u^i du with c=(beta^2+s^2)/s^2 for i=0,1 or 2. Calculated with mathematica 10
		"""
		Int = np.zeros((3,11))

		Int[0,0] = np.sqrt(2.*pi)/np.sqrt(1.+self.beta**2./s**2.)
		Int[0,2] = 2.*(s**2.-self.beta**2.)*np.sqrt(2.*pi)/(np.sqrt(1.+self.beta**2./s**2.)*(s**2.+self.beta**2.))
		Int[0,4] = 12.*(s**2.-self.beta**2.)**2.*np.sqrt(2.*pi)/(np.sqrt(1.+self.beta**2./s**2.)*(s**2.+self.beta**2.)**2.)
		Int[0,6] = 120.*(s**2.-self.beta**2.)**3.*np.sqrt(2.*pi)/(np.sqrt(1.+self.beta**2./s**2.)*(s**2.+self.beta**2.)**3.)
		Int[0,8] = 1680.*(s**2.-self.beta**2.)**4.*np.sqrt(2.*pi)/(np.sqrt(1.+self.beta**2./s**2.)*(s**2.+self.beta**2.)**4.)
		Int[0,10] = 30240.*(s**2.-self.beta**2.)**5.*np.sqrt(2.*pi)/(np.sqrt(1.+self.beta**2./s**2.)*(s**2.+self.beta**2.)**5.)

		Int[1,1] = 2.*np.sqrt(2.*pi)/(1.+self.beta**2./s**2.)**(3./2.)
		Int[1,3] = 12.*np.sqrt(2.*pi)*s**2.*(s**2.-self.beta**2.)/(np.sqrt(1.+self.beta**2./s**2.)*(self.beta**2.+s**2.)**2.)
		Int[1,5] = 120.*np.sqrt(2.*pi)*(s**3.-self.beta**2.*s)**2./(np.sqrt(1.+self.beta**2./s**2.)*(self.beta**2.+s**2.)**3.)
		Int[1,7] = 1680.*np.sqrt(2.*pi)*np.sqrt(1.+self.beta**2./s**2.)*s**4.*(s**2.-self.beta**2.)**3./(self.beta**2.+s**2.)**5.
		Int[1,9] = 30240.*np.sqrt(2.*pi)*np.sqrt(1.+self.beta**2./s**2.)*(s**3.-self.beta**2.*s)**4./(self.beta**2.+s**2.)**6.

		Int[2,0] = np.sqrt(2.*pi)/(1.+self.beta**2./s**2.)**(3./2.)
		Int[2,2] = 2.*np.sqrt(2.*pi)*np.sqrt(1.+self.beta**2./s**2.)*s**4.*(5.*s**2.-self.beta**2.)/(s**2.+self.beta**2.)**3.
		Int[2,4] = 12.*np.sqrt(2.*pi)*np.sqrt(1.+self.beta**2./s**2.)*s**4.*(9.*s**4.-10.*self.beta**2.*s**2.+self.beta**4.)/(s**2.+self.beta**2.)**4.
		Int[2,6] = -120.*np.sqrt(2.*pi)*(self.beta**2.-13.*s**2.)*(self.beta**2.*s-s**3.)**2./(np.sqrt(1.+self.beta**2./s**2.)*(s**2.+self.beta**2.)**4.)
		Int[2,8] = 1680.*np.sqrt(2.*pi)*np.sqrt(1.+self.beta**2./s**2.)*s**4.*(self.beta**2.-17.*s**2.)*(self.beta**2.-s**2.)**3./(s**2.+self.beta**2.)**6.
		Int[2,10] = 30240.*np.sqrt(2.*pi)*np.sqrt(1.+self.beta**2./s**2.)*(21.*s**2.-self.beta**2.)*(s**3.-self.beta**2.*s)**4./(s**2.+self.beta**2.)**7.

		return Int

	def psfmoments(self, i, j, xi, eta, sab, s):
		"""
		Analytical calculation of the PSF moments Q_ij=\int(W(x,y)PSF(x,y) x^i y^j dxdy) at points xi, eta for a gaussian weighting function with size parameter s.
		sab is the coefficients after the summation of the psf map and N is the max orders used for summing the Hermite polynomials.
		see Kuijken+ 2015 Appendix (arXiv:1507.00738)
		"""
		# FIXME: Could this be faster somehow?
		alist, blist = self.create_klab_lists(self.N)
		Sab = np.zeros(len(alist))
		for k in range(len(alist)):
			Sab[k] = self.beta**(i+j+1)*self.I(s)[i,alist[k]]*self.I(s)[j,blist[k]]/np.sqrt(2**(alist[k]+blist[k])*pi*math.factorial(alist[k])*math.factorial(blist[k]))

		#    Sab[k] = beta**(i+j+1)*I(s,beta)[i,alist[k]]*I(s,beta)[j,blist[k]]/np.sqrt(2**(alist[k]+blist[k])*pi*math.factorial(alist[k])*math.factorial(blist[k]))
		return np.sum(sab*Sab)

	def PSF_moments(self, quiet=False, output=None, use_wcs=False):
		"""
		Routine to calculate the PSF moments up to 2nd order
		"""
		if output is None: output=os.path.join(self.output_dir, os.path.basename(self.image)+'_psf_moments.cat')
		self.read_psfmap()

		if not os.path.exists(os.path.dirname(output)): os.makedirs(os.path.dirname(output))
		output_file = open(output, 'w')

		try: len(self.X)
		except:
			self.X = [self.X]
			self.Y = [self.Y]
			self.ID = [self.ID]

		if use_wcs:
		# Write RA,DEC instead of X,Y
			from astropy.wcs import WCS
			w = WCS(self.image)
			ra, dec = w.all_pix2world(self.X, self.Y, 1)

		# start cycle through galaxies
		for obj in range(len(self.X)):
			if not quiet: stdout.write('\r%r out of %r' % (obj+1, len(self.X))); stdout.flush()

			if use_wcs: output_file.write("%s %s %s " % (self.ID[obj], ra[obj], dec[obj]))
			else: output_file.write("%s %s %s " % (self.ID[obj], self.X[obj], self.Y[obj]))
			sab_psf = self.sab(self.X[obj], self.Y[obj])
			for size in self.wfsizes[obj]:
				Q00 = self.psfmoments(0, 0, self.X[obj], self.Y[obj], sab_psf, size)
				Q10 = self.psfmoments(1, 0, self.X[obj], self.Y[obj], sab_psf, size)
				Q01 = self.psfmoments(0, 1, self.X[obj], self.Y[obj], sab_psf, size)
				Q20 = self.psfmoments(2, 0, self.X[obj], self.Y[obj], sab_psf, size)
				Q11 = self.psfmoments(1, 1, self.X[obj], self.Y[obj], sab_psf, size)
				Q02 = self.psfmoments(0, 2, self.X[obj], self.Y[obj], sab_psf, size)
				output_file.write("%s %s %s %s %s %s %s " % (Q00, Q10, Q01, Q20, Q11, Q02, size))
			output_file.write("\n")
		if not quiet: stdout.write("\n")
		output_file.close()

	def Iunw(self):
		'''
		Returns the value of the integral I(i,a)=\int_{-inf}^{inf} e^(-u^2/2) H_a(u) u^i du with for i=0,1 or 2. Calculated with mathematica 10
		'''
		Int = np.zeros((3,11))

		Int[0,0] = np.sqrt(2.*pi)
		Int[0,2] = 2. * np.sqrt(2.*pi)
		Int[0,4] = 12. * np.sqrt(2.*pi)
		Int[0,6] = 120. * np.sqrt(2.*pi)
		Int[0,8] = 1680. * np.sqrt(2.*pi)
		Int[0,10] = 30240. * np.sqrt(2.*pi)

		Int[1,1] = 2. * np.sqrt(2.*pi)
		Int[1,3] = 12. * np.sqrt(2.*pi)
		Int[1,5] = 120. * np.sqrt(2.*pi)
		Int[1,7] = 1680. * np.sqrt(2.*pi)
		Int[1,9] = 30240. * np.sqrt(2.*pi)

		Int[2,0] = np.sqrt(2.*pi)
		Int[2,2] = 10. * np.sqrt(2.*pi)
		Int[2,4] = 108. * np.sqrt(2.*pi)
		Int[2,6] = 1560. * np.sqrt(2.*pi)
		Int[2,8] = 28560. * np.sqrt(2.*pi)
		Int[2,10] = 635040. * np.sqrt(2.*pi)

		return Int

	def psfunweightedmoments(self, i, j, xi, eta, sab, s):
		'''
		Analytical calculation of the PSF unweighted moments Q_ij=\int(PSF(x,y) x^i y^j dxdy) at points xi, eta.
		sab is the coefficients after the summation of the psf map and N is the max orders used for summing the Hermite polynomials.
		see Kuijken+ 2015 Appendix (arXiv:1507.00738)
		'''
		# FIXME: Could this be faster somehow?
		alist, blist = self.create_klab_lists(self.N)
		Sab = np.zeros(len(alist))
		for k in range(len(alist)):
			Sab[k] = self.beta**(i+j+1)*self.Iunw()[i,alist[k]]*self.Iunw()[j,blist[k]]/np.sqrt(2**(alist[k]+blist[k])*pi*math.factorial(alist[k])*math.factorial(blist[k]))
		return np.sum(sab*Sab)

	def PSF_unweighted_moments(self, quiet=False, output=None, use_wcs=False):
		"""
		Routine to calculate the PSF unweighted moments up to 2nd order
		"""
		if output is None: output=os.path.join(self.output_dir, os.path.basename(self.image)+'_psf_moments.cat')
		self.read_psfmap()

		if not os.path.exists(os.path.dirname(output)): os.makedirs(os.path.dirname(output))
		output_file = open(output, 'w')

		try: len(self.X)
		except:
			self.X = [self.X]
			self.Y = [self.Y]
			self.ID = [self.ID]

		if use_wcs:
		# Write RA,DEC instead of X,Y
			from astropy.wcs import WCS
			w = WCS(self.image)
			ra, dec = w.all_pix2world(self.X, self.Y, 1)

		# start cycle through galaxies
		for obj in range(len(self.X)):
			if not quiet: stdout.write('\r%r out of %r' % (obj+1, len(self.X))); stdout.flush()

			if use_wcs: output_file.write("%s %s %s " % (self.ID[obj], ra[obj], dec[obj]))
			else: output_file.write("%s %s %s " % (self.ID[obj], self.X[obj], self.Y[obj]))
			sab_psf = self.sab(self.X[obj], self.Y[obj])
			for size in self.wfsizes[obj]:
				Q00 = self.psfunweightedmoments(0, 0, self.X[obj], self.Y[obj], sab_psf, size)
				Q10 = self.psfunweightedmoments(1, 0, self.X[obj], self.Y[obj], sab_psf, size)
				Q01 = self.psfunweightedmoments(0, 1, self.X[obj], self.Y[obj], sab_psf, size)
				Q20 = self.psfunweightedmoments(2, 0, self.X[obj], self.Y[obj], sab_psf, size)
				Q11 = self.psfunweightedmoments(1, 1, self.X[obj], self.Y[obj], sab_psf, size)
				Q02 = self.psfunweightedmoments(0, 2, self.X[obj], self.Y[obj], sab_psf, size)
				output_file.write("%s %s %s %s %s %s %s " % (Q00, Q10, Q01, Q20, Q11, Q02, size))
			output_file.write("\n")
		if not quiet: stdout.write("\n")
		output_file.close()

"""
USAGE

from astro import PSF_pipeline as psf
import numpy as np
ID, x, y = np.loadtxt('../PAUS/GAaP_no_preGPSFcorr.cat', unpack=True, usecols=(9, 0, 1))
ID = np.arange(1,len(x)+1,1)
wfsizes = np.full(fill_value=5., shape=(len(x),1))
p = psf.psf_pipeline('../PAUS/test_image/red_paucam.6988.0957.0149.FT_NB535_NB605.2345042.std.03.fits', wfsizes, ID, x, y, output_dir='temp')
p.mkpsfcat()
p.gaussianize_psf()
p.PSF_moments()

"""