import os
import subprocess
import shutil
import numpy as np
from sys import stdout
import math
pi = np.pi

# FIXME: Make it work with an external star catalogue
class shapelets_psf(object):

	def __init__(self, image, psfmap, psfmoments_file, ID, x, y, wfsize):
		'''
		Class that calculates the psf model of a .fits image based on the shapelets method. For details see Refregier 2001, Kuijken+ 2015 Appendix (arXiv:1507.00738).
		:param image: image file to be processed
		:param psfmap: output psf map file. If None, it is "image.psf.map"
		:param psfmoments_file: output psf moments file. If None, it is "image.psf.mom"
		:param ID: ID number(s) of the points where the psf moments are calculated
		:param x: x-axis position of the above object(s)
		:param y: y-axis position of the above object(s)
		:param wfsize: Array of shape ( len(ID), N ) where N is the number of different weight function sizes to be computed for every object.
		'''
		self.image = image
		self.imagename = os.path.basename(image)
		self.psfmap = psfmap
		self.psfmoments_file = psfmoments_file
		self.ID, self.x, self.y, self.wfsize = ID, x, y, wfsize
	
	def __getattr__(self, item):
		return None

	def call_script(self, script):
		'''
		Call external mkpsfcat.csh OR gpsf.csh script. 'inimage' is the input image and 'script' is the location of the script.
		'''
		bashCommand = script+" "+self.image+" "+self.imagename+".cat"
		FNULL = open(os.devnull, 'w')
		process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)#, stdin=subprocess.PIPE, stderr=FNULL)
		process.wait()
		# stdout, stderr = process.communicate()
		# print stdout
		
	def create_psf_map(self, starcat='create'):
		if not os.path.exists('./deleteme/'+self.imagename+'.dir'): os.makedirs('./deleteme/'+self.imagename+'.dir')
		os.chdir('./deleteme/'+self.imagename+'.dir')
		realpath = os.path.realpath(__file__)
		class_dir = os.path.dirname(realpath)
		mkpsfcat_script = os.path.join( class_dir, 'scripts', 'gaapsw', 'mkpsfcat.csh')
		gpsf_script = os.path.join( class_dir, 'scripts', 'gaapsw', 'gpsf.csh')

		if not os.path.isfile('default.conv'): shutil.copy(os.path.join(class_dir, 'scripts', 'gaapsw', 'sexconfig', 'default.conv'), os.getcwd())
		if not os.path.isfile('default.nnw'): shutil.copy(os.path.join(class_dir, 'scripts', 'gaapsw', 'sexconfig', 'default.nnw'), os.getcwd())
		if not os.path.isfile('default.param'): shutil.copy(os.path.join(class_dir, 'scripts', 'gaapsw', 'sexconfig', 'default.param'), os.getcwd())
		if not os.path.isfile('default.sex'): shutil.copy(os.path.join(class_dir, 'scripts', 'gaapsw', 'sexconfig', 'default.sex'), os.getcwd())

		if starcat == 'create': self.call_script(mkpsfcat_script)
		self.call_script(gpsf_script)
		shutil.move(self.image+'.psf.map', self.psfmap)
		os.remove(self.image+'.dirtypsfmap.ps')
		os.remove(self.image + '.dirtypsfstars.ps')
		os.remove(self.image + '.dirtypsfsticks.ps')
		os.remove(self.image + '.psf.sh')
		os.chdir('../../')
		shutil.rmtree('./deleteme/'+self.imagename+'.dir')

	def read_psfmap(self):
		'''
		This reads the .psf.map file and stores the various parameters found in there
		'''
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
		'''
		The following two lists contain the numbering of k and l for the sum k,l from 0,0 up to k+l<=Nmax
		This will result in giving the following numbers: k,l=[00,10,01,20,11,02,30,21,12,03,...]
		'''
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
		'''
		Normalised coefficients for the PSF calculation
		(xi, eta)= position on the image. (xi0, eta0)= image's centre. sabkl are the coefficients of the psf map
		klist is the list of k values (so [0,1,0,2...]) and llist the l values (so [0,0,1,0,...])
		'''
		klist, llist = self.create_klab_lists(self.klmax)
		return np.sum(self.sabkl*((2.*xi/(2.*self.xi0)-1.)**klist)*((2.*eta/(2.*self.eta0)-1.)**llist), axis=1)

	def I(self, s):
		'''
		Returns the value of the integral I(i,a)=\int_{-inf}^{inf} e^(-c u^2) H_a(u) u^i du with c=(beta^2+s^2)/s^2 for i=0,1 or 2. Calculated with mathematica 10
		'''
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

	def psfmoments(self, i, j, xi, eta, sab, s):
		'''
		Analytical calculation of the PSF moments Q_ij=\int(W(x,y)PSF(x,y) x^i y^j dxdy) at points xi, eta for a gaussian weighting function with size parameter s.
		sab is the coefficients after the summation of the psf map and N is the max orders used for summing the Hermite polynomials.
		see Kuijken+ 2015 Appendix (arXiv:1507.00738)
		'''
		# FIXME: Could this be faster somehow?
		alist, blist = self.create_klab_lists(self.N)
		Sab = np.zeros(len(alist))
		for k in range(len(alist)):
			Sab[k] = self.beta**(i+j+1)*self.I(s)[i,alist[k]]*self.I(s)[j,blist[k]]/np.sqrt(2**(alist[k]+blist[k])*pi*math.factorial(alist[k])*math.factorial(blist[k]))

		#    Sab[k] = beta**(i+j+1)*I(s,beta)[i,alist[k]]*I(s,beta)[j,blist[k]]/np.sqrt(2**(alist[k]+blist[k])*pi*math.factorial(alist[k])*math.factorial(blist[k]))
		return np.sum(sab*Sab)

	def PSF_moments(self, quiet=False):
		'''
		Routine to calculate the PSF moments up to 2nd order
		'''
		self.read_psfmap()
		
		output_directory = os.path.dirname(self.psfmoments_file)
		if not os.path.exists(output_directory): os.makedirs(output_directory)
		output = open(self.psfmoments_file, 'w')

		try: len(self.x)
		except:
			self.x = [self.x]
			self.y = [self.y]
			self.ID = [self.ID]

		# Hard coded write RA,DEC instead of X,Y
		# from astropy.wcs import WCS
		# w = WCS(self.image)
		# ra, dec = w.all_pix2world(self.x, self.y, 1)

		# start cycle through galaxies
		for obj in range(len(self.x)):
			if not quiet: stdout.write('\r%r out of %r' % (obj+1, len(self.x))); stdout.flush()

			#output.write("%s %s %s " % (self.ID[obj], ra[obj], dec[obj]))
			output.write("%s %s %s " % (self.ID[obj], self.x[obj], self.y[obj]))
			sab_psf = self.sab(self.x[obj], self.y[obj])
			for size in self.wfsize[obj]:
				Q00 = self.psfmoments(0, 0, self.x[obj], self.y[obj], sab_psf, size)
				Q10 = self.psfmoments(1, 0, self.x[obj], self.y[obj], sab_psf, size)
				Q01 = self.psfmoments(0, 1, self.x[obj], self.y[obj], sab_psf, size)
				Q20 = self.psfmoments(2, 0, self.x[obj], self.y[obj], sab_psf, size)
				Q11 = self.psfmoments(1, 1, self.x[obj], self.y[obj], sab_psf, size)
				Q02 = self.psfmoments(0, 2, self.x[obj], self.y[obj], sab_psf, size)
				output.write("%s %s %s %s %s %s %s " % (Q00, Q10, Q01, Q20, Q11, Q02, size))
			output.write("\n")
		if not quiet: stdout.write("\n")
		output.close()

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

	def PSF_unweighted_moments(self, quiet=False):
		'''
		Routine to calculate the PSF unweighted moments up to 2nd order
		'''
		self.read_psfmap()
		
		output_directory = os.path.dirname(self.psfmoments_file)
		if not os.path.exists(output_directory): os.makedirs(output_directory)
		output = open(self.psfmoments_file, 'w')
		
		try:
			len(self.x)
		except:
			self.x = [self.x]
			self.y = [self.y]
			self.ID = [self.ID]
		
		# Hard coded write RA,DEC instead of X,Y
		# from astropy.wcs import WCS
		# w = WCS(self.image)
		# ra, dec = w.all_pix2world(self.x, self.y, 1)
		
		# start cycle through galaxies
		for obj in range(len(self.x)):
			if not quiet: stdout.write('\r%r out of %r' % (obj + 1, len(self.x))); stdout.flush()
			
			# output.write("%s %s %s " % (self.ID[obj], ra[obj], dec[obj]))
			output.write("%s %s %s " % (self.ID[obj], self.x[obj], self.y[obj]))
			sab_psf = self.sab(self.x[obj], self.y[obj])
			for size in self.wfsize[obj]:
				Q00 = self.psfunweightedmoments(0, 0, self.x[obj], self.y[obj], sab_psf, size)
				Q10 = self.psfunweightedmoments(1, 0, self.x[obj], self.y[obj], sab_psf, size)
				Q01 = self.psfunweightedmoments(0, 1, self.x[obj], self.y[obj], sab_psf, size)
				Q20 = self.psfunweightedmoments(2, 0, self.x[obj], self.y[obj], sab_psf, size)
				Q11 = self.psfunweightedmoments(1, 1, self.x[obj], self.y[obj], sab_psf, size)
				Q02 = self.psfunweightedmoments(0, 2, self.x[obj], self.y[obj], sab_psf, size)
				output.write("%s %s %s %s %s %s %s " % (Q00, Q10, Q01, Q20, Q11, Q02, size))
			output.write("\n")
		if not quiet: stdout.write("\n")
		output.close()


# FIXME: This needs to be somewhere else
'''
			if self.sizes[-1]=='kpc':# FIXME: This does not seem to work properly
				#Turn physical size values to pixel size values
				from astropy.cosmology import FlatLambdaCDM
				from astropy import units as u
				self.H0 = 70
				self.Om0 = 0.315
				self.pixel_scale = 0.214 #arcsec/pixel
				cosmo = FlatLambdaCDM(H0=self.H0, Om0=self.Om0)
				D_A = cosmo.angular_diameter_distance(self.z[obj]).to(u.kpc) #this is in kpc
				self.sizerange = 206265*self.sizerange/(2*D_A.value*self.pixel_scale) # 1 rad = 206265 arcsec
			elif self.sizes[-1]=='isophotes':
				self.sizerange = [self.sizes[k]*self.segrad[obj] for k in range(len(self.sizes)-1)]


'''
