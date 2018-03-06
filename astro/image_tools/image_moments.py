import numpy as np
from astropy.io import fits
import os
import subprocess
import matplotlib.pyplot as plt
import pdb
# import colormaps
pi = np.pi

class image_moments(object):

	def __init__(self, image, x, y, image_size=None, s=None, order=2):
		self.image = image
		self.x = x
		self.y = y
		self.s = s
		self.order = order
		self.image_size = image_size
		self.image_data = fits.getdata(self.image)
		if image_size == None: self.image_size = self.image_data.shape[0]

		xrange = np.arange(self.x-self.image_size/2+0.5, self.x+self.image_size/2+0.5, step=1)
		yrange = np.arange(self.y-self.image_size/2+0.5, self.y+self.image_size/2+0.5, step=1)
		self.X, self.Y = np.meshgrid(xrange, yrange, sparse=True)
		self.Xint = np.array(self.X-0.5, dtype=int)
		self.Yint = np.array(self.Y-0.5, dtype=int)

		if self.s==None: self.weight_function = 1.
		else: self.weight_function = np.exp( -((self.X-self.x)**2.+(self.Y-self.y)**2.)/(2.*self.s**2.) )

		self.Q = []
		self.Qarray = []
		for k in range(order+1):
			for i in range(k+1):
				self.Qarray.append(self.weight_function*self.image_data[self.Yint, self.Xint]*(self.X-self.x)**(k-i)*(self.Y-self.y)**(i))
				self.Q.append( np.sum( self.weight_function*self.image_data[self.Yint, self.Xint]*(self.X-self.x)**(k-i)*(self.Y-self.y)**(i) ) )
				print k-i, i, self.Q[-1]
		self.Q = np.array(self.Q)
		if self.order>=2:
			self.Q00, self.Q10, self.Q01, self.Q20, self.Q11, self.Q02 = self.Q[:6]
			self.e1 = (self.Q20-self.Q02)/(self.Q20+self.Q02+2*np.sqrt(self.Q20*self.Q02-self.Q11**2))
			self.e2 = (2*self.Q11)/(self.Q20+self.Q02+2*np.sqrt(self.Q20*self.Q02-self.Q11**2))
		fig = plt.figure()
		plt.subplot(331)
		plt.imshow(self.image_data, interpolation='nearest', aspect='equal')
		plt.title(r'$Q_{00}$')
		plt.colorbar()
		plt.subplot(332)
		plt.imshow(self.Qarray[1], interpolation='nearest', aspect='equal')
		plt.title(r'$Q_{10}$')
		plt.colorbar()
		plt.subplot(333)
		plt.imshow(self.Qarray[2], interpolation='nearest', aspect='equal')
		plt.title(r'$Q_{01}$')
		plt.colorbar()
		plt.subplot(334)
		plt.imshow(self.Qarray[3], interpolation='nearest', aspect='equal')
		plt.title(r'$Q_{20}$')
		plt.colorbar()
		plt.subplot(335)
		plt.imshow(self.Qarray[4], interpolation='nearest', aspect='equal')
		plt.title(r'$Q_{11}$')
		plt.colorbar()
		plt.subplot(336)
		plt.imshow(self.Qarray[5], interpolation='nearest', aspect='equal')
		plt.title(r'$Q_{02}$')
		plt.colorbar()
#		plt.show()
		fig.tight_layout()
		fig.savefig('moments_psd_%s.pdf' %self.image_size)
		plt.close()

