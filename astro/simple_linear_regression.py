import numpy as np
from astropy.io import fits
import os
import subprocess
import matplotlib.pyplot as plt
import pdb
pi = np.pi

class simple_linear_regression_old(object):
	def __init__(self, X, Y):
		self.x = X
		self.y = Y
		if len(self.x) != len(self.y): exit('Arrays do not have equal lengths (N_x=%s, N_y=%s)' %( len(X), len(Y)))
		self.N = len(self.x)
	
	def get_means(self):
		self.xmean = np.mean(self.x)
		self.ymean = np.mean(self.y)
		self.xxmean = np.mean(self.x ** 2)
		self.xymean = np.mean(self.x * self.y)
		self.yymean = np.mean(self.y ** 2)
	
	def get_params(self):
		self.b = (self.xymean - self.xmean * self.ymean) / (self.xxmean - self.xmean ** 2)
		self.c = self.ymean - self.b * self.xmean
		self.err = self.y - self.b * self.x - self.c
		self.stderr = np.sqrt(np.sum(self.err ** 2) / self.N)
		self.sigma_b = np.sqrt((np.sum(self.err ** 2) / (self.N - 2)) / np.sum((self.x - self.xmean) ** 2))
		self.sigma_c = self.sigma_b * np.sqrt(self.xxmean)
		self.rxy = (self.xymean - self.xmean * self.ymean) / np.sqrt((self.xxmean - self.xmean ** 2) * (self.yymean - self.ymean ** 2))
	
	def fit(self):
		self.get_means()
		self.get_params()
		
	def plot(self):
		plt.plot(self.x, self.y, 'o')
		plt.plot(self.x, self.b* self.x+self.c, '-')
		plt.plot(self.x, self.x, '--')
		plt.show()

class simple_linear_regression(object):
	def __init__(self, X, Y, weight=None, yerr=None):
		self.X = X
		self.Y = Y
		if len(self.X) != len(self.Y): exit('Arrays do not have equal lengths (N_x=%s, N_y=%s)' %( len(X), len(Y)))
		self.N = len(self.X)
		if weight is None:
			if yerr is None: self.weight = 2./np.ones(self.N)
			else: self.weight = 2./yerr**2. * np.ones(self.N)
		else: self.weight = weight * np.ones(self.N)
		
	def get_params(self):
		self.alpha = np.sum(self.weight * self.X**2.)
		self.beta = np.sum(self.weight)
		self.gamma = np.sum(self.weight * self.X)
		self.p = np.sum(self.weight * self.X * self.Y)
		self.q = np.sum(self.weight * self.Y)
		
		self.b = (self.beta*self.p - self.gamma*self.q)/(self.alpha*self.beta-self.gamma**2.)
		self.c = (self.alpha*self.q-self.gamma*self.p)/(self.alpha*self.beta-self.gamma**2.)
		
	def get_errors(self):
		self.sigma_b = np.sqrt(2 * self.beta / (self.alpha * self.beta - self.gamma ** 2.))
		self.sigma_c = np.sqrt(2 * self.alpha / (self.alpha * self.beta - self.gamma ** 2.))
		self.cov_bc = np.sqrt(2 * self.gamma / (self.alpha * self.beta - self.gamma ** 2.))
		if np.array_equal(self.weight, 2./np.ones(self.N)):
			S2 = np.sum( (self.b * self.X + self.c - self.Y)**2.)/(self.N-1)
			self.sigma_b *= np.sqrt(S2/2)
			self.sigma_c *= np.sqrt(S2 / 2)
			self.cov_bc *= np.sqrt(S2 / 2)
			#self.weight = 2 * np.ones(self.N)/S
			#self.get_params()
			#self.sigma_b = np.sqrt(2 * self.beta / (self.alpha * self.beta - self.gamma ** 2.))
			#self.sigma_c = np.sqrt(2 * self.alpha / (self.alpha * self.beta - self.gamma ** 2.))
			#self.cov_bc = np.sqrt(2 * self.gamma / (self.alpha * self.beta - self.gamma ** 2.))
			
	def fit(self):
		self.get_params()
		self.get_errors()
		
	def plot(self):
		plt.plot(self.X, self.Y, 'o')
		plt.plot(self.X, self.b* self.X+self.c, '-')
		plt.plot(self.X, self.X, '--')
		plt.show()
	
def test():
	N = 100
	b = 5.3
	c = 0.5
	X = np.linspace(1, 10, N)
	err = 0.2 * np.random.normal(size=N)
	Y = b*np.linspace(1,10,N)+c+err

	regr = simple_linear_regression(X,Y, yerr=1.0)
	regr.fit()
	print regr.b, regr.sigma_b, regr.c, regr.sigma_c, regr.cov_bc
	regrold = simple_linear_regression_old(X,Y)
	regrold.fit()
	print regrold.b, regrold.sigma_b, regrold.c, regrold.sigma_c

'''
from scipy.optimize import curve_fit

def func(x, b, c):
	return b*x+c

popt, pcov = curve_fit(func, X, Y)
perr = np.sqrt(np.diag(pcov))
'''