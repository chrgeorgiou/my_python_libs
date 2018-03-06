import argparse
import os
import subprocess
import numpy as np
import heapq
import matplotlib.pyplot as plt

def mkpsfcat(inimage, sex_script, fcut=20):
	'''
	Find the stars in a fits image by comparing size of objects to flux. Stars are located in the plane where the objects are very bright but also very compact. 
	inimage: the location of the input image
	sex_script: the location of the Source Extractor binary file. 
	fcut: the flux cut to be applied. It represents the number of brightest objects to be discarded (which are due to very saturated pixels and cannot model the psf properly). Default is 20
	'''
	print '\nFinding stars in image: %s\n' % inimage

	bashCommand = sex_script+" "+inimage+" -CATALOG_NAME tmp.cat"
	process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE, stdin=subprocess.PIPE)
	stdout, stderr = process.communicate()

	x, y, f, ferr, r, a, b, theta, num, flag = np.loadtxt('tmp.cat', unpack=True)

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

	np.savetxt('./'+os.path.basename(inimage)+'.cat', starcat, fmt='%s %s %s %s %s %s %s %s %i %i')
	print '\nSources: %s \t Stars: %s\n' % (len(npcat), len(starcat))

	flux_sources, flux_err_sources, radius_sources = f, ferr, r
	flux_stars, flux_err_stars, radius_stars = starcat[:,2], starcat[:,3], starcat[:,4]
	rmin, rmax, fmin, fmax = np.min(radius_stars), np.max(radius_stars), np.min(flux_stars), np.max(flux_stars)

	fig, a = plt.subplots()
	a.errorbar(radius_sources, flux_sources, yerr=flux_err_sources, fmt='bo')
	a.errorbar(radius_stars, flux_stars, yerr=flux_err_stars, fmt='rx')
	a.legend(['Sources','Stars'])
	a.set_xlim([rmin-rmin/10,rmax+rmax/10])
	a.set_ylim([0,fmax+fmax/5])
	a.set_xlabel('FLUX_RADIUS')
	a.set_ylabel('FLUX_AUTO')
#	a.set_title(inimage.replace('.bgsub',''))
	#plt.show()
	plt.savefig('size_plot.pdf')

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Identify stars in an image using an iterative median calculation. --requires SExtractor installed and default files present')
	parser.add_argument('-i', help='Insert the location of the input image.', dest='image')
	parser.add_argument('-f', type=int, help='The nth brightest source to cut from. Default = 20', dest='fcut', default=20)
	args = parser.parse_args()

	mkpsfcat(args.image, 'sex', fcut=args.fcut)

	from ds9_region import ds9_region_file
	x, y, a, b, theta = np.loadtxt('./'+os.path.basename(args.image)+'.cat', unpack=True, usecols=(0,1,5,6,7))
	ds9_region_file(x, y, a, b, theta)

