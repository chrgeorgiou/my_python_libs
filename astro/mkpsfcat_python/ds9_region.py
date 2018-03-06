def ds9_region_file(x, y, a, b=None, theta=None, spatial_units='p', size_units='"', pa_units='d', color='Green', text=None, width=1, filename='./ds9_region_file.reg'):
	'''
	Create a region file that can be loaded to ds9. Regions can be either circles or ellipses.

	x, y: are the spatial coordinates where the region will be created. If b is not given, a coresponds to the radius of the circular region to be drawn. Otherwise, the region will be elliptical and a, b will be the semi-major and semi-minor axes, with theta the position angle. 

	spatial_units: define the units of x and y. They can be 'p' for physical pixels, 'i' for image pixels, 'd' for degrees and 'r' for radians.

	size_units: are the same as before, defining the units of a (and b). They can be 'p', 'i' or 'd' but also " for arcsec and ' for arcmin.

	pa_units: define the units of the position angle theta. The angle is measured from the Y-axis of the image counter-clockwise. 

	color: defines the color of the region (Black, White, Red, Green, Blue, Cyan, Magenta or Yellow).

	text: if defined, the text will appear right on top of the region in the same color as the region.

	width: defines how thick the region boundary is. Can be 1,2,3 or 4.
	'''
	import numpy as np
	# check if input is array. If not, make it a 1-dimension list.
	try: len(x)
	except: x, y, a, b, theta = [x], [y], [a], [b], [theta]
	try: len(color)
	except: color = np.full_like(x, color, dtype=str)

	reg_file = open(filename, 'w')
	for obj in range(len(x)):
		if (b==None) or (b[obj]==None) or (b[0]==None): reg_file.write('circle %s%s %s%s %s%s # color=%s width=%i ' % (x[obj], spatial_units, y[obj], spatial_units, a[obj], size_units, color, width) )
		else: reg_file.write('ellipse %s%s %s%s %s%s %s%s %s%s # color=%s width=%i ' % (x[obj], spatial_units, y[obj], spatial_units, a[obj], size_units, b[obj], size_units, theta[obj], pa_units, color[obj], width) )
		if (text==None) or (text[0]==None): reg_file.write('\n')
		else:
			try: reg_file.write('text={%s}\n' % text[obj] )
			except: reg_file.write('text={%s}\n' % text )


