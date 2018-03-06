import numpy as np
from scipy import spatial

__all__ = ['spatial_matching', 'index_matching']

def spatial_matching(x1, y1, x2, y2, eps ,k=1):
	'''
	:param x1: position 1 of object in catalogue 1 (x or RA)
	:param y1: position 2 of object in catalogue 1 (y or DEC)
	:param x2: position 1 of object in catalogue 2 (x or RA)
	:param y2: position 2 of object in catalogue 2 (y or DEC)
	:param eps: maximum distance within which to look for a match
	:return: an array of dimension 2 x len(x2). The array contains the position of the spatially matched object from catalogue 2 to 1. The first element of the array is the element in catalogue 1 and the second is the element in catalogue 2. If array[i,0]==len(x1) then the object 'i' does not have a match within distance 'eps'.

	USAGE:
	# Match image coordinates to catalogue coordinates with a accepted error of "1.0"
	index = spatial_matching(x_im, y_im, x_cat, y_cat, 1.0)
	# mask containing the matched objects for both coordinates
	matched_objects = index[:,0] != len(x_im) #if the first index is equal to the array's lenght it means the object is not matched
		if len(ra_k[mask]) == 0:
		exit('No matching objects found')

	image_ids = index[mask][:,0] # ids of matched objects in the image
	catalogue_ids = index[mask][:,1] # ids of matched objects in the catalogue

	x_im[image_ids] # x_im for the matched objects
	'''
	c1=np.dstack([x1,y1])[0]
	c2=np.dstack([x2,y2])[0]
	#Creating tree hull and building the tree
	tree=[]
	Tree=spatial.cKDTree(c1,leafsize=128)
	#Querying tree
	for object in c2:
		tree.append(Tree.query(object,k=k,distance_upper_bound=eps))
	#Creating final output array.
	index_list=[]
	#Going through all matches in the tree
	for i in range(0,len(tree)):
		element=[tree[i][1],i]
		index_list.append(element)
	index_list=np.array(index_list)
	return index_list

def index_matching(ID1, ID2):
	'''
	Match the indices of 2 arrays.
	:return: An array of indices. These correspond to the positions of the elements of the smallest array to the largest array.
	
	USAGE
	ID1 = np.array([ 10001.,  10002.,  10003.,  10004.,  10005.,  10006.,  10007., 10008.])
	ID2 = np.array([ 10004.,  10002.,  10001.,  10005.])

	mask = index_matching(ID1,ID2)
	print np.array_equal(ID1[mask], ID2)
	'''
	if len(ID1)<len(ID2): ID1, ID2 = ID2, ID1
	sort_idx = ID1.argsort()
	return sort_idx[np.searchsorted(ID1, ID2, sorter = sort_idx)]

def spatial_matching_improved(x1, y1, x2, y2, eps):
	c1=np.dstack([x1,y1])[0]
	c2=np.dstack([x2,y2])[0]
	Tree=spatial.cKDTree(c1,leafsize=128)
	idx1, idx2, dist = [], [], []
	for i2, object in enumerate(c2):
		d, i1 = Tree.query(object, k=[1], distance_upper_bound=eps)
		if np.isinf(d[0]): continue
		idx1.append(i1[0])
		idx2.append(i2)
		dist.append(d[0])
	idx1 = np.array(idx1, dtype=np.int)
	idx2 = np.array(idx2, dtype=np.int)
	dist = np.array(dist, dtype=np.float)
	return idx1, idx2, dist
