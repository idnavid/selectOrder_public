import numpy as np
import scipy.io
from nilearn import image
from nilearn import plotting
from nilearn import masking
from nilearn import signal 
import nibabel as nib
from nipype.interfaces.spm import SliceTiming

def load_data(filename):
	"""
		Load data stored in matrix format as text file. 
	"""
	return np.loadtxt(filename,delimiter=',')

def plot_fmri(img):
	for t in range(img.shape[-1]):
		img_t = image.index_img(img, t)
		plotting.plot_stat_map(img_t)
		plotting.show()

def load_fmri(filename,maskname='',sigma=3):
	"""
		Reads 4D fmri data. 
		Smooths using 3D Gaussian filter
		Applies mask to data: 
			- If mask not provided, calculates binary mask
		returns fmri data matrix (Time x Voxels)
	"""
	img = nib.load(filename)
	print(img.shape)
	rep_time = img.header.get_zooms()[-1]
	img = image.smooth_img(img,sigma)
	if maskname!='':
		img_mask = nib.load(maskname)
	else:
		print('Mask not provided. Calculating mask ...')
		img_mask = masking.compute_background_mask(img)
	img = masking.apply_mask(img,img_mask)
	print('Mask applied!')
	print('Detrending data!')
	img = signal.clean(img,
					   detrend=True,
					   high_pass=0.01,
					   standardize=False,
					   t_r=rep_time)
	return img

def save_fmri(filename,X):
	"""
		Save fmri data in MATLAB format. 
		This is used for debugging purposes. 
		The output is saved in .mat format 
		under the name filename. 
		The matrix X is n_t by n_s dimensional
		n_t: number time samples
		n_s: number of voxels
	"""
	scipy.io.savemat(filename, 
					 mdict={'out': X}, oned_as='row')
	return 


if __name__=='__main__':
	# fmri_file = '../../data/hcp_data/HCP_rsfMRI_small/100307/rfMRI_REST1_LR.nii.gz'
	fmri_file = '../../data/hcp_data/rfMRI_54.nii.gz'
	load_fmri(fmri_file)
	