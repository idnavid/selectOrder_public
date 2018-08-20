# Create Simulated FMRI-like data as a mixture of spatial and temporal components plus noise. 
# Navid Shokouhi, 
# November 2017

import numpy as np 
import scipy.ndimage 
import pylab 
import time

from tools.estimate_dim import *
from tools.fmri_tools import *

def normalize(X):
	"""
		This normalization hurts performance. 
	"""
	X = X - np.mean(X,axis=1)[:,np.newaxis]
	X = X/np.std(X,axis=0)[np.newaxis,:];
	X = X - np.mean(X,axis=0)[np.newaxis,:]
	return X

def normalize_pca(X):
	n_t,n_s = X.shape
	print(X.shape)
	X = X - np.mean(X,axis=1)[:,np.newaxis]
	V, U = np.linalg.eigh(np.dot(X,X.T)/n_s)
	V = np.abs(V)
	print(V.shape,U.shape)
	X1 = np.dot(U,X)
	X2 = np.dot(np.diag(V**(-0.5)),X1)
	X3 = np.dot(U.T,X2)
	return np.real(X3)

def spatial_comp(n,m,positions = 'uniform'):
	"""
		Creates m matrices of m components, 
		each in an nxn matrix, which is vectorized 
		in the output. 
		comp stands for component. 
		n: number of spatial samples along each axis (x,y) 
		m: number of components 

		NOTE: I'll start out by creating blocks of ones
			  in an all-zero matrix. Later on I might 
			  add blobs insteads of blocks.
	"""
	X = np.zeros((n,n))
	Xvec = np.zeros((m,n*n))+1
	if positions=='random':
		comp_positions = np.random.randint(n, size=(2,m))
	elif positions=='uniform':
		comp_distance = int(n/np.sqrt(m))
		comp_positions = np.array([[int((i)*(comp_distance) + comp_distance/2) ,
						 			int((j)*(comp_distance) + comp_distance/2)]
								   for i in range(int(np.sqrt(m))) 
								   for j in range(int(np.sqrt(m)))]).T
	comp_width = int(n/(2*np.sqrt(m)))
	for i in range(m):
		xloc_i = comp_positions[:,i][0]
		yloc_i = comp_positions[:,i][1]
		X = X*0.
		X[max(0,xloc_i-int(comp_width/2)):min(n,xloc_i+int(comp_width/2)),
		  max(0,yloc_i-int(comp_width/2)):min(n,yloc_i+int(comp_width/2))] = 1
		Xvec[i,:] = X.reshape((1,n*n))
	X = np.sum(Xvec,axis=0)
	return Xvec


def temporal_comp(t,m):
	"""
		Create matrix of temporal components. 
		each column is a function describing
		the temporal behavior of the spatial
		components from spatial_comp. 
		t: number of time samples 
		m: number of components 
	"""
	Y = np.zeros((t,m))
	Amp = 1
	for i in range(m):
		Y[:,i] = np.random.randn(t)
		# Y[:,i] = Amp*np.sin(2*np.pi*np.arange(t)*((i+.1)/30.)) + 0.1*np.random.randn(t)
	return Y

def mix_comp(X_tc,X_sc):
	Y = np.dot(X_tc,X_sc)
	n_t = X_tc.shape[0]
	n_s = int(np.sqrt(X_sc.shape[1]))
	for i in range(n_t):
		Y_i = Y[i,:].reshape((n_s,n_s))
		Y[i,:] = Y_i.reshape(Y[i,:].shape)
	Y_bias = np.mean(Y,axis=0)
	Y = Y + Y_bias[np.newaxis,:]
	return Y

		
	
