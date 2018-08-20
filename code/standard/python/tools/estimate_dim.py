# -*- coding: utf-8 -*-

# Eigenvalue adjustment
import numpy as np 
import pylab 
import scipy.special

def commutation(m,n):
	"""
		Implement commutation matrix, K. 
		By definition, a commutation matrix is 
		a form of permutation matrix that has 
		the following property:
				K*vec(A) = vec(A.T)
	"""
	K = np.zeros((m*n,m*n))
	I = np.eye(m*n)
	k = 0;
	for i in range(n):
		for j in range(i,m*n,n):
			K[k,:] = I[j,:]
			k += 1
	# Checking implementation:
	# A = np.random.rand(m,n)
	# A_vec = np.reshape(A,(m*n,1))
	# At_vec = np.reshape(A.T,(m*n,1))
	# print(np.linalg.norm(np.dot(K,A_vec) - At_vec))
	return K

def adjust_eigenspect(V,n,N):
	"""
		Adjust eigenspectrum. Adopted from the implementation 
		in ICA-TB by Li et al., 2007
		
		[1] Li, Yi‐Ou, Tülay Adalı, and Vince D. Calhoun. 
		 Estimating the number of independent components for functional magnetic resonance 
                 imaging data. Human brain mapping 28.11 (2007): 1251-1266.

	"""
	r = float(n)/N
	pi = np.pi 
	bp = (1+np.sqrt(r))**2
	bm = (1-np.sqrt(r))**2
	
	# QUESTION: where does 5n-1 come from?
	vv = np.arange(bm,bp,(bp-bm)/(5.*n - 1.))
	gv = np.divide(1.,(2*pi*r*vv))*np.sqrt(np.multiply(vv-bm,bp-vv))
	
	gv_cum = np.zeros(gv.shape)
	for i in range(len(gv)):
	    gv_cum[i] = np.sum(gv[:i])
	
	gv_cum = gv_cum/np.max(gv_cum)
	
	lam_emp = np.zeros(np.shape(V))
	for i in range(n):
		# QUESTION: Why do we use i_norm to find the solut
	    i_norm = float(i)/n
	    min_index = np.argmin(np.abs(i_norm-gv_cum))
	    lam_emp[i] = vv[min_index]
	lam_emp = lam_emp[::-1]
	Vadj = np.divide(V,lam_emp)
	return Vadj


def loglikelihood(V,k,n):
	"""
		log-likelihood function, as presented in 
		Wax and Kailath's paper: 
		Wax, M., and Kailath, T. (1985) "Detection of 
		signals by information theoretic criteria." IEEE TASSP. 
	"""
	sigma = np.mean(V[k:]**2)
	return np.log(np.prod(np.power(V[k:],1./(n-k)))/sigma)

def loglikelihood_wk(V,k,n): 
	V = np.abs(V)
	return (n-k)*(np.mean(np.log(V[k:])) - np.log(np.mean(V[k:])))


def loglikelihood_alt(V,k,n):
	V = np.abs(V)
	return n*(np.log(2*np.pi)+1) + (n-k)*np.log(np.mean(V[k:])) + np.sum(np.log(V[:k]))

def loglikelihood_chi2(S,U,K,k,n):
	"""
		S: sample covariance matrix, pxp
		U: eigenvectors of S
		K: pxp commutation matrix
		k: assumed number of components
		n: number of samples
	"""
	p = U.shape[0]
	Uk = U[:,k:]
	Q = np.dot(Uk,Uk.T)
	S_proj = np.matrix.flatten(np.dot(np.dot(Q,S),Q))
	Q_kron = np.kron(Q,Q)
	kappa = np.eye(p**2)+K
	V = np.dot(np.eye(p**2)+K,np.kron(S,S))
	# Omega = np.dot(Q_kron,np.dot(V,Q_kron))
	# Omega_inv = np.linalg.pinv(Omega)
	V_inv = np.linalg.pinv(V)
	Omega_inv = np.dot(Q_kron,np.dot(V_inv,Q_kron))
	a1 = np.dot(S_proj.T,np.dot(Omega_inv,S_proj))
	S_proj = np.dot(np.dot(Q,S),Q)
	tri_idx = np.triu_indices(p)
	X = np.zeros((p,p))
	X[tri_idx] = 1
	S_proj = X*S_proj
	S_proj = np.matrix.flatten(S_proj)
	a2 = np.dot(S_proj.T,np.dot(Omega_inv,S_proj))	
	return a2


def loglikelihood_bayes(V,k,n):
	"""
		log-likelihood term, estimated from the Bayesian 
		approach for PPCA, presented in:
		Seghouane, A.K, Cichocki, A., 2006, "Bayesian estimation of the 
		number of principal components", Elsevier Signal Processing. 
	"""
	sigma = np.mean((V[k:]))
	arithmetic_term = (sigma)**(n-k)
	geometric_term = np.prod(V[:k])
	return np.log(arithmetic_term*geometric_term)



def laplace_approximation(V,n,N):
	"""
	This implementation was derived from 
	https://www.fmrib.ox.ac.uk/datasets/techrep/tr02cb1/tr02cb1/node5.html
	which uses the paper:
		Minka, T. (2000). Automatic choice of dimensionality for PCA.
		Technical Report 514, MIT.

	"""
	V = (V+1e-7)*(N/(N-1))


	V_residuals = np.zeros(V.shape)
	for i in range(n):
		V_residuals[i] = np.sum(np.log(np.abs(V[i] - V[i+1:])))
	V_residuals = np.cumsum(V_residuals)

	V_inv = np.power(V,-1.)
	V_inv = V_inv[:,np.newaxis]
	V_inv_diff = V_inv - V_inv.T 
	V_inv_diff[np.where(V_inv_diff<=0)] = 1
	V_inv_diff = np.log(V_inv_diff)
	V_inv_diff = np.cumsum(V_inv_diff,axis=0)
	V_inv_diff = np.sum(V_inv_diff,axis=1)

	logV = np.log(V)
	cumsum_logV = np.cumsum(logV)

	cumsum_V = np.cumsum(V)

	ks = np.arange(1,n)
	z = np.log(2) + (n-ks+1)/2*np.log(np.pi) - scipy.special.gammaln((n-ks+1)/2);
	cumsum_z = np.cumsum(z)
	p = np.zeros((1,n-1))
	v_prev = 0.
	for k in ks:
		i = k-1
		v = (cumsum_V[-1] - cumsum_V[k])/(n-k) 
		if v<=0:
			v = v_prev
		# v[np.where(v<=0)] = np.min(v[np.where(v>0)])
		p[0,i] = -N/(2.*cumsum_logV[k]) + (-N*(n-k)/2.)*np.log(v)
		p[0,i] = p[0,i] - cumsum_z[i] - (k/2.)*np.log(N)
		h = V_inv_diff[i] + V_residuals[i]
		h = h + (n-k)*np.sum(np.log(np.abs(1./v - V_inv[:k])))
		m = n*k-k*(k+1.)/2
		h = h + m*np.log(N);
		p[0,i] = p[0,i] + (m+k)/2.*np.log(2*np.pi) - h/2.
		p[0,i] = p[0,i] + 1.5*k*np.log(2.)
		p[0,i] = p[0,i] - 0.5*np.log(n-k)
		v_prev = v
	return p

def order_objective_function(V,n,N,penalty=1):
	obj = {'aic':np.zeros([1, n-1]),
		   'kic':np.zeros([1, n-1]),
		   'mdl':np.zeros([1, n-1]),
		   'lap':np.zeros([1,n-1]),
		   'alt-1':np.zeros([1, n-1]),
		   'alt-2':np.zeros([1, n-1]),
		   'alt-3':np.zeros([1, n-1])
		   }
	obj['lap'] = laplace_approximation(V,n,N)
	i = 0
	for k in range(1,n):
		mloglikelihood = 0.5*N*(n-k)*loglikelihood(V,k,n) 
		numfreeparams = 1 + 0.5*k*(2*n-k+1)
		obj['aic'][0,i] =  -2*mloglikelihood + 2*numfreeparams
		obj['kic'][0,i] =  -2*mloglikelihood + 3*numfreeparams
		obj['mdl'][0,i] =  -1*mloglikelihood + .5*numfreeparams
		# NOTE: (N-1.)*np.sum(np.log(V[k:])) also works for the alternative 
		# 		criteria. The idea for using log() came from Lawley's 1956 biometrika paper. 
		obj['alt-1'][0,i] =  (N-1.)*np.sum(V[k:]) - 2*(n-k)*(n-k+1.)/2.		
		obj['alt-2'][0,i] =  (N-1.)*np.sum(V[k:]) - np.log(N)*(n-k)*(n-k+1.)/2.
		obj['alt-3'][0,i] =  (N-1.)*np.sum(V[k:]) - np.log(np.log(N))*(n-k)*(n-k+1.)/2.
		i+=1
	return obj

def order_objective_function_chi2(S,n,N):
	obj = np.zeros([1,n-1])
	K = commutation(n,n)
	V,U = np.linalg.eig(S)
	idx = np.argsort(V)[::-1]
	V = V[idx]
	U = U[:,idx]
	i = 0
	for k in range(1,n):
		sigma2 = np.mean(V[k:])
		S_norm = S - sigma2*np.eye(n)
		mloglikelihood = loglikelihood_chi2(S_norm,U,K,k,N) 
		obj[0,i] =  (N-1.)*mloglikelihood - np.log(N)*(n-k)*(n-k+1.)/2.
		i+=1
	return obj

def find_minimum(obj):
	"""
		Often, it is difficult to find the local minimum of 
		the criterion (objective) function that corresponds
		to the estimated dimension. This algorithm is designed
		to find that peak.  

		The approach here is to first locate all the local 
		minimum using the first difference. and then finding 
		abrupt changes in the objective slope using log. 
	"""
	# pylab.plot(np.diff(np.log(np.abs(np.diff(obj)))))
	# pylab.show()

	delta_obj = np.diff(obj)
	delta_obj = np.log(np.abs(delta_obj))
	new_obj = np.diff(delta_obj)
	# The last big change in the objective is our point 
	# of interest. 
	
	# # Note for real data, I need to ignore the last 100
	new_obj = new_obj[:int(len(obj)/5)]
	new_obj = new_obj[::-1]
	# Assuming the dimension is always greater than 5
	i_start=5
	var_val = np.std(new_obj[:i_start])
	mean_val = np.mean(new_obj[:i_start])
	for i in range(i_start,len(new_obj)):
		if new_obj[i] > mean_val - 10*var_val:#4*var_val:
			var_val = np.std(new_obj[:i])
			mean_val = np.mean(new_obj[:i])
		else:
			return len(new_obj) - i -1

	# In case the above method never gives a response:
	# explanation of +2: 
	# Everytime we apply diff, the vector dimension reduces by 1
	return np.argmin(new_obj) + 2

