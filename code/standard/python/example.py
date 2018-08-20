"""
Simulated an fmri signal by creating two low-rank matrices X_tc, X_sc (for temporal and spatial components). 
The observation signal Y = X_tc*X_sc is the signal we intend to decompose back into X_tc and X_sc. 
For this experiment, we're only interested in the hidden rank of Y. 
"""

from tools.gen_data import *
from tools.estimate_dim import *

if __name__=='__main__':
	n = 100
	p = 30
	q = 5 # number of components
	X_sc = spatial_comp(n,q,'random')
	X_tc = temporal_comp(p,q)


	# NOTE:	
	#   n is the squar-root of the actual number of samples, because
	#   the spatial component X_sc is assumed to be defined in 2D. 
	#	Just like a cross-section of the brain.  
	n_samples = n**2
	noise_sig = .01*np.random.randn(p,n_samples)
	Y = 1.*mix_comp(X_tc,X_sc) + noise_sig
	

	R_y = np.dot(Y,Y.T)/(n_samples)
	Y_eigvals = np.abs(np.linalg.eigvals(R_y))
	Y_eigvals = np.sort(Y_eigvals)[::-1]
	Y_eigvals = adjust_eigenspect(Y_eigvals,p,n_samples)
	
	objective = order_objective_function(Y_eigvals,p,n_samples)
	obj_sum = objective['alt-2']
	q_hat = np.argmin(obj_sum[0])
	q_hat+=1 # NOTE: Numpy indexes start at zero 
	pylab.plot(1+np.arange(len(obj_sum[0])),obj_sum[0],label='Proposed IC')
	pylab.plot(q_hat,np.min(obj_sum[0]),'r*')
	pylab.text(q_hat,np.min(obj_sum[0])*-3,'estimated: '+str(q_hat))
	pylab.text(q_hat,np.min(obj_sum[0])*-2, 'true     : '+str(q))
	pylab.legend()
	pylab.show()
