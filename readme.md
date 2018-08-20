## Select Order
Select the optimal order of latent variables for PCA/ICA/PPCA/PICA/CCA.

This tool-box contains a order selection algorithms for the linear admixture models ([see this link](https://github.com/idnavid/selectOrder/blob/master/notes/disecting_correlation.pdf)).
Most existing methods (based on information-theoretic criteria, see [this document](https://github.com/idnavid/selectOrder/blob/master/notes/deriving_aic.pdf)) rely on the assumption that samples are 
independent and identitically distributed. The proposed algorithm used here is not based on the likelihood function, and therefore 
has been shown to produce more consistent results. 

For more detail, please see our paper in IEEE TMI:<br/>
*Seghouane, Shokouhi, "Consistent Esitmaiton of Dimensionality for Data Driven Methods in fMRI Analysis", IEEE Transactions on Medical Imaging, 2018.*

### Citations:
Please make sure to cite the our paper. 


### Description
This repository contains code for two types of algorithms: 
- single-vector factor analysis (e.g., PCA) [`Python`](here!) and [`Matlab`](here!)
- double-vector analysis (e.g., CCA) [`Matlab`](here!)

### Python Requirements: 
- `numpy`, `scipy`

For experiments on real FMRI data: 
- `nipype` (This automatically installs some pre-req modules, e.g., `nibabel`)
- `nilearn`

Note: Installing `nipype` on Windows Anaconda is a bit tricky. 
This is because pip isn't able to directly install `traits`. 
The solution for me was to download the `traits` wheel from from [here](https://www.lfd.uci.edu/~gohlke/pythonlibs/)
and then run 

`pip install <location of trait wheel file>.whl`

### Matlab Requirements:
- `Matlab`

For direction-of-arrival simulations, I use [doa-tools](https://github.com/morriswmz/doa-tools.git). You can find a forked version [here](here!) in case @morriswmz changes the original. 

NS, 2018

