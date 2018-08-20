## Select Order
Select the optimal order of latent variables for PCA/ICA/PPCA/PICA/CCA.
This is a public version of the original repository selectOrder, which is private. I've added simple demos that use the functions. Please feel free to contact me for help on getting started. I'm happy to help. 

This tool-box contains order selection algorithms for the linear admixture models ([see this link](https://github.com/idnavid/selectOrder/blob/master/notes/disecting_correlation.pdf)).
Most existing methods are based on information-theoretic criteria, see [this document](https://github.com/idnavid/selectOrder/blob/master/notes/deriving_aic.pdf) and typically use the log-likelihood fucntion. The proposed algorithm used here is not based on the likelihood function, and therefore 
has been shown to produce more consistent results. 

### Citations:
Please make sure to cite the our paper:

*Seghouane, Shokouhi, "Consistent Esitmaiton of Dimensionality for Data Driven Methods in fMRI Analysis", IEEE Transactions on Medical Imaging, 2018.*


### Description
This repository contains code for two types of algorithms: 
- single-vector factor analysis models (e.g., PCA) [`Python`](https://github.com/idnavid/selectOrder_public/tree/master/code/standard/python) and [`Matlab`](https://github.com/idnavid/selectOrder_public/tree/master/code/standard/matlab)
- double-vector analysis (e.g., CCA) [`Matlab`](https://github.com/idnavid/selectOrder_public/tree/master/code/cca/matlab)

### Python Requirements: 
- `numpy`, `scipy`

Additionally, for experiments on real FMRI data: 
- `nipype` (This automatically installs some pre-req modules, e.g., `nibabel`)
- `nilearn`

Note: Installing `nipype` on Windows Anaconda is a bit tricky. 
This is because pip isn't able to directly install `traits`. 
The solution for me was to download the `traits` wheel from from [here](https://www.lfd.uci.edu/~gohlke/pythonlibs/)
and then run 

`pip install <location of trait wheel file>.whl`

### Matlab Requirements:
- `Matlab` version at the time of publishing this code was 2017b

For direction-of-arrival simulations in the CCA directory, I use [doa-tools](https://github.com/morriswmz/doa-tools.git). You can find a forked version [here](https://github.com/idnavid/doa-tools) in case @morriswmz changes the original. 

NS, 2018

