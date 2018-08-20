# NOTES
The objective of this project is to find the optimal order fMRI data for 
a successful Independent Component Analysis. 
To understand this, we'll need some context, addressing: 
  - What kind of decomposition is ICA used for? 
  - How does the ICA order impact fMRI analysis? 
  - How do we derive the order of intrinsic components?
    - [What are the assumptions and obstacles](https://github.com/idnavid/selectOrder/blob/master/notes/eigenspectrum_analysis.ipynb)?


## ICA for fMRI (from [Oxford FMRIB lectures](https://fsl.fmrib.ox.ac.uk/fslcourse/lectures/ica.pdf))
fMRI data can be explored in two different ways: 
- **Confirmatory**: compare observations with underlying model, Generalized Linear Model, etc.
- **Exploratory**:  find "interesting" patterns in data, PCA/ICA, etc. 

ICA is a model-free paradigm. Because it is model-free, unlike GLM, ICA doesn't require an experiment paradigm (i.e., event signal). This makes ICA a suitable candidate for resting-state fMRI analysis, where there is no event signal. 

>> **Model-free** analyses such as ICA are useful when there is no explicit time-series model of assumed ‘activity’. 


In ICA, the assumption is that fMRI data can be linearly decomposed into two signals: time-coponents, spatial components. 
The voxel-time representation of fMRI data forms a 2D matrix, where each column represents a voxel (i.e., space) and each row 
represents a time sample. ICA assumes that this 2D matrix, *Y*, can be decomposed into a matrix multiplications, *Y = AS*, where *A* represents the time components and *S* the spacial components. 

In using ICA, one can either assume that: 
- spatial components (rows of *S*) are independent (**spatial ICA**). In this case, each column of *Y* is considered a sample. 
- time components (columns of *A*) are independent (**temporal ICA**). Each row of *Y* is a sample. 

Since there are more voxels than there are time samples, we typically perform spatial ICA. The figure below summarizes the ICA decomposition assumption. 

<p>
    <img src="https://github.com/idnavid/selectOrder/blob/master/notes/figures/spatialICA.png" alt>
    <em>Taken from FMRIB lectures, Oxford University, UK.</em>
</p>

In this context, *A* becomes the unmixing matrix and *S* becomes the matrix of independent sources. 

## Underfitting/Overfitting 
The number of columns of *A* (or equivalently, number of rows of *S*) determine the number of independent components. 
Using fewer ICs than necessary results in underfitted data, in which case some independent spatial components are grouped together (this is in the ideal case where ICA is able to converge to something meaningful). 
Overfitting means that the order of ICA is more than necessary. In this case, some ICs are broken into smaller spatial maps. 
Unless the exact number of ICs is available, underfitting or overfitting is unavoidable. This project finds the correct number of ICs in fMRI data. 

