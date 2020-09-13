This toolbox contains several scripts and functions for MATLAB, to unmix hyperspectral data using the Extended Linear Mixing Model (ELMM).

The contents include:

- ELMM_ADMM.m: Function performing the unmixing with the ELMM
- demo_houston.m: example of use of the function on a real hyperspectral dataset
- real_data_1.mat: real dataset used (crop of the DFC 2013 data)
- endmembers_houston.mat: reference endmember matrix used in the demo
- FCLSU.m : function performing the standard fully constrained least squared unmixing.
- CLSU.m : function performing the standard partially constrained least squared unmixing.
- SCLSU.m: function performing a scaled version of CLSU, which follows a particular case of the ELMM
- soft.m: soft thresholding, proximal operator of the L1 norm
- vector_soft_col.m: vector soft thresholding, proximal operator of the L2 norm
- pca_viz.m: function projecting data and endmembers on the space spanned by the first three principal components of the data and displaying a scatterplot
- rescale.m: function rescaling hyperspectral data between 0 and 1

