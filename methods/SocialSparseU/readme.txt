This toolbox contains code associated to the paper 

Drumetz, L., Meyer, T. R., Chanussot, J., Bertozzi, A. L., & Jutten, C. (2019). Hyperspectral image unmixing with endmember bundles and group sparsity inducing mixed norms. IEEE Transactions on Image Processing, 28(7), 3435-3450.

This code allows to run any of the three algorithms proposed in the paper as well as other tested methods, and to reproduce the results in the paper concerning the Houston Dataset.

The contents include the following files:

- houston_paper.m: example of use of the function on a real hyperspectral dataset
- real_data_1.mat: real dataset used (crop of the DFC 2013 data)
- EIA_VCA: VCA endmember extraction algorithm
- batchvca: function that extract bundles using the Automated Endmember Bundles (AEB) method and clusters the signatures into groups using kmeans.
- bundle2global: function summing abundance maps corresponding to the same group and computing pixelwise endmembers from the global and variant abundances
- social_unmixing.m: function gathering all the algorithms proposed in the paper: group, elitist and fractional.
- FCLSU.m : function performing the standard fully constrained least squared unmixing.
- ADMM_collaborative_unmixing.m : function performing the unmixing using the collaborative 2,1 norm to eliminate candidate endmembers from all the support of the image jointly.
- prox_elitist_group : proximal operator of the l1,2 group norm used in the elitist penalty
- prox_group_lasso : .m: proximal operator of the l2,1 group norm used in the group penalty
- approx_prox_fractional: approximative shrinkage operator for the l1,p/q quasinorm used in the fractional penalty
- vector_soft_col.m: vector soft thresholding, proximal operator of the L2 norm
- pca_viz.m: function projecting data and endmembers on the space spanned by the first three principal components of the data and displaying a scatterplot
- pca_viz_global.m: function projecting data and endmembers on the space spanned by the first three principal components of the data and displaying a scatterplot. This function displays all the results of all algorithms on the same figure.
- rescale.m: function rescaling hyperspectral data between 0 and 1
- subaxis.m: alternative to the "subplot" function allowing to change spacing between figures
- parseArgs.m: subroutine used by subaxis.m

