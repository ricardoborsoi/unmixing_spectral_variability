This projects contain a Matlab spectral unmixing toolbox related to the paper

> Borsoi, R.A., Imbiriba, T., Bermudez, J.C.M., Richard, C., Chanussot, J., Drumetz, L., Tourneret, J.Y., Zare, A. and Jutten, C.  
> __Spectral Variability in Hyperspectral Data Unmixing: A Comprehensive Review.__  
> arXiv preprint arXiv:2001.07307 (2020).  
> ArXiv link: <https://arxiv.org/pdf/2001.07307>  

It contains a comparative evaluation of some spectral unmixing algorithms which consider endmember variability, using realistic synthetically generated data. See Section VI of the reference above for further details.


## Setup and running

Executing `main.m` or `main_montecarlo.m` on Matlab should generate the synthetic data, process it with the unmixing algorithms, and display the results. The `main.m` code runs a single realization and display the estimated endmembers and abundance maps,`main_montecarlo.m` runs a Monte Carlo simulation and only displays the average quantitative metrics.


## Synthetic data generation

The algorithms are evaluated using synthetic endmember variability data, generated according to a simplification of some radiative transfer models. These were:

* A simplification of Hapke's model, to generate *soil* endmembers according to different viewing angles
* The PROSPECT-D model, to generate *vegetation* endmembers according to different biophysical parameters
* A simplified atmospheric compensation model, to generate *water* endmembers according to different viewing angles

See Section VI-A of the reference above for further details.


## Spectral unmixing methods

The unmixing algorithms are contained in the folder 'methods', and can be easily interfaced with the `main.m` code using adaptors. The following algorithms are included:

* __FCLS__ (Fully Constrained Least Squares): Uses fixed endmembers (i.e., do not considers spectral variability)
* __MESMA__ (Multiple Endmember Spectral Mixture Analysis): Needs spectral libraries, perform an exhaustive search for the endmembers therein which best fit the image pixels
* __Sparse Unmixing with Fractional Norms__: Needs spectral libraries, formulates unmixing as a sparse regression problem using a mixed norm regularization
* __ELMM__ (Extended Linear Mixing Model): Employs a parametric endmember model which assumes that the spectral variability is well represented as the linear scaling of reference spectral signatures
* __DeepGUn__ (Deep Generative Endmember Model): Employs a parametric endmember model represented as a neural network that is learned using pure pixels extracted from the observed image
* __RUSAL__ (Robust Unmixing): Addresses the effects of spectral variability by using an additive residual term in the Linear Mixing Model
* __NCM__ (Normal Compositional Model): Bayesian strategy, represents the endmembers as Gaussian random variables
* __BCM__ (Beta Compositional Model): Bayesian strategy, represents the endmembers as Beta random variables 

For more details on these algorithms, see the full references, which are included below.

##  Library extraction

Spectral libraries containing different instances of endmember spectra are required by some of the unmixing algorithms. They were extracted directly from the observed hyperspectral image using the procedure proposed in

> Somers, B., Zortea, M., Plaza, A. and Asner, G.P.
> __Automated extraction of image-based endmember bundles for improved spectral unmixing.__  
> IEEE Journal of Selected Topics in Applied Earth Observations and Remote Sensing, 5(2), pp.396-408 (2012).


## Contributors

Codes for many of the unmixing algorithms used in this toolbox have been generously provided by the the authors of the original papers. These include:

* __Sparse Unmixing with Fractional Norms__: Provided by Lucas Drumetz and collaborators, and can be found [here]()
* __ELMM__: 
* __RUSAL__: 
* __NCM__: 
* __BCM__:  Xiaoxiao Du and Alina Zare <https://github.com/GatorSense/BetaCompositionalModel>
* __Spectral libraries extraction__:
* The Matlab codes simulating the __PROSPECT-D__ vegetation variability model was provided by Jean-Baptiste Feret and collaborators, and can be found [here](http://teledetection.ipgp.jussieu.fr/prosail/)


## References

1. __FCLS__:  
    > Heinz, D.C.  
    > __Fully constrained least squares linear spectral mixture analysis method for material quantification in hyperspectral imagery.__  
    > IEEE transactions on geoscience and remote sensing, 39(3), pp.529-545 (2001).

2. __MESMA__:  
    > Roberts, D.A., Gardner, M., Church, R., Ustin, S., Scheer, G. and Green, R.O.   
    > __Mapping chaparral in the Santa Monica Mountains using multiple endmember spectral mixture models.__   
    > Remote sensing of environment, 65(3), pp.267-279 (1998).

3. __Sparse Unmixing with Fractional Norms__:  
    > Drumetz, L., Meyer, T.R., Chanussot, J., Bertozzi, A.L. and Jutten, C.   
    > __Hyperspectral image unmixing with endmember bundles and group sparsity inducing mixed norms.__   
    > IEEE Transactions on Image Processing, 28(7), pp.3435-3450 (2019).  

4. __ELMM__:  
    > Drumetz, L., Veganzones, M.A., Henrot, S., Phlypo, R., Chanussot, J. and Jutten, C.  
    > __Blind hyperspectral unmixing using an extended linear mixing model to address spectral variability.__ 
    > IEEE Transactions on Image Processing, 25(8), pp.3890-3905 (2016).

5. __DeepGUn__:  
    > Borsoi, R.A., Imbiriba, T. and Bermudez, J.C.M.  
    > __Deep generative endmember modeling: An application to unsupervised spectral unmixing.__  
    > IEEE Transactions on Computational Imaging, 6, pp.374-384 (2019).

6. __RUSAL__:  
    > Halimi, A., Bioucas-Dias, J.M., Dobigeon, N., Buller, G.S. and McLaughlin, S.   
    > __Fast hyperspectral unmixing in presence of nonlinearity or mismodeling effects.__  
    > IEEE Transactions on Computational Imaging, 3(2), pp.146-159 (2016).

7. __NCM__:  
    > Eches, O., Dobigeon, N., Mailhes, C. and Tourneret, J.Y.   
    > __Bayesian estimation of linear mixtures using the normal compositional model. Application to hyperspectral imagery.__    
    > IEEE Transactions on Image Processing, 19(6), pp.1403-1413 (2010).   

8. __BCM__:  
    > Du, X., Zare, A., Gader, P. and Dranishnikov, D.  
    > __Spatial and spectral unmixing using the beta compositional model.__  
    > IEEE Journal of Selected Topics in Applied Earth Observations and Remote Sensing, 7(6), pp.1994-2003 (2014).



