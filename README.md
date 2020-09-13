This projects contain a Matlab software package related to the paper


    Borsoi, R.A., Imbiriba, T., Bermudez, J.C.M., Richard, C., Chanussot, J., Drumetz, L., Tourneret, J.Y., Zare, A. and Jutten, C., 2020.  
    **Spectral Variability in Hyperspectral Data Unmixing: A Comprehensive Review.** arXiv preprint arXiv:2001.07307.  
    ArXiv link: <https://arxiv.org/pdf/2001.07307>

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

* FCLS (Fully Constrained Least Squares): Uses fixed endmembers
* MESMA (Multiple Endmember Spectral Mixture Analysis): Needs spectral libraries, perform an exhaustive search for the endmembers therein which best fit the image pixels
* Sparse Unmixing with Fractional Norms: Needs spectral libraries, formulates unmixing as a sparse regression problem using a mixed norm regularization
* ELMM (Extended Linear Mixing Model): Employs a parametric endmember model which assumes that the spectral variability is well represented as the linear scaling of reference spectral signatures
* DeepGUn (Deep Generative Endmember Model): Employs a parametric endmember model represented as a neural network that is learned using pure pixels extracted from the observed image
* RUSAL (Robust Unmixing): Addresses the effects of spectral variability by using an additive residual term in the Linear Mixing Model
* NCM (Normal Compositional Model): Bayesian strategy, represents the endmembers as Gaussian random variables
* BCM (Beta Compositional Model): Bayesian strategy, represents the endmembers as Beta random variables 

For more details on these algorithms, see the full references, which are included below.

##  Library extraction

Spectral libraries containing different instances of endmember spectra are required by some of the unmixing algorithms. They were extracted directly from the observed hyperspectral image using the procedure proposed in

    Somers, B., Zortea, M., Plaza, A. and Asner, G.P., 2012. **Automated extraction of image-based endmember bundles for improved spectral unmixing.** IEEE Journal of Selected Topics in Applied Earth Observations and Remote Sensing, 5(2), pp.396-408.


## Contributors

Codes for many of the unmixing algorithms used in this toolbox have been generously provided by the the authors of the original papers. These include:

* Sparse Unmixing with Fractional Norms:
* ELMM
* RUSAL
* NCM
* BCM
* Library extraction with batch VCA:

* The Matlab codes simulating the PROSPECT-D vegetation model 


## References
1. FCLS:
asasasa  
asasasa

2. MESMA:
asasas  
asasas

3. Sparse Unmixing with Fractional Norms:
ELMM
DeepGUn
RUSAL
NCM
BCM


