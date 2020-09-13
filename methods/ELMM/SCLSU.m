function [A,psi_maps] = SCLSU(data, S0)
% ELMM_ADMM Unmix hyperspectral data using the Scaled (Partially)
% Constrained Least Squares Unmixing (SCLSU) approach.
% The solutions follow the Extended Linear Mixing Model (ELMM), with the 
% same scaling factor for all materials
%   Inputs:
%   -data: m*n*L image cube, where m is the number of rows, n the number of
%   columns, and L the number of spectral bands.
%   -S0: L*P reference endmember matrix, with P the number of endmembers
%   to consider
%
%   Outputs:
%   -A: P*N abundance matrix
%   -psi_maps: P*N scaling factor matrix (same entries for each material)
%   The algorithm is presented in:
%
%   L. Drumetz, M. A. Veganzones, S. Henrot, R. Phlypo, J. Chanussot and 
%   C. Jutten, "Blind Hyperspectral Unmixing Using an Extended Linear
%   Mixing Model to Address Spectral Variability," in IEEE Transactions on 
%   Image Processing, vol. 25, no. 8, pp. 3890-3905, Aug. 2016.
%
%   Last Revision: 28-October-2016.
%   Revision: 1.0

%% Recover parameters

[m,n,L] = size(data);
P = size(S0,2); % number of endmembers

%% Compute solution of CLSU

A_CLSU = reshape(CLSU(data, S0),m*n,P)';

%% Recover scaling factors and abundances

CLSU_factor = sum(A_CLSU);
psi_maps = repmat(CLSU_factor,[P,1]);

A = A_CLSU./psi_maps;

end

