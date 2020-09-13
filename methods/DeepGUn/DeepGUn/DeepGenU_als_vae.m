function [A,EMs_DeepGen]=DeepGenU_als_vae(data, A_init, mappingAEC, dimAEC, Zref, lambda_zref, lambda_a, flag_useparfor)
% -------------------------------------------------------------------------
% Perform the block coordinate descent optimization routine to estimate the
% endmembers latent representation and the abundance maps.
% 
% 
% INPUTS:
%   data           : nr * nc * L  hyperspectral image cube
%   A_init         : nr * nc * P  abundance maps initialization
%   mappingAEC     : Cell array with generative model of each endmember class
%   dimAEC         : Dimension of the VAEs latent space (recomended = 2)
%   Zref           : Latent representations of the reference endmember
%                    matrix M0
%   lambda_zref    : Regularization parameter
%   lambda_a       : Regularization parameter
%   flag_useparfor : Allows us to disable parfor if necesary
% 
% OUTPUTS:
%   A_deepGen      : Estimated abundances
%   M_DeepGen      : Estimated endmember tensor
%
% -------------------------------------------------------------------------
% This code is part of the DeepGUn algorithm, referring to the following publication:
% 
%   "Deep Generative Endmember Modeling: An Application to Unsupervised Spectral Unmixing"
%   Ricardo Augusto Borsoi, Tales Imbiriba, José Carlos Moreira Bermudez
%   IEEE Transactions on Computational Imaging
% 
% The DeepGUn algorithm performs spectral unmixing with spectral variability
% modeling the endmembers using deep generative models (variational
% autoencoders). 
% -------------------------------------------------------------------------
% Author: Ricardo Borsoi
% last revision: 02/09/2019
%
% -------------------------------------------------------------------------


% use all available cores unless told otherwise
if nargin < 8
    flag_useparfor = inf;
end

% Number of iterations and tolerance
maxiter_als = 8;
eps_tol_als = 1e-3;


% get constants 
[nr,nc,L] = size(data);
[P,N] = size(A_init);


% Parameters for the abundance estimation problem
nnorm = '1,1'; % Use a Total Variation on the abundances
verbose = false; % display
maxiter_anls_glmm = 50;
maxiter_admm = 100;
epsilon_s = 10^(-3);
epsilon_a = 10^(-3);
epsilon_psi = 10^(-3);
epsilon_admm_abs = 10^(-2);
epsilon_admm_rel = 10^(-2);


% initialization --------------------
A        = A_init;
Z_latent = zeros(dimAEC,P,N);

for i=1:maxiter_als
    % fprintf('\n\n ALS iteration number %d \n\n',i)
    
    % "old" variables used to measure the stopping criterion
    Z_latent_old = Z_latent;
    A_old = A;
    
    % Lets estimate the endmembers variables in the latent space
    [EMs_tensor,Z_latent]=DeepGenU_optimization(Z_latent,data,A,mappingAEC,Zref,lambda_zref,flag_useparfor);
    
    % Now lets estimate the abundances
    A = GLMM_fixedEMs(data, A, EMs_tensor,lambda_a,nnorm,verbose,maxiter_anls_glmm,maxiter_admm,epsilon_s,epsilon_a,epsilon_psi,epsilon_admm_abs,epsilon_admm_rel,flag_useparfor);            
    
    
    % stopping criterion ------------------------------------
    if norm(Z_latent(:)-Z_latent_old(:))/norm(Z_latent_old(:)) < eps_tol_als
        if norm(A(:)-A_old(:))/norm(A_old(:)) < eps_tol_als
            break;
        end
    end
end


% Reshape EM tensor -----------------------------------------
EMs_DeepGen = zeros(L,P,N);
for i=1:P
    EMs_DeepGen(:,i,:) = permute(reshape(EMs_tensor{i},nr*nc,L)', [1 3 2]);
end


% Recompute abundances one final time with more iterations ----------------
epsilon_admm_abs = 1e-4;
epsilon_admm_rel = 1e-4;
A = GLMM_fixedEMs(data, A, EMs_tensor,lambda_a,nnorm,verbose,maxiter_anls_glmm,maxiter_admm,epsilon_s,epsilon_a,epsilon_psi,epsilon_admm_abs,epsilon_admm_rel,flag_useparfor);            



