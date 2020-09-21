function [ A ] = g_admm_A(Y,A,M,dM,W,H,alpha,muA,eps_abs,eps_rel,mu,tau_incr,tau_decr,Niter_ADMM)
% Generic abundance update using ADMM steps.
%%
% Code : Pierre-Antoine Thouvenin, February 17th 2015.
%%
%-------------------------------------------------------------------------%
% Inputs:
% > Y        pixels (hyperspectral data cube)(L|N)(lexicographically ordered); 
% > A        initial abundance matrix (K|N);
% > M        initial endmembers matrix (L|K);
% > dM       initial perturbation matrices (L|NK)(cell (1|N));
% > W        width of the image; (N = W*H)
% > H        heigth of th image;
% > alpha    regularisation parameter (spatial piecewise smoothness);
% > muA      regularisation parameter (Augmented Lagrangian (AL));
% > eps_abs  constant related to the stopping criterion;
% > eps_rel     ------------------------------------   ;
% > mu, tau_incr, tau_decr    hyperparameter update values;
% > Niter_ADMM                ADMM maximum iteration number.
%
%Output:
% < A   abundance matrix (L|K).
%-------------------------------------------------------------------------%
%%
N = size(Y,2);

if (numel(muA) < 2)
    muA = muA*ones(1,N);
end

if (alpha == 0)
        A = admm_A(Y,A,M,dM,muA,eps_abs,eps_rel,mu,tau_incr,tau_decr,Niter_ADMM);
else
        A = admm_A_smooth(Y,A,M,dM,W,H,alpha,muA,eps_abs,eps_rel,mu,tau_incr,tau_decr,Niter_ADMM);
end

