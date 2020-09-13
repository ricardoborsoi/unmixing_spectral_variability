function [A,M,alg_time,Yhat] = adaptor_ELMM(Yim,M0,Lib,A_init,opt)


if ~exist('opt','var')
    % tune regularization parameters
    lambda_s   = 0.5; 
    lambda_a   = 0.015;
    lambda_psi = 0.05; 
else
    lambda_s   = opt.lambda_s;
    lambda_a   = opt.lambda_a;
    lambda_psi = opt.lambda_psi;
end


[nr,nc,L] = size(Yim);
Y = reshape(Yim,nr*nc,L)';

cd methods/ELMM

psis_init = ones(size(A_init));



% optional parameters

norm = '1,1'; % Use a Total Variation on the abundances
verbose = false; %true; % display

maxiter_anls = 100;
maxiter_admm = 100;
epsilon_s = 10^(-3);
epsilon_a = 10^(-3);
epsilon_psi = 10^(-3);
epsilon_admm_abs = 10^(-2);
epsilon_admm_rel = 10^(-2);

% run

% simple version with default parameters for the stopping criteria

% tic
% [A_ELMM, psis_ELMM, S_ELMM, optim_struct] = ELMM_ADMM(data, A_init, psis_init, S0,lambda_s,lambda_a,lambda_psi,norm,verbose);
% toc

% complete syntax

tic
[A_ELMM, psis_ELMM, S_ELMM, optim_struct] = ELMM_ADMM(Yim, A_init, psis_init, M0,lambda_s,lambda_a,lambda_psi,norm,verbose,maxiter_anls,maxiter_admm,epsilon_s,epsilon_a,epsilon_psi,epsilon_admm_abs,epsilon_admm_rel);
alg_time = toc;

A = A_ELMM;
M = S_ELMM;

Yhat = zeros(L,nr*nc);
for i=1:nr*nc
    Yhat(:,i) = M(:,:,i) * A(:,i);
end


cd ../..




