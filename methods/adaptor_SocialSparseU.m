function [A,M,time_socialSparseU,Yhat] = adaptor_SocialSparseU(Yim,M0,Lib,A_init,opt)


if ~exist('opt','var')
    % tune regularization parameters
    fraction = 1/10;
    lambda   = 0.4;
else
    fraction = opt.fraction;
    lambda   = opt.lambda;
end


[nr,nc,L] = size(Yim);
Y = reshape(Yim,nr*nc,L)';


P = size(M0,2);
bundle = [];
groups = [];
for i=1:P
    bundle = [bundle Lib{i}];
    groups = [groups; i*ones(size(Lib{i},2),1)];
end

X = Y;
A_init = FCLSU(X,bundle)';

cd methods/SocialSparseU


tic

%% initialize params

% A_init = A_FCLSU_bundle; % abundance initialization for other algorithms
rho = 10;
tol_a = 10^(-6);
maxiter_ADMM = 1000;
fraction = 1/10;
verbose = 0; % display


%{
%% Group penalty
lambda = 2;
disp('group')
type = 'group';
[A_group, optim_struct_group] = social_unmixing(X,bundle,groups,A_init,lambda,rho,maxiter_ADMM,type,fraction,tol_a,verbose);
[A_group_final,S_group] = bundle2global(A_group,bundle,groups); % sum the abundances within each class
%}


%{
%% elitist penalty
lambda = 0.5;
disp('elitist')
type = 'elitist';
[A_elitist, optim_struct_elitist] = social_unmixing(X,bundle,groups,A_init,lambda,rho,maxiter_ADMM,type,fraction,tol_a,verbose);
[A_elitist_final,S_elitist] = bundle2global(A_elitist,bundle,groups); % sum the abundances within each class
%}


%% Fractional penalty
%fraction = 1/10;
%lambda = 0.4;
%disp('fractional')
type = 'fractional';
[A_fractional, optim_struct_fractional] = social_unmixing(X,bundle,groups,A_init,lambda,rho,maxiter_ADMM,type,fraction,tol_a,verbose);
[A_fractional_final,S_fractional] = bundle2global(A_fractional,bundle,groups); % sum the abundances within each class
A = A_fractional_final; 
M = S_fractional;


%{
%% Collaborative sparsity
lambda = 1;
disp('collaborative')
type = 'asc';
[A_collaborative] = ADMM_collaborative_unmixing(X,A_init,bundle,lambda,rho,maxiter_ADMM,type,tol_a,verbose);
[A_collaborative_final, S_collaborative] = bundle2global(A_collaborative,bundle,groups); % sum the abundances within each class
%}


time_socialSparseU = toc;

Yhat = zeros(L,nr*nc);
for i=1:nr*nc
    Yhat(:,i) = M(:,:,i) * A(:,i);
end


cd ../..


