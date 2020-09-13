function [A,M,alg_time,Yhat] = adaptor_DeepGUn(Yim,M0,Lib,A_init,opt)


if ~exist('opt','var')
    % tune regularization parameters
    dimAut = 2;
    lambda_zref = 0.1;
    lambda_a = 0.01;
else
    dimAut = opt.dimAut;
    lambda_zref = opt.lambda_zref;
    lambda_a = opt.lambda_a;
end


[nr,nc,L] = size(Yim);
Y = reshape(Yim,nr*nc,L)';

cd methods/DeepGUn

addpath DeepGUn
addpath(genpath('utils'))

%psis_init = ones(size(A_init));



% optional parameters

% parameters of the algorithm
%dimAut = 2;
%lambda_zref = 0.1;
%lambda_a = 0.01;

% parameters for the endmember bundle extraction
P = size(M0,2);
flag_Npx = true;
vec_Npx = 100*ones(1,P); %[100 100 100 100 100];
flag_useparfor = 0;

tic
[A_deepGen,M_DeepGen]=DeepGUn(Yim, A_init, M0, dimAut, lambda_zref, lambda_a, flag_Npx, vec_Npx, flag_useparfor);
alg_time = toc;



A = A_deepGen;
M = M_DeepGen;

Yhat = zeros(L,nr*nc);
for i=1:nr*nc
    Yhat(:,i) = M(:,:,i) * A(:,i);
end

cd ../..




