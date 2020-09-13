function [A,M,alg_time,Yhat] = adaptor_NCM(Yim,M0,Lib,A_init,opt)


if ~exist('opt','var')
    % tune regularization parameters
    %x   = 0; 
else
    %x   = opt.x;
end


[nr,nc,L] = size(Yim);
Y = reshape(Yim,nr*nc,L)';

cd methods/NCM



tic

% Markov chain length
Nmc = 10000;

% number of burn-in period
Nbi = 1000;

% Unmixing procedure

R = size(M0,2);
N = nr*nc;
AlphasMean = zeros(R,N);

for i=1:N
    y = Y(:,i);
    [TalphaPlus,Tsigma2r] = unmixing(y,M0,Nmc,Nbi);
    AlphasMean(:,i) = mean(TalphaPlus(:,Nbi+1:end),2);
end

alg_time = toc;

A = AlphasMean;
M = [];

Yhat = M0*A;

cd ../..




