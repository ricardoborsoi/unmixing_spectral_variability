function [A,M,alg_time,Yhat] = adaptor_FCLS(Yim,M0,Lib,A_init,opt)

if ~exist('opt','var')
    % no parameters to tune
end

[nr,nc,L] = size(Yim);
Y = reshape(Yim,nr*nc,L)';

[L,N] = size(Y);

tic
A = FCLSU(Y,M0)';

M = []; %repmat(M0,[1,1,N]);
alg_time = toc;

Yhat = M0*A;
