function [A_MESMA,M_MESMA,alg_time,Yhat] = adaptor_MESMA(Yim,M0,Lib,A_init,opt)

if ~exist('opt','var')
    % no parameters to tune
end


[nr,nc,L] = size(Yim);
Y = reshape(Yim,nr*nc,L)';

[L,N]=size(Y);
P = size(M0,2);

addpath methods/MESMA

tic
[A_MESMA,r_err,bestIdx] = MESMA(Y, Lib);

M_MESMA = zeros(L,P,N);
for i=1:N
    for k=1:P
        M_MESMA(:,k,i) = Lib{k}(:,bestIdx(i,k));
    end
end

alg_time = toc;


Yhat = zeros(L,N);
for i=1:N
    Yhat(:,i) = M_MESMA(:,:,i)*A_MESMA(:,i);
end
