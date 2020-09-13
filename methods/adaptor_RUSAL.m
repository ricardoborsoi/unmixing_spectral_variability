function [A,M,alg_time,Yhat] = adaptor_RUSAL(Yim,M0,Lib,A_init,opt)

if ~exist('opt','var')
    % no parameters to tune
    tau  = 1;
    tau2 = 1;
else
    tau = opt.tau;
    tau2 = opt.tau2;
end

[nr,nc,L] = size(Yim);
Y = reshape(Yim,nr*nc,L)';
[L,N] = size(Y);

cd methods/RUSAL

tic

AL_iters = 2000; tol  = 10^(-8);
Model = 'RUSAL';
MPlus = M0;
R = size(M0,2);

[zt,stop_RUSAL,iter_RUSAL] = NUSAL_RUSAL_v1(Yim,MPlus,tau,tau2,AL_iters,tol,Model,A_init);

alg_time   = toc;

alpha_RUSAL    = zt(1:R,:,end);
resid_RUSAL    = idct([zt(R+1:end,:,end); zeros(L-20,nr*nc)]);% LxN
%[RMSE(algo) RE(algo) SAM(algo) SAMn(algo,:) RCest(:,:,algo) RE_RC(algo) SAM_RC(algo) SAMn_RC(algo,:)   ...
%    RMSE_Class(algo,:) RE_Class(algo,:) SAM_Class(algo,:) RE_Class_RC(algo,:) SAM_Class_RC(algo,:) REn(algo,:)]  = ...
%    exploitation_TIP(Y_bloc,MPlus,alpha_RUSAL,alpha0,ones(row*col,1),c_Ev0,resid_RUSAL,RComp,K,Label,3); %%%

cd ../..

A = alpha_RUSAL;
M = [];

c_Ev  = ones(nr*nc,1);
RCest = resid_RUSAL;
Yhat  = (M0*alpha_RUSAL) .*(c_Ev*ones(1,L))' + RCest;    


alg_time = toc;
