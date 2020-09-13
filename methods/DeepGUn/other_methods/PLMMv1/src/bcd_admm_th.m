function [f,A,M,dM,RE,GMSE_A,GMSE_M,GMSE_dM,SA] = bcd_admm_th(Y,Ath,Mth,dMth,A,M,dM,W,H,gamma,eps_abs,eps_rel,epsilon,varargin)
% Resolution by ADMM of the linear unmixing problem in hyperspectral 
% imagery accounting for spectral and spatial variabilities.
%%
% Code : Pierre-Antoine Thouvenin, February 17th 2015.
%
%% Model
%-------------------------------------------------------------------------%
%   Y   =  M    A   + [dM1*a1|...|dMn*an] +   B     (PLMM)
% (L|N)  (L|K)(K|N)                         (L|N)
%
% 0.5*|Y - MA - dMA|² + 0.5*alpha*\Phi(A) + 0.5*beta*\Psi(M) + 0.5*gamma*|dM|²
% 
% s.t.            A >= 0         with N : number of pixels;    
%                 M >= 0              K : number of endmembers;
% forall n  M + dMn >= 0              L : number of spectral bands.
%             A'   1 =  1
%           (N|K)(K|1)(N|1)
%-------------------------------------------------------------------------%
%%
% Inputs:
% >  Y     lexicographically ordered image pixels (L|N);
% >  Ath   true abundance matrix (K|N);
% >  Mth   true endmembers matrix (L|K);
% > dMth   true perturbation matrices (L|NK)[required format : cell (1|N)];
% >  A     initial abundance matrix (K|N);
% >  M     initial endmembers matrix (L|K);
% > dM     initial perturbation matrices (L|NK)[required format : cell (1|N)];
% >  W     width of the image; (N = W*H)
% >  H     heigth of the image;
% > sigma2 constraint on the perturbation matrices;
%
% Optional parameters
% > muA, muM, mu_dM         regularisation parameters (Augmented Lagrangian (AL))
% > mu, tau_incr, tau_decr  hyperparameter update values;
% > Niter_ADMM              ADMM maximum iteration number.
% > alpha                   abundance penalty parameter (spatial smoothness);
% > beta                    endmember penalty parameter ('NONE','DISTANCE','MUTUAL DISTANCE','VOLUME')
%
% >> DISTANCE
% --> M0   reference signatures
%
% >> VOLUME
% --> Y_proj  projected data(PCA)(K-1|N)
% --> Y_bar   data mean (mean(Y,2))(L|1)
% --> U       inverse projector (PCA variables -> M)(L|K-1);
% --> V       projector (M -> PCA variables)(K-1|L);
%
% Outputs:
% < f   objective function value at each BCD iteration;
% < A   abundance matrix (K|N);
% < M   endmember matrix (L|PK;
% < dM  perturbation matrices (K|NP) [contained in a (1|N) cell];
% < mu_A   hyperparamters final values (abundance sub-problems);
% < mu_M   --------------------------- (endmember sub-problems); 
% < mu_dM  --------------------------- (variability sub-problems);
% < eps_abs,eps_rel   subproblems stopping criterion;
% < epsilon   global stopping criterion;
%
% Example:
% [f,A,M,dM] = bcd_admm(Y,Ath,Mth,dMth,A,M,dM,W,H,gamma,eps_abs,eps_rel,epsilon,'HYPERPARAMETERS',{muA,muM,mudM},'PENALTY A',alpha,'PENALTY M',{type,beta,Y_proj,Y_bar,U,V},'AL INCREMENT',{tau_incr,tau_decr,mu},'MAX ADMM STEPS',Niter_ADMM);
%-------------------------------------------------------------------------%
%%
% Number of inputs must be >=minargs and <=maxargs.
narginchk(13, 23);

%--------------------------------------------------------------
% Default parameter values
%--------------------------------------------------------------
% Default Augmented Lagrangian's parameter (hyperparameters)
muA = 1e-4;
muM = 1e-4;
mudM = 1e-4; 
% Hyperparameters' increment
tau_incr = 1.1;
tau_decr = 1.1;
mu = 10;
% Max ADMM steps
Niter_ADMM = 50;
% Default penalty parameters
alpha = 0;
beta = 0;
% Default endmember penalty
type = 'NONE';
% Initial error
err = Inf;

%--------------------------------------------------------------
% Optional parameters
%--------------------------------------------------------------
if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'HYPERPARAMETERS'
                muA = varargin{i+1}{1};
                muM = varargin{i+1}{2};
                mudM = varargin{i+1}{3};
            case 'AL INCREMENT'
                tau_incr = varargin{i+1}{1};
                tau_decr = varargin{i+1}{2};
                mu = varargin{i+1}{3};
            case 'MAX ADMM STEPS'
                Niter_ADMM = varargin{i+1};
            case 'PENALTY M'
                type = varargin{i+1}{1};
                beta = varargin{i+1}{2};
                switch upper(type)
                    case 'DISTANCE'
                        M0 = varargin{i+1}{3};
                        input = {type,beta,M0};
                    case 'VOLUME'
                        Y_proj = varargin{i+1}{3};
                        Y_bar = varargin{i+1}{4};
                        U = varargin{i+1}{5};
                        V = varargin{i+1}{6};
                        input = {type,beta,Y_proj,Y_bar,U,V};
                    otherwise
                        input = {type,beta};
                end
            case 'PENALTY A'
                alpha = varargin{i+1};
            otherwise
                error(['Unrecognized option: ''' varargin{i} '''']);
        end
    end
end

%--------------------------------------------------------------
% Algorithm
%--------------------------------------------------------------
f(1) = objective_p(Y,M,A,dM,H,W,alpha,gamma,'PENALTY',input);
[RE(1),GMSE_A(1),GMSE_M(1),~,GMSE_dM(1),~,SA(1,:)] = th_error(Y,Ath,Mth,dMth,A,M,dM,W,H);
k = 2;

while (err > epsilon)

    % Update A
    A = g_admm_A(Y,A,M,dM,W,H,alpha,muA,eps_abs,eps_rel,mu,tau_incr,tau_decr,Niter_ADMM);
    % Update M
    M = g_admm_M(Y,A,M,dM,beta,muM,eps_abs,eps_rel,mu,tau_incr,tau_decr,Niter_ADMM,'PENALTY',input);
    % Update dM
    dM = admm_dM(Y,A,M,dM,gamma,mudM,eps_abs,eps_rel,mu,tau_incr,tau_decr,Niter_ADMM);

     % Objective function value
    f(k) = objective_p(Y,M,A,dM,H,W,alpha,gamma,'PENALTY',input);
    [RE(k),GMSE_A(k),GMSE_M(k),~,GMSE_dM(k),~,SA(k,:)] = th_error(Y,Ath,Mth,dMth,A,M,dM,W,H);
    
    % Stopping criterion    
    err = abs(f(k-1)-f(k))/f(k-1);
    k = k+1;
end