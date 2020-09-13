function [ M ] = g_admm_M(Y,A,M,dM,beta,muM,eps_abs,eps_rel,mu,tau_incr,tau_decr,Niter_ADMM,varargin)
%-Resolution w.r.t. M by ADMM.
%%
% Code : Pierre-Antoine Thouvenin, February 17th 2015.
%%
%-------------------------------------------------------------------------%
% Inputs:
% Y         pixels (hyperspectral data cube)(L|N)(lexicographically ordered); 
% A         initial abundance matrix (K|N);
% M         initial endmembers matrix (L|K);
% dM        initial perturbation matrices (L|NK)(cell (1|N));
% beta      regularisation parameter;
% muM       regularisation parameter (Augmented Lagrangian (AL));
% eps_abs   constant related to the stopping criterion;
% eps_rel      ------------------------------------   ;
% mu, tu_incr, tau_decr   hyperparmeter update values;
% Niter_ADMM              ADMM maximum iteration number.
%
% Outputs:
% M    endmembers matrix (L|K)(cell (1|N)).
%
% Optional parameters:               
% >> PENALIZATION TYPE
% --> 'NONE','DISTANCE','MUTUAL DISTANCE','VOLUME'  
%
% >> DISTANCE
% --> M0   reference signatures
%
% >> VOLUME
% --> Y_proj  projected data(PCA)(K-1|N)
% --> Y_bar   data mean (mean(Y,2))(L|1)
% --> U       inverse projector (PCA variables -> M)(L|K-1);
% --> V       projector (M -> PCA variables)(K-1|L);
%-------------------------------------------------------------------------%
%%
%--------------------------------------------------------------
% Local variables
%--------------------------------------------------------------
% Size
[K,N] = size(A);

%--------------------------------------------------------------
% Optional parameters
%--------------------------------------------------------------
if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'PENALTY'
                type = varargin{i+1}{1};
                beta = varargin{i+1}{2};
                switch upper(type)
                    case 'DISTANCE'
                        M0 = varargin{i+1}{3};    
                    case 'VOLUME'
                        Y_proj = varargin{i+1}{3};
                        Y_bar = varargin{i+1}{4};
                        U = varargin{i+1}{5};
                        V = varargin{i+1}{6};
                end
            otherwise
                type = 'NONE';
        end
    end
end

switch type
    case 'NONE'
        M = admm_M(Y,A,M,dM,muM,eps_abs,eps_rel,mu,tau_incr,tau_decr,Niter_ADMM);
        
    case 'DISTANCE'
        M = admm_M_d(Y,A,M,M0,dM,beta,muM,eps_abs,eps_rel,mu,tau_incr,tau_decr,Niter_ADMM);
        
    case 'MUTUAL DISTANCE'
        M = admm_M_endm(Y,A,M,dM,beta,muM,eps_abs,eps_rel,mu,tau_incr,tau_decr,Niter_ADMM);
        
    case 'VOLUME'  % row-wise optimization in this case
        if numel(muM) < 2
            muM = muM*ones(1,K-1);
        end
        dT = cell(1,N);
        for n = 1:N
            dT{n} = V*(dM{n} - Y_bar(:,ones(K,1)));
        end
        T = V*(M - Y_bar(:,ones(K,1)));
        T = admm_T_v(Y_proj,Y_bar,A,T,dT,U,V,beta,muM,eps_abs,eps_rel,mu,tau_incr,tau_decr,Niter_ADMM);
        M = U*T + Y_bar(:,ones(K,1)); 
        
    otherwise
        error(['Unrecognized penalty type: ''' type '''']);
end
end
