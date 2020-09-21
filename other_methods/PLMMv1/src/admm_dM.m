function [ dM ] = admm_dM(Y,A,M,dM,gamma,mu_dM,eps_abs,eps_rel,mu,tau_incr,tau_decr,Niter_ADMM)
% Optimization w.r.t dM by ADMM.
%%
% Code : Pierre-Antoine Thouvenin, February 17th 2015.
%%
%-------------------------------------------------------------------------%
% Inputs:
% Y           image pixels (lexicographically ordered)(L|N);
% A           initial abundance matrix (K|N);
% M           initial endmember matrix (L|K);
% dM          initial perturbation matrices (L|NK)(cell (1|N));
% gamma       regularisation parameter;
% mu_dM       regularisation parameter (Augmented Lagrangian (AL));
% eps_abs     constant related to the stopping criterion;
% eps_rel       ------------------------------------   ;
% mu, tu_incr, tau_decr   hyperparmeter update values;
% Niter_ADMM  ADMM maximum iteration number.
%
% Output:
% dM          endmembers matrix (L|NK)[cell (1|N)].
%-------------------------------------------------------------------------%
%%
[L,N] = size(Y);
K = size(M,2);

if numel(mu_dM) < 2
    mu_dM = mu_dM*ones(1,N);
end

for n = 1:N
    count = 0;
    
    %-Stopping criterion initialization
    eps_pri = eps_abs;
    eps_dual = eps_abs;
    norm_pri  = Inf;
    norm_dual = Inf;
    
    %-Lagrange's multiplier initialization 
    Lambda = zeros(L,K);
    Mn = zeros(L,K);
    Mn_prev = zeros(L,K);
    
    while ((norm_pri > eps_pri) || (norm_dual > eps_dual)) && (count < Niter_ADMM)       
        
%         keyboard
        
        %-dM{n} update
        dM{n} = ((Y(:,n) - M*A(:,n))*(A(:,n)') + mu_dM(n)*(Mn - M - Lambda))/(A(:,n)*(A(:,n)') + (gamma + mu_dM(n))*eye(K));
        
        %-Splitting variable update (projection)
        Mn = dM{n} + M + Lambda;
        id = (Mn < 0);
        Mn(id) = 0;
        
        %-Lagrange's mutliplier update
        Lambda = Lambda + M + dM{n} - Mn;  
        
        %-Error computation
        norm_pri = norm(M + dM{n} - Mn, 'fro');
        norm_dual = norm(-mu_dM(n)*(Mn - Mn_prev), 'fro');
        
        eps_pri = sqrt(L*K)*eps_abs + eps_rel*max([norm(dM{n},'fro'),norm(M,'fro'), norm(Mn,'fro')]);
        eps_dual = sqrt(L*K)*eps_abs + eps_rel*norm(mu_dM(n)*Lambda,'fro');
        
        mu_dM_prev = mu_dM(n); % previous mudM saved before potential update
        
        %-Regularization parameter update
        if (norm_pri > mu*norm_dual)
            mu_dM(n) = tau_incr*mu_dM(n);
            %-Lagrange's multiplier update if mudM updated
            Lambda = mu_dM_prev*Lambda/mu_dM(n);
        elseif (norm_dual > mu*norm_pri)
            mu_dM(n) = mu_dM(n)*tau_decr;
            %-Lagrange's multiplier update if mudM updated
            Lambda = mu_dM_prev*Lambda/mu_dM(n);
        end 
        
        Mn_prev = Mn;
        count = count + 1; 
    end
end