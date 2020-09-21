function [ A ] = admm_A(Y,A,M,dM,muA,eps_abs,eps_rel,mu,tau_incr,tau_decr,Niter_ADMM)
% Optimization w.r.t. A by ADMM.
%%
% Code : Pierre-Antoine Thouvenin, February 17th 2015.
%%
%-------------------------------------------------------------------------%
% Inputs:
% Y          image pixels (lexicographically ordered)(L|N); 
% A          initial abundance matrix (K|N);
% M          initial endmember matrix (L|K);
% dM         initial perturbation matrices (L|NK)(cell (1|N));
% muA        regularisation parameter (Augmented Lagrangian (AL));
% eps_abs    constant related to the stopping criterion;
% eps_rel      ------------------------------------   ;
% mu, tu_incr, tau_decr   hyperparmeter update values;
% Niter_ADMM   ADMM maximum iteration number.
%
% Output:
% A          abundance matrix (L|K).
%-------------------------------------------------------------------------%
%%
%-Intermediary variables : constraints------------------------------------%
N = size(Y,2);
K = size(M,2);
A1 = [eye(K);ones(1,K)];
B1 = [-eye(K);zeros(1,K)];
c1 = [zeros(K,1);1];

for n = 1:N
    count = 0;
    %-Selection of the associated perturbation matrix
    M1 = M + dM{n};
    
    %-Stopping criterion initialisation
    eps_pri = eps_abs;
    eps_dual = eps_abs;
    norm_pri  = Inf;
    norm_dual = Inf;
    
    %-Lagrange's multiplier initialisation
    lambda = zeros(K+1,1); % [lambda1;lambda2]
    
    u = zeros(K,1);
    u_prev = zeros(K,1);
    
    while ((norm_pri > eps_pri || norm_dual > eps_dual) && count < Niter_ADMM)
        
        lambda1 = lambda(1:K);
%         keyboard
        
        %-A(:,n) update
        A(:,n) = (M1'*M1 + muA(n)*(A1.')*A1)\(M1'*Y(:,n) + muA(n)*(A1.')*(c1 - B1*u - lambda));
        
        %-Splitting variable update
        u = (A(:,n) + lambda1 >= 0).*(A(:,n) + lambda1);  % projection on the constraint space (positivity)
        
        %-Lagrange's mutliplier update
        lambda = lambda + A1*A(:,n) + B1*u - c1;
        
        %-Error computation / count
        norm_pri = norm(A1*A(:,n) + B1*u - c1, 2);
        norm_dual = norm(muA(n)*(A1.')*B1*(u - u_prev), 2);
        
        eps_pri = sqrt(K+1)*eps_abs + eps_rel*max([norm(A1*A(:,n),2), norm(B1*u,2), norm(c1,2)]);
        eps_dual = sqrt(K)*eps_abs + eps_rel*norm(muA(n)*(A1.')*lambda,2); 
        
        muA_prev = muA(n); % previous muA saved before potential update
        
        %-Regularisation parameter update
        if (norm_pri > mu*norm_dual)
            muA(n) = tau_incr*muA(n);
            %-Lagrange's multiplier update if muA updated
            lambda = muA_prev*lambda/muA(n);
        elseif (norm_dual > mu*norm_pri)
            muA(n) = muA(n)*tau_decr;
            %-Lagrange's multiplier update if muA updated
            lambda = muA_prev*lambda/muA(n);
        end
        
        u_prev = u;
        count = count + 1; 
    end
end