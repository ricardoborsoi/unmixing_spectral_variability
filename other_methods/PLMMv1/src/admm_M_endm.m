function [ M ] = admm_M_endm(Y,A,M,dM,beta,muM,eps_abs,eps_rel,mu,tau_incr,tau_decr,Niter_ADMM)
% Optimization w.r.t M by ADMM. (mutual distance penalization)
%%
% Code : Pierre-Antoine Thouvenin, February 17th 2015.
%%
%-------------------------------------------------------------------------%
% Inputs:
% Y         image pixels (lexicographically ordered)(L|N); 
% A         initial abundance matrix (P|N);
% M         initial endmember matrix (L|P);
% dM        initial perturbation matrices (L|NP)(cell (1|N));
% beta      regularisation parameter;
% muM       regularisation parameter (Augmented Lagrangian (AL));
% eps_abs   constant related to the stopping criterion;
% eps_rel      ------------------------------------   ;
% mu, tu_incr, tau_decr    hyperparmeter update values;
% Niter_ADMM               ADMM maximum iteration number.
%
% Output:
% M         endmembers matrix (L|P)(cell (1|N)).
%-------------------------------------------------------------------------%
%%
[L,N] = size(Y);
P = size(M,2);

% Computation of the constraint-related matrices (distance between
% endmembers)
R = zeros(P,P*P);
for l = 0:P-1
    el = zeros(P,1);
    el(l+1) = 1;
    R(:,1+l*P:(l+1)*P) = -eye(P) + el*ones(1,P);
end

A1 = repmat(eye(P),[1,N+1]);

C1 = zeros(L,(N+1)*P);
D = zeros(L,N);
for n = 1:N
    D(:,n) = dM{n}*A(:,n);
    C1(:,n*P+1:(n+1)*P) = dM{n};
end

count = 0;

%-Stopping criterion initialisation
eps_pri = eps_abs;
eps_dual = eps_abs;
norm_pri  = Inf;
norm_dual = Inf;

%-Lagrange's multiplier initialisation 
Lambda = zeros(L,(N+1)*P);
V = zeros(L,(N+1)*P);
V_prev = zeros(L,(N+1)*P);

while ((norm_pri > eps_pri) || (norm_dual > eps_dual)) && (count < Niter_ADMM)

    %-M update
    M = ((Y - D)*A' + muM*(V-C1-Lambda)*(A1.'))/(A*(A') + beta*R*(R.') + muM*(A1*(A1.')));
    %-Splitting variable update (projection)
    V = M*A1 + C1 + Lambda;
    id = (V < 0);
    V(id) = 0;

    %-Lagrange's mutliplier update
    Lambda = Lambda + M*A1 - V + C1;  

    %-Error computation / count
    norm_pri = norm(M*A1 - V + C1, 'fro');
    norm_dual = norm(-muM*(V - V_prev)*(A1.'), 'fro');

    eps_pri = sqrt(L*(N+1)*P)*eps_abs + eps_rel*max([norm(M*A1,'fro'), norm(V,'fro'), norm(C1,'fro')]); 
    eps_dual = sqrt(L*P)*eps_abs + eps_rel*norm(muM*Lambda*(A1'),'fro');
    
    muM_prev = muM; % previous muM saved before potential update

    %-Regularisation parameter update
    if (norm_pri > mu*norm_dual)
        muM = tau_incr*muM;
        %-Lagrange's multiplier update if muM updated
        Lambda = muM_prev*Lambda/muM;
    elseif (norm_dual > mu*norm_pri)
        muM = muM*tau_decr;
        %-Lagrange's multiplier update if muM updated
        Lambda = muM_prev*Lambda/muM;
    end 

    V_prev = V;
    count = count + 1; 
end