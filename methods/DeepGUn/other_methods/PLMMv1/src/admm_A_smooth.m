function [ A ] = admm_A_smooth(Y,A,M,dM,W,H,alpha,muA,eps_abs,eps_rel,mu,tau_incr,tau_decr,Niter_ADMM)
% Optimization w.r.t A by ADMM.
%%
% Code : Pierre-Antoine Thouvenin, February 17th 2015.
%%
%-------------------------------------------------------------------------%
% Inputs:
% Y         image pixels (lexicographically ordered)(L|N);
% A         initial abundance matrix (K|N);
% M         initial endmember matrix (L|K);
% dM        initial perturbation matrices (L|NK)(cell (1|N));
% W         width of the image; (N = W*H)
% H         heigth of th image;
% alpha     regularisation parameter (spatial piecewise smoothness);
% muA       regularisation parameter (Augmented Lagrangian (AL));
% eps_abs   constant related to the stopping criterion;
% eps_rel      ------------------------------------   ;
% mu        regularisation update parameter (muA);
% tau_incr      -------------------------    (muA);
% tau_decr      -------------------------    (muA);
% mu, tu_incr, tau_decr   hyperparmeter update values;
% Niter_ADMM   ADMM maximum iteration number.
%
% Output:
% A         abundance matrix (L|K).
%-------------------------------------------------------------------------%
%%
N = size(Y,2);
K = size(M,2);
%-Intermediary variables : constraints------------------------------------%
A1 = [eye(K);ones(1,K)];
B1 = [-eye(K);zeros(1,K)];
c1 = [zeros(K,1);1];

%-Intermediary variables
hUpA    = [zeros(1,W),ones(1,N-W)];
hDownA  = [ones(1,N-W),zeros(1,W)];

a = [0,ones(1,W-1)];
a = a(ones(H,1),:);  % line "a" duplicated H times
b = [ones(1,W-1),0];
b = b(ones(H,1),:);  % line "b" duplicated H times
hLeftA  = reshape(a',1,N);
hRightA = reshape(b',1,N);

clear a b;

for n = 1:N
    
    %-Smoothness constraint
    UpA = zeros(K,N);
    A2 = A;
    A2(:,n) = zeros(K,1); % suppression des termes liés à la colonne n
    for k = 1:N-W
        UpA(:,W+k) = A2(:,W+k) - A2(:,k);
    end
    DownA  = [-UpA(:,W+1:end),zeros(K,W)];
    LeftA  = zeros(K,N);
    RightA = zeros(K,N);
    for k = 0:H-1
        LeftA(:,1+k*W:(k+1)*W)  = [zeros(K,1), diff(A2(:,1+k*W:(k+1)*W),1,2)];
        RightA(:,1+k*W:(k+1)*W) = [-LeftA(:,2+k*W:(k+1)*W), zeros(K,1)];
    end

    %-Intermediary variables (smoothness constraint)
    constA = hLeftA(n)^2 +  hRightA(n)^2 +  hDownA(n)^2 +  hUpA(n)^2;
    c = (hLeftA(n)*LeftA(:,n) +  hRightA(n)*RightA(:,n) +  hDownA(n)*DownA(:,n) +  hUpA(n)*UpA(:,n));
       
    count = 0;
    
    %-Selection of the associated perturbation matrix
    M_var = dM{n};
    M1 = M + M_var;
    
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
        A(:,n) = (M1'*M1 + muA(n)*(A1.')*A1 + (alpha*constA)*eye(K))\(M1'*Y(:,n) - alpha*c + muA(n)*(A1.')*(-lambda - B1*u + c1)); % correction
        
        %-Splitting variable update
        u = (A(:,n) + lambda1 >= 0).*(A(:,n) + lambda1); % projection on the constraint space (positivity)
        
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