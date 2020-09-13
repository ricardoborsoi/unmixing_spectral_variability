function [ T ] = admm_T_v(Y_proj,Y_bar,A,T,dT,U,V,beta,muM,eps_abs,eps_rel,mu,tau_incr,tau_decr,Niter_ADMM)
% Optimization w.r.t T by ADMM. (volume penalization)
%%
% Code : Pierre-Antoine Thouvenin, February 17th 2015.
%%
%-------------------------------------------------------------------------%
% Inputs:
% Y_proj  projected data(PCA)(lexicographically ordered)(K-1|N);
% Y_bar   true data mean (mean(Y,2))(L|1);
% A       initial abundance matrix (K|N);
% T       initial projected endmember matrix(PCA)(K-1|K);
% dT      initial projected endmember perturbation matrices(PCA)(K-1|K);
% U       inverse projection matrix (PCA variables -> M)(L|K-1);
% V       projection matrix (M -> PCA variables)(K-1|L);
% beta    regularisation parameter;
% muM     regularisation parameter (Augmented Lagrangian (AL));
% eps_abs constant related to the stopping criterion;
% eps_rel    ------------------------------------   ;
% mu, tu_incr, tau_decr   hyperparmeter update values;
% Niter_ADMM              ADMM maximum iteration number.
%
% Output:
% T       projected endmember matrix (K-1|K).
%-------------------------------------------------------------------------%
%%
%-Intermediary variables
N = size(Y_proj,2);
K = size(T,2);

delta = zeros(K-1,N);
for n = 1:N
    delta(:,n) = dT{n}*A(:,n);
end
S = delta - (U.')*(Y_bar(:,ones(N,1)) - 2*Y_bar(:,ones(K,1))*A);
Z = V*Y_bar(:,ones(K,1));

%-Intermediary variable (volume constraint)
d = zeros(K,1);
t_minus = zeros(1,K);
t_plus = zeros(1,K);
t_minus_N = zeros(N,K);
t_plus_N  = zeros(N,K);
var = zeros(N,K); % variable related to g(T)

%-Unecessary constraint detection
empty_set_plus = find(~any(U > 0)); % index of empty positive sets
empty_set_minus = find(~any(U < 0));% index of empty negative sets

for l = 1:K-1
    
    count = 0;
    
    %-Stopping criterion initialisation
    eps_pri = eps_abs;
    eps_dual = eps_abs;
    norm_pri  = Inf;
    norm_dual = Inf;
        
    %-Positivity constraint on M expressed for T
    id_diff = setdiff(1:K-1,l);
    set_Up = (U(:,l) > 0);
    set_Um = (U(:,l) < 0);

%-Useless constraints removal---------------------------------------------%    
    if (ismember(l,empty_set_plus))   
        t_minus_N = zeros(N,K);
        t_minus = zeros(1,K);
        
        for r = 1:K
            vect_e = -(Y_bar + U(:,id_diff)*T(id_diff,r))./U(:,l);
            t_plus(1,r) =  min(vect_e(set_Um)); 
        end      
    
        for n = 1:N
            for r = 1:K
                vect_e = -(Y_bar + U(:,id_diff)*(T(id_diff,r) + dT{n}(id_diff,r) + Z(id_diff,r)))./U(:,l);
                t_plus_N(n,r)  = min(vect_e(set_Um));
            end 
            var(n,:) = dT{n}(l,:);
        end
        A1 = [-1;-ones(N,1)];
        C1 = [t_plus; -ones(N,1)*Z(l,:) - var + t_plus_N];
%-------------------------------------------------------------------------%        
    elseif (ismember(l,empty_set_minus))
        t_plus_N = zeros(N,K);
        t_plus = zeros(1,K);
        
        for r = 1:K
            vect_e = -(Y_bar + U(:,id_diff)*T(id_diff,r))./U(:,l);
            t_minus(1,r) = max(vect_e(set_Up)); 
        end      
    
        for n = 1:N
            for r = 1:K
                vect_e = -(Y_bar + U(:,id_diff)*(T(id_diff,r) + dT{n}(id_diff,r) + Z(id_diff,r)))./U(:,l);
                t_minus_N(n,r) = max(vect_e(set_Up)); 
            end 
            var(n,:) = dT{n}(l,:);
        end
        A1 = [1;ones(N,1)];
        C1 = [-t_minus; ones(N,1)*Z(l,:) + var - t_minus_N];
%-------------------------------------------------------------------------%        
    else
        for r = 1:K
           vect_e = -(Y_bar + U(:,id_diff)*T(id_diff,r))./U(:,l);
           t_minus(1,r) = max(vect_e(set_Up));
           t_plus(1,r) =  min(vect_e(set_Um));
        end      
    
        for n = 1:N
            for r = 1:K
                vect_e = -(Y_bar + U(:,id_diff)*(T(id_diff,r) + dT{n}(id_diff,r) + Z(id_diff,r)))./U(:,l);
                t_minus_N(n,r) = max(vect_e(set_Up));
                t_plus_N(n,r)  = min(vect_e(set_Um));
            end 
            var(n,:) = dT{n}(l,:);
        end
        A1 = [1;-1;ones(N,1);-ones(N,1)];
        C1 = [-t_minus; t_plus; ones(N,1)*Z(l,:) + var - t_minus_N; -ones(N,1)*Z(l,:) - var + t_plus_N];
    end
%-------------------------------------------------------------------------%
    %-Lagrange's multiplier initialization 
    Lambda = zeros(length(A1),K);
    W = zeros(length(A1),K);
    W_prev = zeros(length(A1),K);

    %-Volume
    H = [T;ones(1,K)];
    H(l,:) = []; % row l removed
    H1 = H;      % H saved in an intermediary variable
    for k = 1:K
        H1(:,k) = [];                % column k removed
        d(k) = ((-1)^(k+l))*det(H1); % computation of the minor
        H1 = H;                      % return to the initial matrix
    end   
%-------------------------------------------------------------------------%    
    while ((norm_pri > eps_pri) || (norm_dual > eps_dual)) && (count < Niter_ADMM)
   
%         keyboard
        
        %-T(l,:) update
        T(l,:) = ((Y_proj(l,:) - S(l,:))*A' + muM(l)*(t_minus + t_plus + sum(t_minus_N + t_plus_N)) ...
                    - muM(l)*(N*Z(l,:) + sum(var)) + muM(l)*(A1')*(W-Lambda) ...
                   )/(A*(A') + muM(l)*eye(K) + beta*d*(d.')/(factorial(K-1)^2)); % V²(T) implemented to tackle the case det(H) < 0.
        
        %-Splitting variable update (projection)        
        G  = A1*T(l,:) + C1;
        W = G + Lambda;
        id3 = (W < 0);
        W(id3) = 0;
        
        %-Lagrange's mutliplier update
        Lambda = Lambda + G - W;  
        
        %-Error computation / count
        norm_pri = norm(G - W, 'fro');
        norm_dual = norm(-muM(l)*(A1')*(W - W_prev), 'fro');
        
        eps_pri = sqrt(2*(N+1)*K)*eps_abs + eps_rel*max([norm(A1*T(l,:),'fro'), norm(W,'fro'), norm(C1,'fro')]);
        eps_dual = sqrt(K)*eps_abs + eps_rel*norm(muM(l)*(A1')*Lambda,'fro');

        muM_prev = muM(l); % previous muM saved before potential update
        
        %-Regularisation parameter update
        if (norm_pri > mu*norm_dual)
            muM(l) = tau_incr*muM(l);
            %-Lagrange's multiplier update if muM updated
            Lambda = muM_prev*Lambda/muM(l);
        elseif (norm_dual > mu*norm_pri)
            muM(l) = muM(l)*tau_decr;
            %-Lagrange's multiplier update if muM updated
            Lambda = muM_prev*Lambda/muM(l);
        end  
        
        W_prev = W;
        count = count + 1; 
    end
end