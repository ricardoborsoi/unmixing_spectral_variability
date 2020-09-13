function [A, optim_struct] = social_unmixing(data,sources,groups,A_init,lambda,rho,maxiter_ADMM,type,fraction,tol_a, verbose)

%   We try to minimize the following cost function:
%   J(A) = 1/2 * ||X - BA||_{F}^{2} + \lambda ||A||_{G,p,q}
%
%   with B a collection of endmember candidates, A the
%   abundances in each pixel and for each candidate.
%   G is the group structure used, and p and q can be:
%   (p,q) = (2,1) for the group penalty
%   (p,q) = (1,2) for the elitist penalty
%   (p,q) = (1,fraction) for the fractional penalty (0 < fraction <= 1)
%
%   The abundances are subject to the usual nonnegativity and sum to one
%   constraints.
%
% Inputs:
%
% -data = LxN data matrix with L the number of spectral bands, and N the
% number of pixels. 
% -sources = LxQ endmember matrix, with Q the total number of endmember 
% candidates.
% -groups = Qx1 vector indicating the group structure of the abundance
% matrix. Values have to range between 1 and P, the number of groups 
% (endmembers).
% -A_init: initial QxN abundance matrix: especially useful for the 
% fractional case, where the optimization problem is not convex.
% -lambda: regularization parameter for the sparsity inducing terms
% -type: string indicating which penalty to use: 'group' for group lasso,
% 'elitist' for elitist lasso, and 'fractional' for fractional lasso.
% -fraction: fraction to use for the fractional case (not used otherwise).
% 0 < fraction <= 1 (if fraction = 1, this is a regular lasso, otherwise
% the problem becomes nonconvex but allows to use the sum to one constraint
% while promoting both intra and inter group sparsity)
% -maxiter_ADMM: maximum number of iterations before the algorithm stops.
% -tol_a : tolerance on the relative variation of the norm of the abundance
% matrix under which the algorithm is considered converged
% -verbose: flag for display in console. Display if true, no display
% otherwise 
%
% Outputs: 
% -A: PxN abundance maps
% -optim_struct: structure containing the values of the cost function and
% its two terms (data fit and sparsity) for all the iterations.
%
% reference: 
%
%  Drumetz, L., Meyer, T. R., Chanussot, J., Bertozzi, A. L., & Jutten, C. 
%  (2019). Hyperspectral image unmixing with endmember bundles and group 
%  sparsity inducing mixed norms. IEEE Transactions on Image Processing,
%  28(7), 3435-3450.
%
% Author: Lucas Drumetz
% Latest Revision: 26-July-2019
% Revision: 1.3

%% INIT

P = length(groups);
G = length(unique(groups));
N = size(data,2);

objective = zeros(maxiter_ADMM,1);
norm_sparse = zeros(maxiter_ADMM,1);
data_fit = zeros(maxiter_ADMM,1);

U = A_init;
W = A_init;
E = zeros(size(U));

% define simplex projection

proj_simplex_array = @(y) max(bsxfun(@minus,y,max(bsxfun(@rdivide,cumsum(sort(y,1,'descend'),1)-1,(1:size(y,1))'),[],1)),0);

% precomputing

StS = (sources'*sources);
StSp2rhoIp = StS +2*rho*eye(P);


%% ALGORITHMS

if strcmp(type,'group')
    
    V = U;
    D = E;
    
    for i = 1:maxiter_ADMM
        
        U_old = U;
        
        U = StSp2rhoIp\ (sources'*data + rho *(V - D + W - E));
        
        V = prox_group_lasso(U+D,groups,lambda/rho);
        
        W = proj_simplex_array(U+E);
        
        D = D + U - V;
        E = E + U - W;
        
        group_norm = 0;
        
        for j = 1:G
            group_norm  = group_norm + sum((sqrt(sum(U(groups == j,:).^2))));
        end
        
        objective(i) = 1/2 * norm(data-sources*U,'fro')^2 + lambda* group_norm ;
        norm_sparse(i) = lambda * group_norm ;
        data_fit(i) = 1/2 * norm(data-sources*U,'fro')^2;
        
        rel_A = norm(U_old-U,'fro')/norm(U_old,'fro')^2;
        
        if verbose
            fprintf('iter %d of %d, objective = %f, fitting = %f, norm = %f, rel_A = %f\n',i,maxiter_ADMM,objective(i),data_fit(i),norm_sparse(i),rel_A)
        end
        
        if i>1 && rel_A < tol_a
            break
        end
        
    end
    
elseif strcmp(type,'elitist')
    
    V = U;
    D = E;
    
    for i = 1:maxiter_ADMM
        
        U_old = U;
        
        U = StSp2rhoIp\ (sources'*data + rho *(V + D + W + E));
        
        V = prox_elitist_group(U-D,groups,lambda/rho);    

        W = proj_simplex_array(U-E);
        
        D = D - U + V;
        E = E - U + W;
        
        partial_elitist_norm = 0;
        
        for j = 1:G
            partial_elitist_norm = partial_elitist_norm + (sum(abs(U(groups == j,:)))).^2;
        end
        
        elitist_norm  = sum(sqrt(partial_elitist_norm));
        
        objective(i) = 1/2 * norm(data-sources*U,'fro')^2 + lambda* elitist_norm;
        norm_sparse(i) = lambda * elitist_norm;
        data_fit(i) = 1/2 * norm(data-sources*U,'fro')^2;
        
        rel_A = norm(U_old-U,'fro')/norm(U_old,'fro')^2;
        
        if verbose
            fprintf('iter %d of %d, objective = %f, fitting = %f, norm = %f, rel_A = %f\n',i,maxiter_ADMM,objective(i),data_fit(i),norm_sparse(i),rel_A)
        end
        
        if i>1 && rel_A < tol_a
            break
        end
        
        
    end
    
                    
elseif strcmp(type,'fractional')
    
    M = zeros(G,P);
    
    for i = 1:G
        M(i,groups == i) = 1;
    end
    
    U_old = U;
    
    V = M*U;
    D = zeros(G,N);
    MtM = M'*M;
    
    StspMtMp2rhoIp = (StS+rho* MtM + rho*eye(P));
    
    for i = 1:maxiter_ADMM
        
        U = StspMtMp2rhoIp\ (sources'*data + rho*M'*(V + D) + rho*(W + E));
        V = approx_prox_fractional(M*U-D,lambda/rho,fraction);
        W = proj_simplex_array(U-E);
        
        D = D - M*U + V;
        E = E - U + W;
        
        fractional_norm = sum(sum(abs(U).^(fraction),1));
        
        objective(i) = 1/2 * norm(data-sources*U,'fro')^2 + lambda * fractional_norm;
        norm_sparse(i) = lambda * fractional_norm;
        data_fit(i) = 1/2 * norm(data-sources*U,'fro')^2;
        
        rel_A = norm(U_old-U,'fro')/norm(U_old,'fro')^2;
        
        if verbose
            fprintf('iter %d of %d, objective = %f, fitting = %f, norm = %f, rel_A = %f\n',i,maxiter_ADMM,objective(i),data_fit(i),norm_sparse(i),rel_A)
        end
        
        if i>1 && rel_A < tol_a
            break
        end
        
    end
    
    
end

A = U;
optim_struct = struct;
optim_struct.objective = objective;
optim_struct.norm = norm_sparse;
optim_struct.data_fit = data_fit;

end

