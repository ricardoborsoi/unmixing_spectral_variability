function [A,psi,cost] = ADMM_collaborative_unmixing(X,A_init,S,lambda,rho,maxiter_admm,type,tol_a,verbose)


%   We try to minimize the following cost function:
%   J(A) = 1/2 * ||X - BA||_{F}^{2} + \lambda ||A||_{2,1}^{2})
%
%   with B a collection of endmember candidates, A the
%   abundances in each pixel and for each candidate. This approach is
%   referred to as "collaborative sparsity", e.g. in 
%
%   Iordache, M. D., Bioucas-Dias, J. M., & Plaza, A. (2013). 
%   Collaborative sparse regression for hyperspectral unmixing.
%   IEEE Transactions on Geoscience and Remote Sensing, 52(1), 341-354.
%
%   The abundances are subject to the usual nonnegativity and sum to one
%   constraints. 
%
% Inputs:
%
% -data = LxN data matrix with L the number of spectral bands, and N the
% number of pixels. 
% -A_init: initial QxN abundance matrix: especially useful for the 
% -lambda: regularization parameter for the sparsity inducing terms
% -maxiter_admm: maximum number of iterations before the algorithm stops.
% -type: string indicating whether or not the abundance sum-to-one
% constraint (asc) should be enforced. Input 'asc' to enforce it.
% - tol_a : tolerance on the relative variation of the norm of the abundance
% matrix under which the algorithm is considered converged
% -verbose: flag for display in console. Display if true, no display
% otherwise 
%
% Outputs: 
% -A: PxN abundance maps
% -optim_struct: 
%
% reference: 
%
%  Drumetz, L., Meyer, T. R., Chanussot, J., Bertozzi, A. L., & Jutten, C. 
%  (2019). Hyperspectral image unmixing with endmember bundles and group 
%  sparsity inducing mixed norms. IEEE Transactions on Image Processing,
%  28(7), 3435-3450.
%
% Author: Lucas Drumetz
% Latest Revision: 22-July-2019
% Revision: 1.2

% Init and precomputing

[P,N] = size(A_init);

phi = A_init;
V = A_init;
U = A_init;

tol = tol_a;

C = zeros(P,N);
D = zeros(P,N);

% define simplex projection

proj_simplex_array = @(y) max(bsxfun(@minus,y,max(bsxfun(@rdivide,cumsum(sort(y,1,'descend'),1)-1,(1:size(y,1))'),[],1)),0);

SSplusI = (S'*S+2*rho*speye(P));
Sx = S'*X;

cost = zeros(maxiter_admm,1);

   for i = 1:maxiter_admm
   
   phi_old = phi;    
   
   phi = SSplusI\(Sx + rho *(U+V+C+D));
   U = vector_soft_col((phi-C)', lambda/rho)'; 
   if strcmp(type,'asc')
     V = proj_simplex_array(phi-D);  
   else
     V = max((phi-D),0) ;
   end
  
   C = C + U - phi;
   D = D + V - phi;
   
   norm_p = zeros(P,1);
   
   for p = 1:P
      norm_p(p) = norm(phi(p,:)); 
   end
   
   cost(i) = 1/2*norm(X-S*phi,'fro')^2 + lambda* sum(norm_p);
   
   rel_phi = norm(phi_old-phi,'fro')/norm(phi_old,'fro');
   
   if verbose
       fprintf('iteration %d of %d, cost = %d, rel_phi = %d\n ',i,maxiter_admm, cost(i), rel_phi)
   end
   
   if i > 2 & rel_phi < tol
      break 
   end
   
   end

%% get scaling factors and abundances

if strcmp(type,'SCLSU')
    psi = sum(phi);
    A = zeros(size(phi));
    
    for i = 1:P
        A(i,:) = phi(i,:)./psi;
    end
else
    A = phi;
end
  
    
    
cost(cost == 0) = [];

end