function Y = multrandn(m,Sigma,N)


%------------------------------------------------------------------
% This function allows ones to sample according to a multivariate 
% Gaussian distribution
% 
% % INPUT
%         m         : mean vector
%         Sigma     : covariance matrix
%         N         : dimension of the output
% 
%
% OUTPUT
%       Y  :  generated vectors
%
%------------------------------------------------------------------

n = length(m);

B = m;
A = chol(Sigma)';

X = randn(n,N);

Y = A*X+kron(B,ones(1,N));
