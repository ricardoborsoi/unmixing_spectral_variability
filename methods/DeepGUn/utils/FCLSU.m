function out = FCLSU(HIM,M)
% Fully Constrained Linear Spectral Unmixing
% Perform a Linear least squares with nonnegativity constraints.
% --------------------------------------------------------------------
% Input:   HIM : hyperspectral image cube [nbands x nsamples]
%          M   : set of p endmembers [nbands x p].
% 
% Output:  out : fractions [p x nsamples] 
%
% 
% Copyright (2007) GRNPS group @ University of Extremadura, Spain. 


[~,ns] = size(HIM);
[l,p] = size(M);

Delta = 1/1000; % should be an small value

N = zeros(l+1,p);
N(1:l,1:p) = Delta*M;
N(l+1,:) = ones(1,p);
s = zeros(l+1,1);

out = zeros(ns,p);

%disp('Please wait...')
for i = 1:ns
    s(1:l) = Delta*HIM(:,i);
    s(l+1) = 1;
    %Abundances = M73lsqnonneg(N,s);
    Abundances = lsqnonneg(N,s);
    out(i,:) = Abundances;
end
% disp('End')