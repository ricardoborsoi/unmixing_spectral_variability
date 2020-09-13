function out = CLSU(HIM,M)
% Constrained Linear Spectral Unmixing
% Perform a Linear least squares with nonnegativity constraints.
% --------------------------------------------------------------------
% Input:   HIM : hyperspectral image cube [nrows x ncols x nchannels]
%          M   : set of p endmembers [nchannels x p].
% 
% Output:  out : fractions [nrows x ncols x p] 
%
% 
% Copyright (2007) GRNPS group @ University of Extremadura, Spain. 


[ns,nl,nb] = size(HIM);
[l,p] = size(M);

Delta = 1/1000; % should be an small value

N = Delta*M;

OutputImage = zeros(ns,nl,p);

%disp('Please wait...')
for i = 1:ns
    for j = 1:nl
        s = Delta*squeeze(HIM(i,j,:));
        %Abundances = M73lsqnonneg(N,s);
        Abundances = lsqnonneg(N,s);
        OutputImage(i,j,:) = Abundances;
    end
 end
% disp('End')

out = OutputImage;
