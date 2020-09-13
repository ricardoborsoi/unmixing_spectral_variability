function [A_global,sources_global] = bundle2global(A_bundle,bundle,groups)

% This function sums abundance maps corresponding to the same group and 
% computes pixelwise endmembers from the variant abundances
%
% Inputs:
% -A_bundle is the Q*N abundance matrix of all considered signatures
% - bundle is the LxQ matrix containing these signatures
% - groups is the Qx1 vector indicating the group structure of the abundance
% matrix. Values have to range between 1 and P, the number of groups 
% (endmembers).
% 
% Author: Lucas Drumetz
% Latest Revision: 22-July-2019
% Revision: 1.2

N = size(A_bundle,2);
L = size(bundle,1);
nbg = max(groups);

A_global = zeros(nbg,N);
sources_global = zeros(L,nbg,N);


threshold = 10^(-2);

% test

A_bundle_new = A_bundle;
A_bundle_new(abs(A_bundle) < threshold) = 0;

for p = 1:nbg
    A_global(p,:) = sum(A_bundle_new(groups == p,:),1);

    
    for i = 1:N
        if A_global(p,i) ~=0
            sources_global(:,p,i) = sum(repmat(A_bundle_new(groups == p,i),1,L)'.*bundle(:,groups ==p),2)/A_global(p,i);
        else
            sources_global(:,p,i) = mean(bundle(:,groups == p),2);
        end
    end
end



end

