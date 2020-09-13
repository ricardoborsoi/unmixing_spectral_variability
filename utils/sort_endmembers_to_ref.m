function [M0, idxPerm] = sort_endmembers_to_ref(Mref,Mest)
% input:
% Mref - Desired/reference endmember matrix
% Mest - estimated endmemebrs
%
% output:
% M0 - reordered endmember matrix
% idxPerm - index used to permute 'Mest' into 'M0'
% Ricardo Borsoi, 2019

if size(Mref,2) ~= size(Mest,2) || size(Mref,1) ~= size(Mest,1)
    error('matrix sizes do not agree!')
end

P = size(Mref,2);


M0 = Mest;
Mth = Mref;

% Sort M0 with respect to real/desired EM signatures to ease the comparison of estimated abundance maps
id = zeros(P,1);
for k = 1:P
    for l = 1:P
        s(l) = 180*acos( (Mth(:,k).')*M0(:,l) /(norm(Mth(:,k))*norm(M0(:,l))) )/pi; 
    end
    [~, id(k)] = min(s);
end
M0 = M0(:,id);

idxPerm = id;
