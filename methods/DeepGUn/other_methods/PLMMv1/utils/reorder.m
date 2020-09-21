function [M2,id] = reorder(M1,M2)
% Reorder the endmembers contained in M2.
%%
% Code : Pierre-Antoine Thouvenin, February 16th 2015.
%%
%-------------------------------------------------------------------------%
% Inputs:
% > M1     endmembers 1 (reference order)
% > M2     endmembers 2
% Outputs:
% < M2        reordered endmembers.
%-------------------------------------------------------------------------%
%%
K = size(M1,2);

%--------------------------------------------------------------
% Spectral angle computation, data reordering
%-------------------------------------------------------------- 
s = zeros(K,1);
SAM = zeros(K,1);
id = zeros(K,1);

for k = 1:K
    for l = 1:K
        s(l) = 180*acos( (M1(:,k).')*M2(:,l) /(norm(M1(:,k))*norm(M2(:,l))) )/pi; 
    end
    [SAM(k), id(k)] = min(s);
end

M2 = M2(:,id);