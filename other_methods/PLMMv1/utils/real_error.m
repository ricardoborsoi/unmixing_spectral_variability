function [RE,A,var_map] = real_error(Y,A,M,dM,W,H)
% Error computation for rela data.
%%
% Code : Pierre-Antoine Thouvenin, February 16th 2015.
%%
%-------------------------------------------------------------------------%
% Inputs:
% >   Y     lexicographically ordered hyperspectral data (L|N);
% >   A     abundance matrix (K|N);
% >   M     endmember matrix (L|K);
% >  dM     variability matrices (L|NK)(cell 1|N);
% >   W     width of the image; (N = W*H)
% >   H     heigth of th image.
%
% Outputs:
% < RE       reconstruction mean square error; || Y - MA - Delta || / sqrt(LN)
% <  A       re-shaped abundance matrix (K|N);
% < var_map  variability energy map associated to each endmember (H|W|K).
%-------------------------------------------------------------------------%
N = W*H;
[L,K] = size(M);

%--------------------------------------------------------------
% Abundance data (lexicographical order)
%--------------------------------------------------------------
if (size(A,3) > 1)
    A = (reshape(permute(A,[2 1 3]),H*W,K))';  %(K|N) : abundance map
end

%--------------------------------------------------------------
% RE computation
%--------------------------------------------------------------
if isempty(dM)
    RE = (norm(Y-M*A,'fro')^2)/(L*N);
else
    dMA = zeros(L,N);
    for n = 1:N
        dMA(:,n) = dM{n}*A(:,n);
    end
    RE = (norm(Y-M*A-dMA,'fro')^2)/(L*N);
end

% Abundance reshape
if (size(A,3) <= 1)
    A = permute(reshape(A',W,H,K),[2 1 3]); %(H|W|K) : abundance cube    
end

%--------------------------------------------------------------
% Variability map
%--------------------------------------------------------------
if ~isempty(dM)
    dM = reshape(dM,W,H)';    %(H|W) : perturbation map
    var_map = zeros(H,W,K);
    for p = 1:K 
        for k = 1:H
            for l = 1:W
                var_map(k,l,p) = (norm(dM{k,l}(:,p),2))/sqrt(L);
            end
        end
    end
else
    var_map = [];
end