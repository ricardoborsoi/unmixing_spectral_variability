function [RE,GMSE_A,GMSE_M,MSE_dM,GMSE_dM,Error_A,SA,A,M,dM,var_map] = th_error(Y,Ath,Mth,dMth,A,M,dM,W,H)
% Error computation.
%%
% Code : Pierre-Antoine Thouvenin, February 16th 2015.
%%
%-------------------------------------------------------------------------%
% Inputs:
% >   Y     lexicographically ordered hyperspectral data (L|N);
% >  Ath    theoretical abundance matrix (K|N);
% >  Mth    theoretical endmember matrix (L|K);
% > dMth    theoretical perturbation matrices (L|NK)(cell 1|N);
% >   A     abundance matrix (K|N);
% >   M     endmember matrix (L|K);
% >  dM     variability matrices (L|NK)(cell 1|N);
% >   W     width of the image; (N = W*H)
% >   H     heigth of th image.
%
% Outputs:
% < RE        reconstruction mean square error; || Y - MA - Delta || / sqrt(LN)
% < GMSE_A    reconstruction mean square error; || Ath - A || / sqrt(KN)
% < GMSE_M    reconstruction mean square error; || Mth - M || / sqrt(LK)
% < MSE_dM    reconstruction mean square error (H|W); || dM_th - dM || / sqrt(LK)
% < Error_A   normalized abundance error;
% <   SA      spectral angle;
% <   A       re-ordered abundance matrix (K|N);
% <   M       re-ordered endmember matrix (L|K);
% <  dM       re-ordered variability matrix (L|K)(cell 1|N);
% < var_map   perturbation energy map associated to each endmember (H|W|K).
%-------------------------------------------------------------------------%
%%
N = W*H;
[L,K] = size(Mth);

%--------------------------------------------------------------
% Abundance data cube (lexicographical order)
%--------------------------------------------------------------

if (size(A,3) <= 1)
    A = permute(reshape(A',W,H,K),[2 1 3]);      %(H|W|K)
end

if (size(Ath,3) <= 1)
    Ath1 = permute(reshape(Ath',W,H,K),[2 1 3]); %(H|W|K)
end


%--------------------------------------------------------------
% Spectral angle computation, data reordering
%-------------------------------------------------------------- 
s = zeros(K,1);
SA = zeros(K,1);
id = zeros(K,1);

for k = 1:K
    for l = 1:K
        s(l) = 180*acos( (Mth(:,k).')*M(:,l) /(norm(Mth(:,k))*norm(M(:,l))) )/pi; 
    end
    [SA(k), id(k)] = min(s);
end

A = A(:,:,id);
M = M(:,id);
SA = SA(id);

Error_A = zeros(K,1); 
for k = 1:K
    Error_A(k) = norm(Ath1(:,:,k) - A(:,:,k),'fro')/norm(Ath1(:,:,k),'fro');
end

A = (reshape(permute(A,[2 1 3]),H*W,K))'; %(K|N)

%--------------------------------------------------------------
% Error computation
%--------------------------------------------------------------
GMSE_A = (norm(Ath - A,'fro')^2)/(K*N);
GMSE_M = (norm(Mth - M,'fro')^2)/(L*K);

dMA = zeros(L,N);
if isempty(dMth)||isempty(dM)
    MSE_dM = [];
    GMSE_dM = [];
else
    MSE_dM = zeros(1,N);
    for n = 1:N
        dM{n} = dM{n}(:,id); % perturbation reordering
        dMA(:,n) = dM{n}*A(:,n);
        MSE_dM(n) = (norm(dMth{n} - dM{n},'fro')^2)/(L*K);
    end
    GMSE_dM = mean(MSE_dM);
    MSE_dM = reshape(MSE_dM,W,H)';
end
RE = (norm(Y-M*A-dMA,'fro')^2)/(L*N);

A = permute(reshape(A',W,H,K),[2 1 3]);

%--------------------------------------------------------------
% Variability map
%--------------------------------------------------------------
if ~isempty(dM)
    dM1 = reshape(dM,W,H)';    %(H|W) : variability map
    var_map = zeros(H,W,K);
    for p = 1:K 
        for k = 1:H
            for l = 1:W
                var_map(k,l,p) = (norm(dM1{k,l}(:,p),2))/sqrt(L);
            end
        end
    end
end