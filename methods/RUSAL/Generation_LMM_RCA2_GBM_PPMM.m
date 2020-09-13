function [alpha0 gamma0 s b Y Y_bloc Y_bloc_SN sigma2 label_NL Lev_NL gamEqRCA] ...
    = Generation_LMM_RCA2_GBM_PPMM(R,N,Trunc,MPlus,SNR,noise,eta);

%% Generation_LMM_RCA2_GBM_PPMM(R,N,Trunc,MPlus,SNR,noise,eta);
%  Paper :  A. Halimi, P. Honeine, J. Bioucas-Dias, "Hyperspectral Unmixing
%          in Presence of Endmember Variability, Nonlinearity or Mismodelling
%          Effects", IEEE Trans. Image Process., 2016.
%  Code  :  Generation with four classes: LMM, RCA, GBM and PPNMM
%  Version (April, 2016) by Abderrahim Halimi (a.halimi@hw.ac.uk)
%  For any comments contact the author
%
%% --------------- Inputs --------------------------------------------
%%% R      : number of endmembers
%%% N      : number of pixels in rows (row=column=N)
%%% Trunc  : abundance truncature
%%% MPlus  : endmembers of size (L x R)
%%% SNR    : Noise level
%%% noise  : Noise type
%%%          0: iid noise     1: Gauss. shape of the variance
%%% eta    : Gauss. Width of the variance

%%%
%% --------------- Outputs --------------------------------------------
%%% alpha0    : Actual abundance values  (R x N^2)
%%% gamma0    : Actual nonlinear coeffcients for the RCA class (D x N^2)
%%% s         : Variances or energies of the gamma0  (N^2 x 1)
%%% b         : Nonlinearity coefficient for PPNMM (N^2 x 1)
%%% Y         : pixels of size (N x N x L)
%%% Y_bloc    : pixels of size (L x N^2)
%%% Y_bloc_SN : noiseless pixels of size (L x N^2)
%%% sigma2    : noise variance (L x 1)
%%% label_NL  : label for the spatial classes (1 x N^2)
%%% Lev_NL    : Energy of the NL term (N^2 x 1)
%%% gamEqRCA  : NL expressed using RCA coefficients (D x N^2)
%% --------------------------------------------------------------------


%%%% Modification: abondance K=4 -> K=1
%%% Class 1: LMM
%%% Class 2: RCA-2     s=0.005
%%% Class 3: GBM-Fan
%%% Class 4: PPMM      b=0.5 or rand(0,1)
if(N==100)
    load Labels_10000_4Classes.mat label
elseif(N==50)
    load Labels_2500_4Classes.mat label
end
label_NL = label;

%%%%%%%%%% Non-linearity
b             = zeros(1,N^2); % PPMM
p_GBM         = zeros(R,R);   % GBM

%%% RCA
M_NL2    = [];M_NL3    = [];
for i=1:R
    for j=i+1:R
        M_NL2   = [M_NL2  MPlus(:,i).*MPlus(:,j)];
    end
end
M_NL2           = [sqrt(2)*M_NL2 MPlus.^2]; % LxD
M_NL            = M_NL2;
D               = size(M_NL,2);
L               = size(MPlus,1);
s               = 0.005;
pos             = find(label_NL==2);
gamma0          = zeros(D,N^2);
gamma0(:,pos)   = abs(sqrt(s)* randn(D,length(pos)));
gamEqRCA        = zeros(D,N^2);
gamEqRCA(:,pos) = gamma0(:,pos);


%%%%%%%%%%%%%%   Generation  %%%%%%%%%%%%%
Lev_NL = zeros(N^2,1);
Pdir0  = ones(R,1);% For random and non-uniform abundances
                   % in the simplex use: randi(20,R,1);
for i = 1:N
    for j = 1:N
        n           = i+(j-1)*N;
        alpha0(:,n) = [1 ;zeros(R-1,1)];
        while(sum(alpha0(:,n)>Trunc*ones(R,1))) 
            for r=1:R
                alpha0(r,n) = gamrnd(Pdir0(r,1),1);
            end
            alpha0(:,n) = alpha0(:,n)/sum(alpha0(:,n));
        end 
        
        switch label_NL(n)
            case 1  % LMM
                Y_bloc_SN(:,n) =  MPlus *alpha0(:,n);
                
            case 2  % RCA
                Y_bloc_SN(:,n) =  MPlus *alpha0(:,n)+  M_NL * gamma0(:,n);
                Lev_NL(n)      = sum((M_NL * gamma0(:,n)).^2);
                
            case 3  % GBM
                Y_bloc_SN(:,n) =  MPlus *alpha0(:,n) ;f_NL=0;
                pii=1;
                for ii=1:R-1
                    for jj= ii+1:R
                        p_GBM(ii,jj)    = 0.8+0.2*rand;
                        gamEqRCA(pii,n) = p_GBM(ii,jj)*alpha0(ii,n) *alpha0(jj,n);
                        Y_bloc_SN(:,n)  = Y_bloc_SN(:,n) +  p_GBM(ii,jj)*(MPlus(:,ii).*MPlus(:,jj))*alpha0(ii,n) *alpha0(jj,n) ;
                        f_NL            = f_NL+ (( p_GBM(ii,jj)*(MPlus(:,ii).*MPlus(:,jj))*alpha0(ii,n) *alpha0(jj,n)));
                        pii             = pii+1;
                    end
                end
                Lev_NL(n)  = sum(f_NL.^2); 
                
            case 4  % PPMM
                b(n)           = 0.5; %rand;
                Y_bloc_SN(:,n) = MPlus *alpha0(:,n);
                Y_bloc_SN(:,n) = Y_bloc_SN(:,n) +  b(n)*Y_bloc_SN(:,n).^2;
                pii=1;
                for ii=1:R-1
                    for jj= ii+1:R
                        gamEqRCA(pii,n) = b(n)*2*alpha0(ii,n) *alpha0(jj,n);
                        pii             = pii+1;
                    end
                end
                gamEqRCA(pii:pii+R-1,n) = b(n)*alpha0(:,n).^2;
                Lev_NL(n)               = sum((b(n)*Y_bloc_SN(:,n).^2).^2);
        end
    end
end


variance  = sum(Y_bloc_SN(:).^2)/10^(SNR/10) /L/N^2 ;
if noise==0
    quad_term = ones(1,L);
elseif(noise==1)
    quad_term = exp(-((1:L)-L/2).^2/eta^2/2);
end
sigma2    = variance*L*quad_term/sum(quad_term);
Y_bloc    = Y_bloc_SN + repmat(sqrt(sigma2(:)),1,N^2).*randn(L,N^2);
Y         = reshape(Y_bloc',N,N,L);
