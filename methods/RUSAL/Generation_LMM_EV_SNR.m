function [alpha0 Y Y_bloc Y_bloc_SN sigma2 c_EV] = Generation_LMM_EV_SNR(R,N,Trunc,MPlus,SNR,eta,noise,centre);
%% Generation_LMM_EV_SNR(R,N,Trunc,MPlus,SNR,eta,noise_type,centre);
%  Paper :  A. Halimi, P. Honeine, J. Bioucas-Dias, "Hyperspectral Unmixing
%          in Presence of Endmember Variability, Nonlinearity or Mismodelling
%          Effects", IEEE Trans. Image Process., 2016.
%  Code  :  Generation according to LMM 
%  Version (April, 2016) by Abderrahim Halimi (a.halimi@hw.ac.uk)
%  For any comments contact the author
%
%% --------------- Inputs --------------------------------------------
%%% R      : number of endmembers
%%% N      : number of pixels in rows (row=column=N)
%%% Trunc  : abundance truncature
%%% MPlus  : endmembers of size (L x R)
%%% SNR    : Noise level
%%% noise : Noise type
%%%          0: iid noise     1: Gauss. shape of the variance
%%% eta    : Gauss. Width of the variance
%%% centre : Centre of Gauss. shape of the variance

%%%
%% --------------- Outputs --------------------------------------------
%%% alpha0    : Actual abundance values  (R x N^2)  
%%% Y         : pixels of size (N x N x L)
%%% Y_bloc    : pixels of size (L x N^2)
%%% Y_bloc_SN : noiseless pixels of size (L x N^2)
%%% sigma2    : noise variance (L x 1)
%%% c_EV      : illumination coefficient (1 x N^2) 
%% --------------------------------------------------------------------


%%%%%%%%%% Illumination EV   
c_EV          = ones(N,1)*linspace(0.9,1.15,N);  
c_EV          = reshape(c_EV,1,N^2); 

L = size(MPlus,1);

%%%%%%%%%%%%%%   Generation  %%%%%%%%%%%%%
for i = 1:N
    for j = 1:N
        n           = i+(j-1)*N;
        alpha0(:,n) = [1 ;zeros(R-1,1)];
        while(sum(alpha0(:,n)>Trunc*ones(R,1)))
            alpha0(:,n) = gamrnd(1,1,R,1);
            alpha0(:,n) = alpha0(:,n)/sum(alpha0(:,n));
        end
    end
end


Y_bloc_SN = MPlus *(alpha0.*repmat(c_EV,R,1));
c_EV      =   repmat(MPlus,[1 1 N^2]) .* permute(repmat(c_EV-1,[R,1, L]),[3 1 2]);
variance  = sum(Y_bloc_SN(:).^2)/10^(SNR/10) /L/N^2 ;

if(noise==0)
    sigma2    = variance;
    Y_bloc    = Y_bloc_SN + sqrt(variance)*randn(L,N^2);
    Y         = reshape(Y_bloc',N,N,L);
    sigma2    = variance*ones(1,L);
elseif(noise==1) 
    quad_term = exp(-((1:L)-centre).^2/eta^2/2); 
    sigma2    = variance*L*quad_term/sum(quad_term);
    Y_bloc    = Y_bloc_SN + chol(diag(sigma2))'*randn(L,N^2);
    Y         = reshape(Y_bloc',N,N,L);
end

sigma2 = sigma2(:);

 
 