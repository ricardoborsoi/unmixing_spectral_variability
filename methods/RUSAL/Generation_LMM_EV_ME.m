function [MPlus alpha0 epsi_ME epsi_EV Y Y_bloc Y_bloc_SN sigma2 label_EV gamEqRCA MPlusVar] ...
    = Generation_LMM_EV_ME(R,N,Trunc,MPlus,SNR,noise,eta);


%% Generation_LMM_EV_ME(R,N,Trunc,MPlus,SNR,noise,eta) 
%  Paper :  A. Halimi, P. Honeine, J. Bioucas-Dias, "Hyperspectral Unmixing
%          in Presence of Endmember Variability, Nonlinearity or Mismodelling
%          Effects", IEEE Trans. Image Process., 2016.
%  Code  :  Generation with three classes: LMM, EV, ME
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
%%% MPlus     : Corrected endmembers (to avoid negative values) 
%%% alpha0    : Actual abundance values  (R x N^2)
%%% epsi_ME   : Variances or energies of the ME  (1 x 1)
%%% epsi_EV   : Variances or energies of the EV  (1 x 1)  
%%% Y         : pixels of size (N x N x L)
%%% Y_bloc    : pixels of size (L x N^2)
%%% Y_bloc_SN : noiseless pixels of size (L x N^2)
%%% sigma2    : noise variance (L x 1)
%%% label_EV  : label for the spatial classes (1 x N^2)
%%% Lev_NL    : Energy of the NL term (N^2 x 1)
%%% gamEqRCA  : ME expressed using RCA coefficients (D x N^2)
%%% MPlusVar  : Endmember of each pixel (L x R x N^2)
%% --------------------------------------------------------------------

 
%%% Class 1: LMM
%%% Class 2: EV       
%%% Class 3: ME 
if(N==100)
    load Labels_10000_3Classes.mat label
%     label = [ones(1,7500),2*ones(1,1250),3*ones(1,1250)];
elseif(N==50)
    load Labels_2500_3Classes.mat label
end
label_EV = label; clear label
   
%%% ME
L             = size(MPlus,1); 
Hdir          = toeplitz(exp(-[1:L].^2/100^2));
Kinv          = inv(Hdir);
Hdir          = chol(Kinv);
Kdir          = inv(Kinv); 
epsi_ME       = 0.002;  
gamEqRCA      = zeros(L,N^2);  

%%%%%%%%%%%%%%   Generation  %%%%%%%%%%%%%
Lev_NL = zeros(N^2,1);
Pdir0  = ones(R,1);% For random and non-uniform abundances
                   % in the simplex use: randi(20,R,1);
for i = 1:N
    for j = 1:N
        n           = i+(j-1)*N;
        alpha0(:,n) = [1 ;zeros(R-1,1)];
        while(sum(alpha0(:,n)>Trunc*ones(R,1)))
            % %             alpha0(:,n) = gamrnd(1,1,R,1);
            for r=1:R
                alpha0(r,n) = gamrnd(Pdir0(r,1),1);
            end
            alpha0(:,n) = alpha0(:,n)/sum(alpha0(:,n));
        end
    end
end
        
for i=1:3
        switch i
            case 1  % LMM
                ind              = find(label_EV==1); 
                nk               = length(ind);
                Y_bloc_SN(:,ind) =  MPlus *alpha0(:,ind); 
                gamEqRCA(:,ind)  = 0;clear ind
            case 2  % EV
                ind                   = find(label_EV==2); 
                nk                    = length(ind); 
                [MPlusVar DM epsi_EV] = generation_variabilite(MPlus,nk);
                Y_bloc_SN(:,ind)      = squeeze(sum(MPlusVar.*permute(repmat(alpha0(:,ind),[1 1 L]),[3 1 2]),2)) ;  % LxN
                gamEqRCA(:,ind)       = squeeze(sum(DM.*permute(repmat(alpha0(:,ind),[1 1 L]),[3 1 2]),2)) ;  % LxN
                clear ind
            case 3  % ME
          ind              = find(label_EV==3); 
          nk               = length(ind);
          gamEqRCA(:,ind)  = sqrt(epsi_ME)* chol(Kdir)'*randn(L,nk) ; 
          Y_bloc_SN(:,ind) =  MPlus *alpha0(:,ind) + gamEqRCA(:,ind);  
          indNeg           = find(sum( Y_bloc_SN(:,ind)<0));
          while(length(indNeg) >1)
              gamEqRCA(:,ind(indNeg))  =  gamEqRCA(:,ind(indNeg))/2;
              Y_bloc_SN(:,ind(indNeg)) =  MPlus *alpha0(:,ind(indNeg)) + gamEqRCA(:,ind(indNeg));
              indNeg                   = find(sum( Y_bloc_SN(:,ind)<0));
          end
          clear ind 
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

if(min(min(Y_bloc))<0) 
    MPlus     = MPlus - min(min(Y_bloc));
    Y_bloc    = Y_bloc - min(min(Y_bloc));
end  
Y         = reshape(Y_bloc',N,N,L);

 
