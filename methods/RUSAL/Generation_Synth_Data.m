%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%     Initialization     %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load spectres.txt;
spectres              = spectres-min(spectres(:));
[N_bandes N_spectres] = size(spectres);
N_spectres            = N_spectres-1;
wavelength            = spectres(:,1);
spectres(:,2:7)       = spectres(:,2:7)+abs(min(min(spectres(:,2:7))));
MPlus0                = spectres(1:4:end,[2:7]);
L                     = size(MPlus0,1);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%   Synthetic images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(R==3)
    interv = [4 5 6];
elseif(R==6)
    interv = randperm(6);interv = interv(1:R);
end
MPlus  = MPlus0(:,interv);
if(eta < size(MPlus,1)/6) 
        eta = size(MPlus,1)/6;
        disp(['corrected eta = ' num2str(eta)])
end

switch Synth_Im
    case 1
        %%% LMM %%%%%%%
        [alpha0 Y Y_bloc Y_bloc_SN sigma20 c_Ev0] = Generation_LMM_EV_SNR(R,N,Trunc,MPlus,SNR,eta,noise,L/2);
        gamma0   = zeros(D,N^2);
        s0       = 10^(-50)*ones(N);   % NxN
        RComp    =  zeros(L,N^2);      % LxN
        K        = 1;
        Label    = ones(1,N^2);
        Y        = reshape(Y_bloc',N,N,L);
    case 2
        %%% LMM + RCA-NL + GBM + PPM %%%%%%%
        [alpha0 gamma0 s b Y Y_bloc Y_bloc_SN sigma20 label_NL Lev_NL gamEqRCA]= Generation_LMM_RCA2_GBM_PPMM(R,N,Trunc,MPlus,SNR,noise,eta);
        s0       = reshape(Lev_NL,N,N); % NxN
        c_Ev0    = ones(N^2,1);
        RComp    = Y_bloc_SN - (MPlus*alpha0).*(c_Ev0*ones(1,L))' ; %  LxN
        K        = 4;
        Label    = label_NL;
        Y        = reshape(Y_bloc',N,N,L);
        
    case 3
        %%% LMM + RCA-EV + RCA-ME %%%%%%%
        [MPlus alpha0 epsi_ME epsi_EV Y Y_bloc Y_bloc_SN sigma20 label_EV gamEqRCA MPlusVar] ...
            = Generation_LMM_EV_ME(R,N,Trunc,MPlus,SNR,noise,eta);
        s0         = 10^(-50)*ones(N);   % NxN
        c_Ev0      = ones(N^2,1);
        gamma0     = zeros(D,N^2);
        RComp      = Y_bloc_SN - (MPlus*alpha0)  ; % LxN
        K          = 3;
        Label      = label_EV;
        %
end
row   = N;
col   = N;
L0    = L;
Bands = 1:L;