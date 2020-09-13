% % % % % % % % % % % % clc
clear all
global betag

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%        Constants       %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Synth_Im  = 3;          %  1: LMM  2: LMM+NL   3: ME+EV
noise     = 0;          %  For synthetic data
%  0: iid noise     1: Gauss. shape of the variance
eta       = 50;         %  Gauss. Width of the variance
%  eta > L/6,  L being the number of bands
SNR       = 25;         %  Noise level
R         = 3;          %  Number of endmembers (R in {3,6})
N         = 50;         %  Number of pixels for rows and columns
%  available: 50x50  or 100x100 pixels
Trunc     = 0.9999;     %  Maximum value for the abundances
D         = R*(R+1)/2;  %  Number of bilinear terms
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%  Synthetic/Real Data   %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Generation_Synth_Data

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%  Unmixing algorithms   %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%        Sunsal      %%%%%%%%%%%%%%%
disp('Sunsal')
algo = 1;
tic
lambda = 0;
alpha_sunsal  = sunsal(MPlus,Y_bloc,'lambda',lambda,'ADDONE','yes','POSITIVITY','yes', ...
    'AL_iters',200,'TOL', 1e-4, 'verbose','no');
time(algo) = toc;
[RMSE(algo) RE(algo) SAM(algo) SAMn(algo,:) RCest(:,:,algo) RE_RC(algo) SAM_RC(algo) SAMn_RC(algo,:)   ...
    RMSE_Class(algo,:) RE_Class(algo,:) SAM_Class(algo,:) RE_Class_RC(algo,:) SAM_Class_RC(algo,:) REn(algo,:)]  = ...
    exploitation_TIP(Y_bloc,MPlus,alpha_sunsal,alpha0,ones(row*col,1),c_Ev0,zeros(L,row*col),RComp,K,Label,1);
 
 
%%%%%%%%%%%%%         NUSAL      %%%%%%%%%%%%%%%
disp('NUSAL2')
Model = 'NUSAL2';   %NUSAL3

disp('------------------------------------')
disp('Tuning the regularization parameters')
disp('------------------------------------')
algo = 2;
vectTau    = [10^(-3)  10^(-2) 10^(-1)];%[10^(-6) 10^(-5) 10^(-4)];
AL_iters = 2000; tol  = 10^(-6);
k_NUSAL2=1;
for tau = vectTau
    for tau2 = 10*vectTau
        tic
        [zt,stop_NUSAL2,iter_NUSAL2] = NUSAL_RUSAL_v1(Y,MPlus,tau,tau2,AL_iters,tol,Model);
        time_NUSAL2(k_NUSAL2)   = toc;
        alpha_NUSAL2   = zt(1:R,:,end);
        resid_NUSAL2   = zt(R+1:end,:,end) ;% DxN
        [RMSE_NUSAL2(k_NUSAL2) RE_NUSAL2(k_NUSAL2) SAM_NUSAL2(k_NUSAL2) SAMn_NUSAL2(k_NUSAL2,:) RCest_NUSAL2(:,:,k_NUSAL2) RE_RC_NUSAL2(k_NUSAL2) SAM_RC_NUSAL2(k_NUSAL2) SAMn_RC_NUSAL2(k_NUSAL2,:)   ...
            RMSE_Class_NUSAL2(k_NUSAL2,:) RE_Class_NUSAL2(k_NUSAL2,:) SAM_Class_NUSAL2(k_NUSAL2,:) RE_Class_RC_NUSAL2(k_NUSAL2,:) SAM_Class_RC_NUSAL2(k_NUSAL2,:) REn_NUSAL2(k_NUSAL2,:)]  = ...
            exploitation_TIP(Y_bloc,MPlus,alpha_NUSAL2,alpha0,ones(row*col,1),c_Ev0,resid_NUSAL2,RComp,K,Label,2); %%%
        taut(k_NUSAL2,:) = [tau tau2];
        k_NUSAL2=k_NUSAL2+1;
    end
end

disp('------------------------------------')
disp('Runing the algorithm')
disp('------------------------------------')
[val pos_NUSAL2] = min(RMSE_NUSAL2);
[tau] = taut(pos_NUSAL2,1);[tau2 ] = taut(pos_NUSAL2,2);

AL_iters = 2000; tol  = 10^(-6);
tic
[zt,stop_NUSAL2,iter_NUSAL2] = NUSAL_RUSAL_v1(Y,MPlus,tau,tau2,AL_iters,tol,Model);
time(algo)  = toc;
alpha_NUSAL2   = zt(1:R,:,end);
resid_NUSAL2   = zt(R+1:end,:,end) ;% DxN
[RMSE(algo) RE(algo) SAM(algo) SAMn(algo,:) RCest(:,:,algo) RE_RC(algo) SAM_RC(algo) SAMn_RC(algo,:)   ...
    RMSE_Class(algo,:) RE_Class(algo,:) SAM_Class(algo,:) RE_Class_RC(algo,:) SAM_Class_RC(algo,:) REn(algo,:)]  = ...
    exploitation_TIP(Y_bloc,MPlus,alpha_NUSAL2,alpha0,ones(row*col,1),c_Ev0,resid_NUSAL2,RComp,K,Label,2); %%%



%%%%%%%%%%%%%       RUSAL      %%%%%%%%%%%%%%%
disp('RUSAL')
Model = 'RUSAL';

disp('------------------------------------')
disp('Tuning the regularization parameters')
disp('------------------------------------')
algo     = 3;
AL_iters = 2000; tol  = 10^(-6);
k_RUSAL=1;   vectTauRus = [10^(-3) 10^(-2) 10^(-1)];
for tau = vectTauRus
    for tau2 = 10*vectTauRus
        tic
        [zt,stop_RUSAL,iter_RUSAL] = NUSAL_RUSAL_v1(Y,MPlus,tau,tau2,AL_iters,tol,Model);
        time_RUSAL(k_RUSAL)   = toc;
        alpha_RUSAL    = zt(1:R,:,end);
        resid_RUSAL    = idct([zt(R+1:end,:,end); zeros(L-20,row*col)]);% LxN
        [RMSE_RUSAL(k_RUSAL) RE_RUSAL(k_RUSAL) SAM_RUSAL(k_RUSAL) SAMn_RUSAL(k_RUSAL,:) RCest_RUSAL(:,:,k_RUSAL) RE_RC_RUSAL(k_RUSAL) SAM_RC_RUSAL(k_RUSAL) SAMn_RC_RUSAL(k_RUSAL,:)   ...
            RMSE_Class_RUSAL(k_RUSAL,:) RE_Class_RUSAL(k_RUSAL,:) SAM_Class_RUSAL(k_RUSAL,:) RE_Class_RC_RUSAL(k_RUSAL,:) SAM_Class_RC_RUSAL(k_RUSAL,:) REn_RUSAL(k_RUSAL,:)]  = ...
            exploitation_TIP(Y_bloc,MPlus,alpha_RUSAL,alpha0,ones(row*col,1),c_Ev0,resid_RUSAL,RComp,K,Label,3); %%%
        taut(k_RUSAL,:) = [tau tau2];
        k_RUSAL=k_RUSAL+1;
    end
end
disp('------------------------------------')
disp('Runing the algorithm')
disp('------------------------------------')
[val pos_RUSAL] = min(RMSE_RUSAL);

[tau] = taut(pos_RUSAL,1);[tau2 ] = taut(pos_RUSAL,2);
AL_iters = 2000; tol  = 10^(-8);
tic
[zt,stop_RUSAL,iter_RUSAL] = NUSAL_RUSAL_v1(Y,MPlus,tau,tau2,AL_iters,tol,Model);

time(algo)   = toc;
alpha_RUSAL    = zt(1:R,:,end);
resid_RUSAL    = idct([zt(R+1:end,:,end); zeros(L-20,row*col)]);% LxN
[RMSE(algo) RE(algo) SAM(algo) SAMn(algo,:) RCest(:,:,algo) RE_RC(algo) SAM_RC(algo) SAMn_RC(algo,:)   ...
    RMSE_Class(algo,:) RE_Class(algo,:) SAM_Class(algo,:) RE_Class_RC(algo,:) SAM_Class_RC(algo,:) REn(algo,:)]  = ...
    exploitation_TIP(Y_bloc,MPlus,alpha_RUSAL,alpha0,ones(row*col,1),c_Ev0,resid_RUSAL,RComp,K,Label,3); %%%
clear taut 






%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%        Display results      %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CodeName = ['SUNSAL'; 'NUSAL2'; 'RUSAL ';];

disp(' ')
disp('**************************')
disp('Processing time (seconds)')
disp('**************************')
for i=1:3
    disp([num2str(CodeName(i,:)) ':  '  num2str(time(i))]);
end
disp('**************************')
disp(' ')
disp('****************')
disp('RMSE (x 10^(-2))')
disp('****************')
for i=1:3
    disp([num2str(CodeName(i,:)) ':  '  num2str(RMSE(i)*100)]);
end
disp('****************')
disp(' ')
disp('****************')
disp('RE (x 10^(-2))')
disp('****************')
for i=1:3
    disp([num2str(CodeName(i,:)) ':  '  num2str(RE(i)*100)]);
end
disp('****************')
disp(' ')
disp('****************')
disp('SAM (x 10^(-2))')
disp('****************')
for i=1:3
    disp([num2str(CodeName(i,:)) ':  '  num2str(SAM(i)*100)]);
end
disp('****************')



