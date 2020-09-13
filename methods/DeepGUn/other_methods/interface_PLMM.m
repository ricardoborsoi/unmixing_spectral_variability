function [A_mtx_colord,dM_lexico2, M_tot] = interface_PLMM(data, A_init, M0, alpha, beta, gamma)
% PLMM wrapper
% 
% Input: data - data cube of observed image
%        A_init - initialization
%        
% Output:
% 
% 


[H,W,L] = size(data); % size of the hyperspectral datacube

Y = (reshape(permute(data,[2 1 3]),H*W,L))'; % pixel lexicographical ordering

K = size(M0,2); % desired endmember number
N = W*H; % number of pixels


% We use only the distance to M0 as the regularization
type = 'DISTANCE';



% % % % Regularization 
% % % type = 'DISTANCE';'NONE';'VOLUME';           % regularization type ('NONE','MUTUAL DISTANCE','VOLUME','DISTANCE')
% % % alpha = 0;1e-4;0;                 % abundance regularization parameter
% % % beta = 5.4e-4;             % endmember regularization parameter
% % % gamma = 1;                 % variability regularization parameter




% Convert row-wise lexocographic order back to cube
temp = permute(reshape(A_init',W,H,K),[2 1 3]);
% Convert cube to column-wise lexocographic order
A0 = reshape(temp, W*H, K)'; 



%--------------------------------------------------------------
% Unmixing parameters
%--------------------------------------------------------------
% Stopping criteria
% epsilon = 1e-2; % BCD iteration
% maxit_bcd = 20;
% % epsilon = 1e-3; 
% % maxit_bcd = 100;

% Stopping criteria
epsilon = 1e-3; % BCD iteration
maxit_bcd = 20;
Niter_ADMM = 100; 



eps_abs = 1e-2; % primal residual
eps_rel = 1e-4; % dual residual




eps_abs = 1e-2; % primal residual
eps_rel = 1e-4; % dual residual

% ADMM parameters
% Niter_ADMM = 30;          % maximum number of ADMM subiterations
% Niter_ADMM = 20;
muA = (1e-4)*L/(K+1);     % hyperparameters (AL constants)
muM = (1e-8)*N/(2*(N+1)*K); % Volume (1e-8)*N/(2*(N+1)*K) % Distance/None (1e-4)*N/((N+1)*K)
mudM = (1e-4)/(2*K);      % LN/(2LKN)
tau_incr = 1.1;
tau_decr = 1.1;
mu = 10;


%--------------------------------------------------------------
% Initialization
%--------------------------------------------------------------
% Endmember initialization (VCA [1])
% [M0, V, U, Y_bar, endm_proj, Y_proj] = find_endm(Y,K,'vca');
% % Abundance initialization (SUNSAL [2])
% A0 = sunsal(M0,Y,'POSITIVITY','yes','ADDONE','yes');
% % Abundance projection onto the unit simplex to strictly satisfy the constraints [3]
% for n = 1:N
%     A0(:,n) = ProjectOntoSimplex(A0(:,n),1);
% end

% Perturbation matrices initialization
dM0 = eps*ones(L,K*N);
dM0 = mat2cell(dM0,L,K*ones(N,1));    % [dM_1 | dM_2 | ... | dM_N]
% input = {type,beta,Y_proj,Y_bar,U,V}; % input = {type,beta} for 'MUTUAL DISTANCE' | input = {type,beta,M0} for 'DISTANCE'
                                      % input = {type,beta,Y_proj,Y_bar,U,V} for 'VOLUME'    
% input = {type,beta,M0};
                                      
                                      
%--------------------------------------------------------------
% BCD/ADMM unmixing (based on the PLMM)
%--------------------------------------------------------------
disp(['ADMM processing (M : ', type,', alpha = ',num2str(alpha),', beta = ',num2str(beta),', gamma = ',num2str(gamma),')...'])
tic
[f,A,M,dM] = bcd_admm(Y,A0,M0,dM0,W,H,gamma,eps_abs,eps_rel,epsilon,'HYPERPARAMETERS',{muA,muM,mudM},'PENALTY A',alpha,'PENALTY M',{type,beta,M0},'AL INCREMENT',{tau_incr,tau_decr,mu},'MAX ADMM STEPS',Niter_ADMM,'MAXITER BCD',maxit_bcd);
% [f,A,M,dM] = bcd_admm(Y,A0,M0,dM0,W,H,gamma,eps_abs,eps_rel,epsilon,'HYPERPARAMETERS',{muA,muM,mudM},'PENALTY A',alpha,'PENALTY M',{type,beta,Y_proj,Y_bar,U,V},'AL INCREMENT',{tau_incr,tau_decr,mu},'MAX ADMM STEPS',Niter_ADMM);
% [f,A,M,dM,RE,GMSE_A,GMSE_M,GMSE_dM,SA,mu_A,mu_M,mu_dM] = bcd_admm_th(Y,Ath,Mth,dMth,A0,M0,dM0,W,H,gamma,eps_abs,eps_rel,epsilon,'HYPERPARAMETERS',{muA,muM,mudM},'PENALTY A',alpha,'PENALTY M',{type,beta,Y_proj,Y_bar,U,V},'AL INCREMENT',{tau_incr,tau_decr,mu},'MAX ADMM STEPS',Niter_ADMM);
time =  toc;         

%--------------------------------------------------------------
% Error computation
%--------------------------------------------------------------
% [RE,A,var_map] = real_error(Y,A,M,dM,W,H);
% [RE_ADMM,GMSE_A_ADMM,GMSE_M_ADMM,MSE_dM_ADMM,GMSE_dM_ADMM,Error_A_ADMM,SA_ADMM,A,M,dM,var_map] = th_error(Y,Ath,Mth,dMth,A,M,dM,W,H);
disp('---------------------------------------------------------------------------');



% =========================================================================
% convert lexico ordering

% Convert row-wise lexocographic order back to cube
Aest_cube = permute(reshape(A',W,H,K),[2 1 3]);

% Convert cube to column-wise lexocographic order
A_mtx_colord = reshape(Aest_cube, W*H, K)'; 



dM_lexico2 = zeros(L,K,N);
temp = zeros(L,K,H,W);

n=0;
for i=1:H
    for j=1:W
        n = n+1;
        temp(:,:,i,j) = dM{n};
    end
end

n=0;
for j=1:W
    for i=1:H
        n = n+1;
        dM_lexico2(:,:,n) = temp(:,:,i,j);
    end
end


M_tot = zeros(L,K,N);
for n=1:N
    M_tot(:,:,n) = M + dM_lexico2(:,:,n);
end



