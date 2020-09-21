% -------------------------------------------------------------------------
% This code executes example 2 of the DeepGUn algorithm, referring to 
% the following publication:
% 
%   "Deep Generative Endmember Modeling: An Application to Unsupervised Spectral Unmixing"
%   Ricardo Augusto Borsoi, Tales Imbiriba, José Carlos Moreira Bermudez
%   IEEE Transactions on Computational Imaging, 2019
% 
% The DeepGUn algorithm performs spectral unmixing with spectral variability
% modeling the endmembers using deep generative models (variational autoencoders). 
% -------------------------------------------------------------------------

clear all
close all
clc


addpath(genpath('utils'))
addpath(genpath('other_methods'))
addpath(genpath('DeepGUn'))


clus = gcp('nocreate'); % If no pool, do not create new one.
if isempty(clus)
    c = parcluster('local');
    c.NumWorkers = 1; 5; 8; 7;
    parpool(c, c.NumWorkers);
end


% selects whether to use parfor or not
flag_useparfor = true; false;
if flag_useparfor == true
    flag_useparfor = inf;
else
    flag_useparfor = 0;
end



rng(5,'twister')





% -------------------------------------------------------------------------
% load data
load('data/data_ex2.mat')

[P,N]   = size(alphas);
[m,n,~] = size(alphas_cube);
L       = size(r,1);



%% Endmember Extraction

M0 = vca(r,'Endmembers',P);
id = zeros(P,1);
for k = 1:P
    for l = 1:P
        s(l) = 180*acos( (M(:,k).')*M0(:,l) /(norm(M(:,k))*norm(M0(:,l))) )/pi;
    end
    [~, id(k)] = min(s);
end
M0 = M0(:,id);




%% Fully Constrained Least Squares Unmixing (FCLSU)

disp('FCLSU...')
tic
A_FCLSU = FCLSU(r,M0)';
time_fcls = toc;

A_FCLSU_cube = matrixToHCube(A_FCLSU,m,n,1);


%% Full Extended Linear Mixing Model

disp('ELMM')

A_init = hCubeToMatrix(A_FCLSU_cube);
psis_init = ones(size(A_init));

% optional parameters
nnorm = '1,1'; % Use a Total Variation on the abundances
verbose = true; % display
maxiter_anls = 20;
maxiter_admm = 100;
epsilon_s = 10^(-3);
epsilon_a = 10^(-3);
epsilon_psi = 10^(-3);
epsilon_admm_abs = 10^(-2);
epsilon_admm_rel = 10^(-2);

    
% parameters example 2
    lambda_s = 1;
    lambda_a = 0.05;
    lambda_psi = 1e3;
    
tic
[A_ELMM, psis_ELMM, S_ELMM, ~] = ELMM_ADMM(r_cube, A_init, psis_init, M0,lambda_s,lambda_a,lambda_psi,nnorm,verbose,maxiter_anls,maxiter_admm,epsilon_s,epsilon_a,epsilon_psi,epsilon_admm_abs,epsilon_admm_rel,flag_useparfor);
time_elmm = toc;

A_ELMM = row2col_lexico_order(A_ELMM,m,n);
S_ELMM = row2col_lexico_order(S_ELMM,m,n);


%% Full Generalized Extended Linear Mixing Model

disp('GELMM')

maxiter_anls = 20;
maxiter_admm = 100;
epsilon_s = 10^(-3);
epsilon_a = 10^(-3);
epsilon_psi = 10^(-3);
epsilon_admm_abs = 10^(-2);
epsilon_admm_rel = 10^(-2);


% parameters example 2
    lambda_s   = 1;
    lambda_a   = 0.01;
    lambda_psi = 1e3;
    

psis_init = ones(L,P,N);

tic
[A_GELMM, psis_GELMM, S_GELMM, optim_struct] = GLMM_ADMM(r_cube, A_init, psis_init, M0,lambda_s,lambda_a,lambda_psi,nnorm,verbose,maxiter_anls,maxiter_admm,epsilon_s,epsilon_a,epsilon_psi,epsilon_admm_abs,epsilon_admm_rel,flag_useparfor);
time_glmm = toc;
A_GELMM = row2col_lexico_order(A_GELMM,m,n);
S_GELMM = row2col_lexico_order(S_GELMM,m,n);



%% PLMM

% parameters example 2
    t_alpha = 0.1;
    t_beta  = 1e-5;
    t_gamma = 10;

tic
[A_PLMM,dM_PLMM, M_PLMM] = interface_PLMM(r_cube, A_init,M0,t_alpha, t_beta,t_gamma);
time_plmm = toc;

A_PLMM = row2col_lexico_order(A_PLMM,m,n);
M_PLMM = row2col_lexico_order(M_PLMM,m,n);



%% Deep generative model
% ====================================================
disp('Test Deep Generative Model...')

% parameters of the algorithm
dimAut = 2;
lambda_zref = 0.1;
lambda_a = 0.01;

% parameters for the endmember bundle extraction
flag_Npx = true;
vec_Npx = [100 100 100];

tic
[A_deepGen,M_DeepGen]=DeepGUn(r_cube, A_init, M0, dimAut, lambda_zref, lambda_a, flag_Npx, vec_Npx, flag_useparfor);
time_deepGen = toc;





%% Show Results

A_FCLSU_im   = reshape(row2col_lexico_order(A_FCLSU,m,n)',m,n,P);
A_ELMM_im    = reshape(row2col_lexico_order(A_ELMM,m,n)',m,n,P);
A_PLMM_im    = reshape(row2col_lexico_order(A_PLMM,m,n)',m,n,P);
A_GELMM_im   = reshape(row2col_lexico_order(A_GELMM,m,n)',m,n,P);
A_deepGen_im = reshape(A_deepGen',m,n,P);


rmse_A_FCLS    = sqrt(norm(A_FCLSU_im(:)-alphas_cube(:))^2/(norm(alphas_cube(:))^2));
rmse_A_PLMM    = sqrt(norm(A_PLMM_im(:)-alphas_cube(:))^2/(norm(alphas_cube(:))^2));
rmse_A_ELMM    = sqrt(norm(A_ELMM_im(:)-alphas_cube(:))^2/(norm(alphas_cube(:))^2));
rmse_A_GELMM   = sqrt(norm(A_GELMM_im(:)-alphas_cube(:))^2/(norm(alphas_cube(:))^2));
rmse_A_DeepGen = sqrt(norm(A_deepGen_im(:)-alphas_cube(:))^2/(norm(alphas_cube(:))^2));

fprintf('\nRMSE for Abundances\n')
fprintf('FCLS.....: %.4f \n', rmse_A_FCLS)
fprintf('PLMM.....: %.4f \n', rmse_A_PLMM)
fprintf('ELMM.....: %.4f \n', rmse_A_ELMM)
fprintf('GLMM.....: %.4f \n', rmse_A_GELMM)
fprintf('DeepGen..: %.4f \n', rmse_A_DeepGen)



rmse_M_PLMM    = sqrt(norm(Mvs(:)-M_PLMM(:))^2/(norm(Mvs(:))^2));
rmse_M_ELMM    = sqrt(norm(Mvs(:)-S_ELMM(:))^2/(norm(Mvs(:))^2));
rmse_M_GELMM   = sqrt(norm(Mvs(:)-S_GELMM(:))^2/(norm(Mvs(:))^2));
rmse_M_Deepgen = sqrt(norm(Mvs(:)-M_DeepGen(:))^2/(norm(Mvs(:))^2));

fprintf('\nRMSE for Endmembers\n')
fprintf('PLMM.....: %.4f \n', rmse_M_PLMM)
fprintf('ELMM.....: %.4f \n', rmse_M_ELMM)
fprintf('GLMM.....: %.4f \n', rmse_M_GELMM)
fprintf('DeepGen..: %.4f \n', rmse_M_Deepgen)



acos_M_PLMM = 0;
acos_M_ELMM = 0;
acos_M_GELMM = 0;
acos_M_DeepGen = 0;

R_PLMM    = zeros(size(r));
R_ELMM    = zeros(size(r));
R_GELMM   = zeros(size(r));
R_DeepGen = zeros(size(r));

R_FCLS = M0*A_FCLSU;
A_deepGen2 = row2col_lexico_order(A_deepGen,m,n);
for i=1:N
    for j=1:P
        acos_M_PLMM    = acos_M_PLMM + acos(Mvs(:,j,i)'*M_PLMM(:,j,i)/(norm(Mvs(:,j,i))*norm(M_PLMM(:,j,i))));
        acos_M_ELMM    = acos_M_ELMM + acos(Mvs(:,j,i)'*S_ELMM(:,j,i)/(norm(Mvs(:,j,i))*norm(S_ELMM(:,j,i))));
        acos_M_GELMM   = acos_M_GELMM + acos(Mvs(:,j,i)'*S_GELMM(:,j,i)/(norm(Mvs(:,j,i))*norm(S_GELMM(:,j,i))));
        acos_M_DeepGen = acos_M_DeepGen + acos(Mvs(:,j,i)'*M_DeepGen(:,j,i)/(norm(Mvs(:,j,i))*norm(M_DeepGen(:,j,i))));
    end
    R_PLMM(:,i)    = squeeze(M_PLMM(:,:,i))*A_PLMM(:,i);
    R_ELMM(:,i)    = squeeze(S_ELMM(:,:,i))*A_ELMM(:,i);
    R_GELMM(:,i)   = squeeze(S_GELMM(:,:,i))*A_GELMM(:,i);
    R_DeepGen(:,i) = squeeze(M_DeepGen(:,:,i))*A_deepGen2(:,i);
end

acos_M_PLMM    = acos_M_PLMM/(N*P);
acos_M_ELMM    = acos_M_ELMM/(N*P);
acos_M_GELMM   = acos_M_GELMM/(N*P);
acos_M_DeepGen = acos_M_DeepGen/(N*P);

fprintf('\nSAM for Endmembers\n')
fprintf('PLMM.....: %.4f \n', acos_M_PLMM)
fprintf('ELMM.....: %.4f \n', acos_M_ELMM)
fprintf('GLMM.....: %.4f \n', acos_M_GELMM)
fprintf('DeepGen..: %.4f \n', acos_M_DeepGen)



rmse_r_FCLS    = sqrt(norm(r(:) - R_FCLS(:))^2/(norm(r(:))^2));
rmse_r_PLMM    = sqrt(norm(r(:) - R_PLMM(:))^2/(norm(r(:))^2));
rmse_r_ELMM    = sqrt(norm(r(:) - R_ELMM(:))^2/(norm(r(:))^2));
rmse_r_GELMM   = sqrt(norm(r(:) - R_GELMM(:))^2/(norm(r(:))^2));
rmse_r_DeepGen = sqrt(norm(r(:) - R_DeepGen(:))^2/(norm(r(:))^2));

fprintf('\nRMSE for R\n')
fprintf('FCLS.....: %.4f \n', rmse_r_FCLS)
fprintf('PLMM.....: %.4f \n', rmse_r_PLMM)
fprintf('ELMM.....: %.4f \n', rmse_r_ELMM)
fprintf('GLMM.....: %.4f \n', rmse_r_GELMM)
fprintf('DeepGen..: %.4f \n', rmse_r_DeepGen)


fprintf('\n\n')
fprintf('\nTime FCLS....: %.2f',time_fcls)
fprintf('\nTime PLMM....: %.2f',time_plmm)
fprintf('\nTime ELMM....: %.2f',time_elmm)
fprintf('\nTime GLMM....: %.2f',time_glmm)
fprintf('\nTime Proposed: %.2f',time_deepGen)
fprintf('\n\n')




%%



fh = figure;
[ha, pos] = tight_subplot(6, P, 0.01, 0.1, 0.1);
for i=1:P
    axes(ha(i)); j=1;
    imagesc(alphas_cube(:,:,i), [0 1]), set(gca,'ytick',[],'xtick',[])%, axis square
    axes(ha(i+j*P)); j=j+1;
    imagesc(A_FCLSU_im(:,:,i), [0 1]), set(gca,'ytick',[],'xtick',[])%, axis square
    axes(ha(i+j*P)); j=j+1;
    imagesc(A_PLMM_im(:,:,i), [0 1])
    set(gca,'ytick',[],'xtick',[]) %, axis square
    axes(ha(i+j*P)); j=j+1;
    imagesc(A_ELMM_im(:,:,i), [0 1]); 
    set(gca,'ytick',[],'xtick',[]) %, axis square
    axes(ha(i+j*P)); j=j+1;
    imagesc(A_GELMM_im(:,:,i), [0 1]) 
    set(gca,'ytick',[],'xtick',[])%, axis square
    axes(ha(i+j*P)); j=j+1;
    imagesc(A_deepGen_im(:,:,i), [0 1]); 
    set(gca,'ytick',[],'xtick',[])%, axis square
end

set(fh, 'Position', [0 0 550 700])
axes(ha(1));
title('EM \# 1','interpreter','latex')
axes(ha(2));
title('EM \# 2','interpreter','latex')
axes(ha(3));
title('EM \# 3','interpreter','latex')

j=0;
axes(ha(j*P+1)); j=j+1;
ylabel('True','interpreter','latex')
axes(ha(j*P+1)); j=j+1;
ylabel('FCLS','interpreter','latex')
axes(ha(j*P+1)); j=j+1;
ylabel('PLMM','interpreter','latex')
axes(ha(j*P+1)); j=j+1;
ylabel('ELMM','interpreter','latex')
axes(ha(j*P+1)); j=j+1;
ylabel('GLMM','interpreter','latex')
axes(ha(j*P+1)); j=j+1;
ylabel('DeepGUn','interpreter','latex')
colormap jet










