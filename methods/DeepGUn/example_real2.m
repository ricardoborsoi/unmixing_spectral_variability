% -------------------------------------------------------------------------
% This code executes the Samson example of the DeepGUn algorithm,  
% referring to the following publication:
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
    c.NumWorkers = 1; 8; 7;
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
% load data and initial reference endmembers
load('data/real_data_samson.mat')

[m,n,L] = size(r_cube);
r = reshape(r_cube,m*n,L)';
P = size(M0,2);
N = m*n;

% identify materials
materials{1} = 'Water';
materials{2} = 'Vegetation';
materials{3} = 'Soil';



%% Fully Constrained Least Squares Unmixing (FCLSU)

disp('FCLSU...')
tic
A_FCLSU = FCLSU(r,M0)';
time_fcls = toc;

A_FCLSU_cube = matrixToHCube(A_FCLSU,m,n,1);


%% Full Extended Linear Mixing Model

disp('ELMM')

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

    
% parameters example Samson
    lambda_s = 0.4;
    lambda_a = 0.005;
    lambda_psi = 0.001;
    
tic
[A_ELMM, psis_ELMM, S_ELMM, ~] = ELMM_ADMM(r_cube, A_init, psis_init, M0,lambda_s,lambda_a,lambda_psi,nnorm,verbose,maxiter_anls,maxiter_admm,epsilon_s,epsilon_a,epsilon_psi,epsilon_admm_abs,epsilon_admm_rel,flag_useparfor);
time_elmm = toc;



%% Full Generalized Extended Linear Mixing Model

disp('GELMM')

maxiter_anls = 20;
maxiter_admm = 100;
epsilon_s = 10^(-3);
epsilon_a = 10^(-3);
epsilon_psi = 10^(-3);
epsilon_admm_abs = 10^(-2);
epsilon_admm_rel = 10^(-2);


% parameters example Samson
    lambda_s   = 1;
    lambda_a   = 0.01;
    lambda_psi = 1e-6;

psis_init = ones(L,P,N);

tic
[A_GELMM, psis_GELMM, S_GELMM, optim_struct] = GLMM_ADMM(r_cube, A_init, psis_init, M0,lambda_s,lambda_a,lambda_psi,nnorm,verbose,maxiter_anls,maxiter_admm,epsilon_s,epsilon_a,epsilon_psi,epsilon_admm_abs,epsilon_admm_rel,flag_useparfor);
time_glmm = toc;




%% PLMM

% parameters example Samson
    t_alpha = 0.0014;
    t_beta  = 500;
    t_gamma = 1;

tic
[A_PLMM,dM_PLMM, M_PLMM] = interface_PLMM(r_cube, A_init,M0,t_alpha, t_beta,t_gamma);
time_plmm = toc;





%% Deep generative model
% ====================================================
disp('Test Deep Generative Model...')

% parameters of the algorithm
dimAut = 2;
lambda_zref = 0.1;
lambda_a = 0.005;

% parameters for the endmember bundle extraction
flag_Npx = true;
vec_Npx = [100 100 100];

tic
[A_deepGen,M_DeepGen]=DeepGUn(r_cube, A_init, M0, dimAut, lambda_zref, lambda_a, flag_Npx, vec_Npx, flag_useparfor);
time_deepGen = toc;





%% Show Results

A_FCLSU_im   = reshape(A_FCLSU',m,n,P);
A_ELMM_im    = reshape(A_ELMM',m,n,P);
A_PLMM_im    = reshape(A_PLMM',m,n,P);
A_GELMM_im   = reshape(A_GELMM',m,n,P);
A_deepGen_im = reshape(A_deepGen',m,n,P);


R_PLMM    = zeros(size(r));
R_ELMM    = zeros(size(r));
R_GELMM   = zeros(size(r));
R_DeepGen = zeros(size(r));

R_FCLS = M0*A_FCLSU;
for i=1:N
    R_PLMM(:,i)    = squeeze(M_PLMM(:,:,i))*A_PLMM(:,i);
    R_ELMM(:,i)    = squeeze(S_ELMM(:,:,i))*A_ELMM(:,i);
    R_GELMM(:,i)   = squeeze(S_GELMM(:,:,i))*A_GELMM(:,i);
    R_DeepGen(:,i) = squeeze(M_DeepGen(:,:,i))*A_deepGen(:,i);
end

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
[ha, pos] = tight_subplot(5, P, 0.01, 0.1, 0.1);
for i=1:P
    axes(ha(i)); j = 1;
    imagesc(A_FCLSU_im(:,:,i),[0 1]), set(gca,'ytick',[],'xtick',[])
    axes(ha(i+j*P)); j = j+1;
    imagesc(A_PLMM_im(:,:,i),[0 1]), set(gca,'ytick',[],'xtick',[])
    axes(ha(i+j*P)); j = j+1;
    imagesc(A_ELMM_im(:,:,i),[0 1])
    set(gca,'ytick',[],'xtick',[])
    axes(ha(i+j*P)); j = j+1;
    imagesc(A_GELMM_im(:,:,i),[0 1])
    set(gca,'ytick',[],'xtick',[])
    axes(ha(i+j*P)); j = j+1;
    imagesc(A_deepGen_im(:,:,i),[0 1])
    set(gca,'ytick',[],'xtick',[])
end


set(fh, 'Position', [0 0 550 700])
for i=1:P
    axes(ha(i));
    title(materials{i},'interpreter','latex')
end
axes(ha(1)); j = 1;
ylabel('FCLS','interpreter','latex')
axes(ha(1+j*P)); j = j+1;
ylabel('PLMM','interpreter','latex')
axes(ha(1+j*P)); j = j+1;
ylabel('ELMM','interpreter','latex')
axes(ha(1+j*P)); j = j+1;
ylabel('GLMM','interpreter','latex')
axes(ha(1+j*P)); j = j+1;
ylabel('DeepGUn','interpreter','latex')


colormap(jet)










