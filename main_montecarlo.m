% =========================================================================
% This code compares different specral unmixing algorithms that consider
% spectral variability using realistic synthetically generated endmember 
% data.
%  *** This code follows the same structure of main.m, but uses a Monte 
%      Carlo simulation to evaluate the algorithms more reliably ***
% 
% The code is provided as part of the following publication:
%     Spectral Variability in Hyperspectral Data Unmixing: A 
%     Comprehensive Review
% 
% =========================================================================

clear all
close all
clc

addpath(genpath('synthetic_data_generation'))
addpath(genpath('library_extraction'))
addpath(genpath('utils'))
addpath('methods')


% number of Monte Carlo realizations
num_runs = 10;


% initialize the vectors to store the quantitative metrics
RMSEA(1) = nan;
RMSEY(1) = nan;
RMSEM(1) = nan;
SAMM(1)  = nan;
TIMES(1) = nan;


for mc_run = 1:num_runs
rng(mc_run,'twister');

% select the signal to noise ratio:
SNR = 30;

% generate synthetic data (abundance, endmember signatures, mixed image)
[Y,Yim,A,A_cube,Mth,M_avg]=generate_image(SNR);
P = size(M_avg,2);
[nr,nc,L] = size(Yim);

% extract reference endmember and spectral library from the image
M0 = vca(Y,'Endmembers',P);
M0 = sort_endmembers_to_ref(M_avg,M0);

% spectral library extraction
bundle_nbr = 5; % number of VCA runs
percent = 20; % percentage of pixels considered in each run
Lib = extractbundles_batchVCA(Y, M0, bundle_nbr, percent);



% *********************************************************************** %
% perform unmixing with the different methods
% *********************************************************************** %
clear algnames

% ----------------------------------------
[A_FCLS,M_FCLS,time_fcls,Yhat_FCLS] = adaptor_FCLS(Yim,M0,Lib,[]);
% saves the performance metrics:
algnames{1} = 'FCLS';
RMSEA(end+1) = norm((A-A_FCLS)/numel(A),'fro');
RMSEY(end+1) = norm((Y-Yhat_FCLS)/numel(Y),'fro');
RMSEM(end+1) = nan;
SAMM(end+1)  = nan;
TIMES(end+1) = time_fcls;


% ----------------------------------------
% Set the initialization of some of the other algorithms:
A_init = A_FCLS;



% ----------------------------------------
[A_MESMA,M_MESMA,time_MESMA,Yhat_MESMA] = adaptor_MESMA(Yim,M0,Lib,A_init);
% saves the performance metrics:
algnames{end+1} = 'MESMA';
RMSEA(end+1) = norm((A-A_MESMA)/numel(A),'fro');
RMSEY(end+1) = norm((Y-Yhat_MESMA)/numel(Y),'fro');
RMSEM(end+1) = norm((Mth(:)-M_MESMA(:))/numel(Mth),'fro');
temp = 0; for i=1:nr*nc, for j=1:P, temp = temp + subspace(Mth(:,j,i),M_MESMA(:,j,i))/(nr*nc); end, end; SAMM(end+1) = temp; 
TIMES(end+1) = time_MESMA;



% ----------------------------------------
opt_SocialSparseU.fraction = 1/10;
opt_SocialSparseU.lambda = 0.1;
[A_SocialSparseU,M_SocialSparseU,time_socialSparseU,Yhat_SocialSparseU] = adaptor_SocialSparseU(Yim,M0,Lib,A_init,opt_SocialSparseU);
% saves the performance metrics:
algnames{end+1} = 'Fractional SU';
RMSEA(end+1) = norm((A-A_SocialSparseU)/numel(A),'fro');
RMSEY(end+1) = norm((Y-Yhat_SocialSparseU)/numel(Y),'fro');
RMSEM(end+1) = norm((Mth(:)-M_SocialSparseU(:))/numel(Mth),'fro');
temp = 0; for i=1:nr*nc, for j=1:P, temp = temp + subspace(Mth(:,j,i),M_SocialSparseU(:,j,i))/(nr*nc); end, end; SAMM(end+1) = temp; 
TIMES(end+1) = time_socialSparseU;



% ----------------------------------------
opt_elmm.lambda_s = 1;
opt_elmm.lambda_a = 0.05;
opt_elmm.lambda_psi = 0.01;
[A_ELMM,M_ELMM,time_elmm,Yhat_ELMM] = adaptor_ELMM(Yim,M0,Lib,A_init,opt_elmm);
% saves the performance metrics:
algnames{end+1} = 'ELMM';
RMSEA(end+1) = norm((A-A_ELMM)/numel(A),'fro');
RMSEY(end+1) = norm((Y-Yhat_ELMM)/numel(Y),'fro');
RMSEM(end+1) = norm((Mth(:)-M_ELMM(:))/numel(Mth),'fro');
temp = 0; for i=1:nr*nc, for j=1:P, temp = temp + subspace(Mth(:,j,i),M_ELMM(:,j,i))/(nr*nc); end, end; SAMM(end+1) = temp; 
TIMES(end+1) = time_elmm;



% ----------------------------------------
opt_DeepGUn.dimAut = 2;
opt_DeepGUn.lambda_zref = 0.001;
opt_DeepGUn.lambda_a = 0.01;
[A_DeepGUn,M_DeepGUn,time_DeepGUn,Yhat_DeepGUn] = adaptor_DeepGUn(Yim,M0,Lib,A_init,opt_DeepGUn);
% saves the performance metrics:
algnames{end+1} = 'DeepGUn';
RMSEA(end+1) = norm((A-A_DeepGUn)/numel(A),'fro');
RMSEY(end+1) = norm((Y-Yhat_DeepGUn)/numel(Y),'fro');
RMSEM(end+1) = norm((Mth(:)-M_DeepGUn(:))/numel(Mth),'fro');
temp = 0; for i=1:nr*nc, for j=1:P, temp = temp + subspace(Mth(:,j,i),M_DeepGUn(:,j,i))/(nr*nc); end, end; SAMM(end+1) = temp; 
TIMES(end+1) = time_DeepGUn;



% ----------------------------------------
opt_RUSAL.tau = 0.001;
opt_RUSAL.tau2 = 0.001;
[A_RUSAL,M_RUSAL,time_RUSAL,Yhat_RUSAL] = adaptor_RUSAL(Yim,M0,Lib,A_init,opt_RUSAL);
% saves the performance metrics:
algnames{end+1} = 'RUSAL';
RMSEA(end+1) = norm((A-A_RUSAL)/numel(A),'fro');
RMSEY(end+1) = norm((Y-Yhat_RUSAL)/numel(Y),'fro');
RMSEM(end+1) = nan;
SAMM(end+1)  = nan; 
TIMES(end+1) = time_RUSAL;



% ----------------------------------------
opt_NCM = [];
[A_NCM,M_NCM,time_NCM,Yhat_NCM] = adaptor_NCM(Yim,M0,Lib,A_init,opt_NCM);
% saves the performance metrics:
algnames{end+1} = 'NCM';
RMSEA(end+1) = norm((A-A_NCM)/numel(A),'fro');
RMSEY(end+1) = norm((Y-Yhat_NCM)/numel(Y),'fro');
RMSEM(end+1) = nan;
SAMM(end+1)  = nan; 
TIMES(end+1) = time_NCM;



% ----------------------------------------
[A_BCM,M_BCM,time_BCM,Yhat_BCM] = adaptor_BCM(Yim,M0,Lib,A_init);
% saves the performance metrics:
algnames{end+1} = 'BCM';
RMSEA(end+1) = norm((A-A_BCM)/numel(A),'fro');
RMSEY(end+1) = norm((Y-Yhat_BCM)/numel(Y),'fro');
RMSEM(end+1) = nan;
SAMM(end+1)  = nan; 
TIMES(end+1) = time_BCM;

end


% remove the 'nan' initialization from the vectors
RMSEA = RMSEA(2:end);
RMSEY = RMSEY(2:end);
RMSEM = RMSEM(2:end);
SAMM  = SAMM(2:end);
TIMES = TIMES(2:end);

% reorder
RMSEA = reshape(RMSEA,[],num_runs);
RMSEY = reshape(RMSEY,[],num_runs);
RMSEM = reshape(RMSEM,[],num_runs);
SAMM  = reshape(SAMM, [],num_runs);
TIMES = reshape(TIMES,[],num_runs);


% *********************************************************************** %
% show average results
% *********************************************************************** %
scaleV = 10000;

fprintf('\n\n Abundance estimation results: \n')
for i=1:length(algnames)
    fprintf([pad(algnames{i},15,'right','.') ': %f    +-%f \n'], scaleV*mean(RMSEA(i,:)), scaleV*std(RMSEA(i,:)) )
end

fprintf('\n\n RMSEM results: \n')
for i=1:length(algnames)
    fprintf([pad(algnames{i},15,'right','.') ': %f    +-%f \n'], scaleV*mean(RMSEM(i,:)), scaleV*std(RMSEM(i,:)) )
end

fprintf('\n\n SAMM results: \n')
for i=1:length(algnames)
    fprintf([pad(algnames{i},15,'right','.') ': %f    +-%f \n'], scaleV*mean(SAMM(i,:)), scaleV*std(SAMM(i,:)) )
end

fprintf('\n\n RMSEY results: \n')
for i=1:length(algnames)
    fprintf([pad(algnames{i},15,'right','.') ': %f    +-%f \n'], scaleV*mean(RMSEY(i,:)), scaleV*std(RMSEY(i,:)) )
end

fprintf('\n\n Times results: \n')
for i=1:length(algnames)
    fprintf([pad(algnames{i},15,'right','.') ': %f    +-%f \n'], mean(TIMES(i,:)), std(TIMES(i,:)) )
end




