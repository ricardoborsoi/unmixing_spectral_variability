% =========================================================================
% This code compares different specral unmixing algorithms that consider
% spectral variability using realistic synthetically generated endmember 
% data.
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

% rng(1000,'twister'); % for reproducibility, if desired

% select the signal to noise ratio:
SNR = 30;

% generate synthetic data (abundance, endmember signatures, mixed image)
[Y,Yim,A,A_cube,Mth,M_avg]=generate_image(SNR);
P = size(M_avg,2);
[nr,nc,L] = size(Yim);

% extract reference endmember and library from the image
M0 = vca(Y,'Endmembers',P);
M0 = sort_endmembers_to_ref(M_avg,M0);

% spectral library extraction
bundle_nbr = 5; % number of VCA runs
percent = 20; % percentage of pixels considered in each run
Lib = extractbundles_batchVCA(Y, M0, bundle_nbr, percent);

% uncomment to plot the extracted spectral library:
%{
figure
subplot(1,3,1)
plot(linspace(0.4,2.45,198), Lib{1})
xlim([0.4 2.45]), ylim([0 0.5])
xlabel('Wavelength [$\mu$m]','interpreter','latex','fontsize',14)
ylabel('Reflectance','interpreter','latex','fontsize',14)
subplot(1,3,2)
plot(linspace(0.4,2.45,198), Lib{2})
xlim([0.4 2.45]), ylim([0 0.9])
xlabel('Wavelength [$\mu$m]','interpreter','latex','fontsize',14)
subplot(1,3,3)
plot(linspace(0.4,2.45,198), Lib{3})
xlim([0.4 2.45]), ylim([0 0.15])
xlabel('Wavelength [$\mu$m]','interpreter','latex','fontsize',14)
%}



%%
% *********************************************************************** %
% perform unmixing with the different algorithms
% *********************************************************************** %

% ----------------------------------------
% FCLS (baseline)
[A_FCLS,M_FCLS,time_fcls,Yhat_FCLS] = adaptor_FCLS(Yim,M0,Lib,[]);
% saves the performance metrics:
algnames{1} = 'FCLS';
RMSEA(1) = norm((A-A_FCLS)/numel(A),'fro');
RMSEY(1) = norm((Y-Yhat_FCLS)/numel(Y),'fro');
RMSEM(1) = nan;
SAMM(1)  = nan;
TIMES(1) = time_fcls;
A_est{1} = A_FCLS;
M_est{1} = nan;


% ----------------------------------------
% Set the initialization of some of the other algorithms:
A_init = A_FCLS;



% ----------------------------------------
% MESMA
[A_MESMA,M_MESMA,time_MESMA,Yhat_MESMA] = adaptor_MESMA(Yim,M0,Lib,A_init);
% saves the performance metrics:
algnames{end+1} = 'MESMA';
RMSEA(end+1) = norm((A-A_MESMA)/numel(A),'fro');
RMSEY(end+1) = norm((Y-Yhat_MESMA)/numel(Y),'fro');
RMSEM(end+1) = norm((Mth(:)-M_MESMA(:))/numel(Mth),'fro');
temp = 0; for i=1:nr*nc, for j=1:P, temp = temp + subspace(Mth(:,j,i),M_MESMA(:,j,i))/(nr*nc); end, end; SAMM(end+1) = temp; 
TIMES(end+1) = time_MESMA;
A_est{end+1} = A_MESMA;
M_est{end+1} = M_MESMA;



% ----------------------------------------
% Fractional Sparse SU
opt_SocialSparseU.fraction = 1/10;
opt_SocialSparseU.lambda = 0.1;
[A_SocialSparseU,M_SocialSparseU,time_socialSparseU,Yhat_SocialSparseU] = adaptor_SocialSparseU(Yim,M0,Lib,A_init,opt_SocialSparseU);
% saves the performance metrics:
algnames{end+1} = 'Fractional';
RMSEA(end+1) = norm((A-A_SocialSparseU)/numel(A),'fro');
RMSEY(end+1) = norm((Y-Yhat_SocialSparseU)/numel(Y),'fro');
RMSEM(end+1) = norm((Mth(:)-M_SocialSparseU(:))/numel(Mth),'fro');
temp = 0; for i=1:nr*nc, for j=1:P, temp = temp + subspace(Mth(:,j,i),M_SocialSparseU(:,j,i))/(nr*nc); end, end; SAMM(end+1) = temp; 
TIMES(end+1) = time_socialSparseU;
A_est{end+1} = A_SocialSparseU;
M_est{end+1} = M_SocialSparseU;



% ----------------------------------------
% ELMM
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
A_est{end+1} = A_ELMM;
M_est{end+1} = M_ELMM;



% ----------------------------------------
% DeepGUn
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
A_est{end+1} = A_DeepGUn;
M_est{end+1} = M_DeepGUn;



% ----------------------------------------
% RUSAL
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
A_est{end+1} = A_RUSAL;
M_est{end+1} = nan;



% ----------------------------------------
% Normal Compositional Model
opt_NCM = [];
[A_NCM,M_NCM,time_NCM,Yhat_NCM] = adaptor_NCM(Yim,M0,Lib,A_init,opt_NCM);
% saves the performance metrics:
algnames{end+1} = 'NCM';
RMSEA(end+1) = norm((A-A_NCM)/numel(A),'fro');
RMSEY(end+1) = norm((Y-Yhat_NCM)/numel(Y),'fro');
RMSEM(end+1) = nan;
SAMM(end+1)  = nan; 
TIMES(end+1) = time_NCM;
A_est{end+1} = A_NCM;
M_est{end+1} = nan;



% ----------------------------------------
% Beta Compositional Model
[A_BCM,M_BCM,time_BCM,Yhat_BCM] = adaptor_BCM(Yim,M0,Lib,A_init);
% saves the performance metrics:
algnames{end+1} = 'BCM';
RMSEA(end+1) = norm((A-A_BCM)/numel(A),'fro');
RMSEY(end+1) = norm((Y-Yhat_BCM)/numel(Y),'fro');
RMSEM(end+1) = nan;
SAMM(end+1)  = nan; 
TIMES(end+1) = time_BCM;
A_est{end+1} = A_BCM;
M_est{end+1} = nan;



%%
% *********************************************************************** %
% show quantitative results
% *********************************************************************** %
scaleV = 10000; % scale the numbers to make vizualization easier

fprintf('\n\n Abundance estimation results: \n')
for i=1:length(algnames)
    fprintf([pad(algnames{i},15,'right','.') ': %f \n'], scaleV*RMSEA(i))
end

fprintf('\n\n RMSY results: \n')
for i=1:length(algnames)
    fprintf([pad(algnames{i},15,'right','.') ': %f \n'], scaleV*RMSEY(i))
end

fprintf('\n\n RMSM results: \n')
for i=1:length(algnames)
    fprintf([pad(algnames{i},15,'right','.') ': %f \n'], scaleV*RMSEM(i))
end

fprintf('\n\n SAMM results: \n')
for i=1:length(algnames)
    fprintf([pad(algnames{i},15,'right','.') ': %f \n'], scaleV*SAMM(i))
end

fprintf('\n\n Times results: \n')
for i=1:length(algnames)
    fprintf([pad(algnames{i},15,'right','.') ': %f \n'], TIMES(i))
end




%%
% *********************************************************************** %
% show visual results
% *********************************************************************** %
FSize = 14; % fontsize for the labels

% plot the abundances ---------------------------------
fh = figure;
[ha, pos] = tight_subplot(length(algnames)+1, 3, 0.01, 0.1, 0.1);

axes(ha(1));
imagesc(A_cube(:,:,1),[0 1]), set(gca,'ytick',[],'xtick',[])%, axis square
axes(ha(2));
imagesc(A_cube(:,:,2),[0 1]), set(gca,'ytick',[],'xtick',[])%, axis square
axes(ha(3));
imagesc(A_cube(:,:,3),[0 1]), set(gca,'ytick',[],'xtick',[])%, axis square

for i=1:length(algnames)
    axes(ha(1 + (i-0)*3));
    imagesc(reshape(A_est{i}(1,:),[nr nc]),[0 1]), set(gca,'ytick',[],'xtick',[])%, axis square
    axes(ha(2 + (i-0)*3));
    imagesc(reshape(A_est{i}(2,:),[nr nc]),[0 1]), set(gca,'ytick',[],'xtick',[])%, axis square
    axes(ha(3 + (i-0)*3));
    imagesc(reshape(A_est{i}(3,:),[nr nc]),[0 1]), set(gca,'ytick',[],'xtick',[])%, axis square
    
    axes(ha(1 + (i-0)*3));
    ylabel(algnames{i},'interpreter','latex','FontSize',FSize)
end

axes(ha(1));
ylabel('Reference','interpreter','latex','FontSize',FSize)
title('Vegetation','interpreter','latex','FontSize',FSize)
axes(ha(2));
title('Soil','interpreter','latex','FontSize',FSize)
axes(ha(3));
title('Water','interpreter','latex','FontSize',FSize)

colormap(jet)

set(fh, 'Position', [0 0 500 800])

% uncomment to save the figure:
% savefig('RESULTS/abundTests.fig')
% print('RESULTS/abundTests.pdf','-dpdf')
% [~,~]=system(['pdfcrop RESULTS/abundTests.pdf RESULTS/abundTests.pdf']);




% plot the estimayed endmembers  -----------------------------
K_decim_spec = 10;
wavelenths = linspace(0.4,2.4,L);
Nplots = 1;
idxM   = [];
for i=1:length(algnames)
    if ~isnan(M_est{i}), Nplots = Nplots + 1; idxM = [idxM i]; end
end

ff = figure;
%[ha, pos] = tight_subplot(Nplots, 3, 0.01, 0.1, 0.1);
%[ha, pos] = tight_subplot(Nplots, 3, 0.05, 0.1, 0.1);
[ha, pos] = tight_subplot(Nplots, 3, [0.01 0.06], 0.1, 0.1);

for j=1:3
    axes(ha(j));
    plot(wavelenths,squeeze(Mth(:,j,1:K_decim_spec:end)))
    xlim([min(wavelenths) max(wavelenths)])
    ylim([0 1.15*max(max(squeeze(Mth(:,j,1:K_decim_spec:end))))])
    set(gca,'xtick',[])
end

for i=1:Nplots-1
    for j=1:3
        axes(ha(j + (i-0)*3));
        plot(wavelenths,squeeze(M_est{idxM(i)}(:,j,1:K_decim_spec:end)))
        xlim([min(wavelenths) max(wavelenths)])
        ylim([0 1.15*max(max(squeeze(M_est{idxM(i)}(:,j,1:K_decim_spec:end))))])
        if i<Nplots-1, set(gca,'xtick',[]), end
    end 
    axes(ha(1 + (i-0)*3));
    ylabel(algnames{idxM(i)},'interpreter','latex','FontSize',FSize)
end

axes(ha(1));
ylabel('Reference','interpreter','latex','FontSize',FSize)
title('Vegetation','interpreter','latex','FontSize',FSize)
axes(ha(2));
title('Soil','interpreter','latex','FontSize',FSize)
axes(ha(3));
title('Water','interpreter','latex','FontSize',FSize)

axes(ha(end-2))
xlabel('Wavelength [$\mu$m]','interpreter','latex','FontSize',FSize)
axes(ha(end-1))
xlabel('Wavelength [$\mu$m]','interpreter','latex','FontSize',FSize)
axes(ha(end))
xlabel('Wavelength [$\mu$m]','interpreter','latex','FontSize',FSize)

colormap(jet)

set(ff, 'Position', [0 0 500 600])

% uncomment to save the figure:
% savefig('RESULTS/estimEMsTests.fig')
% print('RESULTS/estimEMsTests.pdf','-dpdf','-painters')
% [~,~]=system(['pdfcrop RESULTS/estimEMsTests.pdf RESULTS/estimEMsTests.pdf']);



