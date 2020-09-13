% This script reproduces the results of the algorithms presented in
%
%  Drumetz, L., Meyer, T. R., Chanussot, J., Bertozzi, A. L., & Jutten, C.
%  (2019). Hyperspectral image unmixing with endmember bundles and group
%  sparsity inducing mixed norms. IEEE Transactions on Image Processing,
%  28(7), 3435-3450.
%
% to unmix a real dataset acquired above the University of Houston Campus,
% in June 2012. We use a small part of the hyperspectral dataset provided
% courtesy of the Hyperspectral Image Analysis group and the NSF Funded
% Center for Airborne Laser Mapping at the University of Houston, and used
% in the 2013 Data Fusion Contest (DFC):
%
%  C. Debes et al., "Hyperspectral and LiDAR Data Fusion: Outcome of the
%  2013 GRSS Data Fusion Contest," in IEEE Journal of Selected Topics in
%  Applied Earth Observations and Remote Sensing, vol. 7, no. 6,
%  pp. 2405-2418, June 2014.
%
% Details on the dataset can also be found at:
% http://hyperspectral.ee.uh.edu/?page_id=459
%
% The blind unmixing is
% performed using the group sparsity inducing norms tested in the paper
% above, which incorporates spectral variability.
% Several illustrative results are then displayed and quantitative metrics
% are computed.
%
% You can choose to extract the endmembers yourself or to load the ones
% used in the paper (flag 'extract_ems' -- default value = 0). In the
% former case, the values of the regularization parameters may not be
% optimal, and the materials in the visualization part may not correspond
% to what is indicated in this code. 5 instances for each material will be 
% randomly selected for display. To reproduce exactly the results of the 
% paper, choose extract_ems = 0. If you choose extract_ems = 1, and there
% is a bundle with less than 5 signatures in it, the visual results
% are not displayed (but you can still do it on your own if you want).
%
%
% Author: Lucas Drumetz
% Latest Revision: 26-July-2019
% Revision: 1.3

%% DEMO group sparsity inducing norms

close all
clear

extract_ems = 0; % flag to perform endmember bundle extraction (1)
% or to load the results to reproduce the results of the paper.
%% load data

load real_data_1

EMs = figure;
imshow(data(:,:,[57,30,20]))
hold on

[m,n,L] = size(data);

%% Extract endmember candidates and cluster them into bundles

X = reshape(data,m*n,L)';

P = 5; % number of endmembers

bundle_nbr = 10; % number of VCA runs
percent = 10; % percentage of pixels considered in each run
clustering = 'kmeans';

if extract_ems == 1
    seed = 100; % fix seed
    rng(seed)
    [groups, bundle] = batchvca(X, P, bundle_nbr, percent); % extract endmembers and cluster them into groups
else
    load bundles
end


pca_viz(X,bundle) % display dataset scatterplot and extracted endmembers on first three PCA axes


%% unmix (empirical tuning of regularization parameters)

% regular FCLSU on the bundles

tic
disp('FLCSU bundle')
A_FCLSU_bundle = FCLSU(X,bundle)';
toc



%% initialize params

A_init = A_FCLSU_bundle; % abundance initialization for other algorithms
rho = 10;
tol_a = 10^(-6);
maxiter_ADMM = 1000;
fraction = 1/10;

verbose = 0; % display
%% Group penalty

lambda = 2;

tic
disp('group')
type = 'group';
[A_group, optim_struct_group] = social_unmixing(X,bundle,groups,A_init,lambda,rho,maxiter_ADMM,type,fraction,tol_a,verbose);
toc


%% elitist penalty

lambda = 0.5;

tic
disp('elitist')
type = 'elitist';
[A_elitist, optim_struct_elitist] = social_unmixing(X,bundle,groups,A_init,lambda,rho,maxiter_ADMM,type,fraction,tol_a,verbose);
toc


%% Fractional penalty

fraction = 1/10;
lambda = 0.4;

tic
disp('fractional')
type = 'fractional';
[A_fractional, optim_struct_fractional] = social_unmixing(X,bundle,groups,A_init,lambda,rho,maxiter_ADMM,type,fraction,tol_a,verbose);
toc



%% Collaborative sparsity

lambda = 1;
tic
disp('collaborative')
type = 'asc';
tic
[A_collaborative] = ADMM_collaborative_unmixing(X,A_init,bundle,lambda,rho,maxiter_ADMM,type,tol_a,verbose);
toc

%% sum the abundances within each class (bundle2global function)

[A_FCLSU_bundle_final,S_FCLSU_bundle] = bundle2global(A_FCLSU_bundle,bundle,groups);
[A_collaborative_final, S_collaborative] = bundle2global(A_collaborative,bundle,groups);
[A_group_final,S_group] = bundle2global(A_group,bundle,groups);
[A_elitist_final,S_elitist] = bundle2global(A_elitist,bundle,groups);
[A_fractional_final,S_fractional] = bundle2global(A_fractional,bundle,groups);

%% Compute reconstruction errors

H_FCLSU_bundle = bundle*A_FCLSU_bundle; % reconstruction for FCLSU
H_group = bundle*A_group; % reconstruction for FCLSU
H_elitist = bundle*A_elitist; % reconstruction for FCLSU
H_fractional = bundle*A_fractional; % reconstruction for FCLSU
H_collaborative = bundle*A_collaborative;

RMSE_FCLSU_bundle = sqrt(1/L*sum((H_FCLSU_bundle-X).^2,1));
RMSE_group = sqrt(1/L*sum((H_group-X).^2,1));
RMSE_elitist = sqrt(1/L*sum((H_elitist-X).^2,1));
RMSE_fractional = sqrt(1/L*sum((H_fractional-X).^2,1));
RMSE_collaborative =  sqrt(1/L*sum((H_collaborative-X).^2,1));

N = m*n;

SAM_FCLSU = zeros(N,1);
SAM_group = zeros(N,1);
SAM_elitist = zeros(N,1);
SAM_fractional = zeros(N,1);
SAM_collaborative = zeros(N,1);

for k = 1:N
    SAM_FCLSU(k) = 180/pi*real(acos((X(:,k)'*H_FCLSU_bundle(:,k))...
        /(norm(X(:,k))*norm(H_FCLSU_bundle(:,k)))));
end

for k = 1:N
    SAM_group(k) = 180/pi*real(acos((X(:,k)'*H_group(:,k))...
        /(norm(X(:,k))*norm(H_group(:,k)))));
end

for k = 1:N
    SAM_elitist(k) = 180/pi*real(acos((X(:,k)'*H_elitist(:,k))...
        /(norm(X(:,k))*norm(H_elitist(:,k)))));
end

for k = 1:N
    SAM_fractional(k) = 180/pi*real(acos((X(:,k)'*H_fractional(:,k))...
        /(norm(X(:,k))*norm(H_fractional(:,k)))));
end

for k = 1:N
    SAM_collaborative(k) = 180/pi*real(acos((X(:,k)'*H_collaborative(:,k))...
        /(norm(X(:,k))*norm(H_fractional(:,k)))));
end


%% display RMSEs and SAMs

mean(SAM_FCLSU(:));
mean(SAM_collaborative(:)) ;
mean(SAM_group(:));
mean(SAM_elitist(:));
mean(SAM_fractional(:));

mean(RMSE_FCLSU_bundle(:));
mean(RMSE_collaborative(:));
mean(RMSE_group(:));
mean(RMSE_elitist(:));
mean(RMSE_fractional(:));

%% Global abundance maps

A_FCLSU_bundle_im = reshape(A_FCLSU_bundle_final',m,n,P);
A_group_im = reshape(A_group_final',m,n,P);
A_elitist_im = reshape(A_elitist_final',m,n,P);
A_fractional_im = reshape(A_fractional_final',m,n,P);
A_collaborative_im =reshape(A_collaborative_final',m,n,P);

figure,
for p = 1:P
    
    subaxis(5,P,p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
    imshow(A_FCLSU_bundle_im(:,:,p),[],'colormap', jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('FCLSU','fontname','times','fontsize',15)
    end
    colormap jet
    
    subaxis(5,P,P+p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
    imshow(A_collaborative_im(:,:,p),[],'colormap', jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('Collaborative','fontname','times','fontsize',15)
    end
    colormap jet
    
    subaxis(5,P,2*P+p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
    imshow(A_group_im(:,:,p),[],'colormap', jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('Group','fontname','times','fontsize',15)
    end
    colormap jet
    
    subaxis(5,P,3*P+p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
    imshow(A_elitist_im(:,:,p),[],'colormap', jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('Elitist','fontname','times','fontsize',15)
    end
    colormap jet
    
    subaxis(5,P,4*P+p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
    imshow(A_fractional_im(:,:,p),[],'colormap', jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('Fractional','fontname','times','fontsize',15)
        xlabel('Concrete','fontname','times','fontsize',15)
    elseif p == 2
        xlabel('Red Roofs','fontname','times','fontsize',15)
    elseif p == 3
        xlabel('Vegetation','fontname','times','fontsize',15)
    elseif p == 4
        xlabel('Asphalt','fontname','times','fontsize',15)
    else
        xlabel('Colored Structures','fontname','times','fontsize',15)
    end
    
    
    
end
set(gcf,'color', 'white')

%% Show a few endmember variants abundance maps

A_FCLSU_bundle_full_im = reshape(A_FCLSU_bundle',m,n,P*bundle_nbr);
A_group_full_im = reshape(A_group',m,n,P*bundle_nbr);
A_elitist_full_im = reshape(A_elitist',m,n,P*bundle_nbr);
A_fractional_full_im = reshape(A_fractional',m,n,P*bundle_nbr);
A_collaborative_full_im = reshape(A_collaborative',m,n,P*bundle_nbr);


%% Vegetation

EMs_cons = find(groups == 4); % find corresponding group

if extract_ems == 0
    ems_sel = [1,3,6,7,9]; % select instances to display
elseif sum(groups == 4) > 4
    ems_sel = randi(sum(groups == 4),5);
else
    disp('less than 5 elements in bundle -- results are not displayed')
    return
end

figure,
for p = 1:length(ems_sel)
    
    
    subaxis(5,length(ems_sel),p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
    imshow(A_FCLSU_bundle_full_im(:,:,EMs_cons(ems_sel(p))),[],'colormap',jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('FCLSU','fontname','times','fontsize',15)
    end
    
    subaxis(5,length(ems_sel),length(ems_sel)+p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
    imshow(A_collaborative_full_im(:,:,EMs_cons(ems_sel(p))),[],'colormap',jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('Collaborative','fontname','times','fontsize',15)
    end
    
    subaxis(5,length(ems_sel),2*length(ems_sel)+p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
    imshow(A_group_full_im(:,:,EMs_cons(ems_sel(p))),[],'colormap',jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('Group','fontname','times','fontsize',15)
    end
    
    
    subaxis(5,length(ems_sel),3*length(ems_sel)+p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
    imshow(A_elitist_full_im(:,:,EMs_cons(ems_sel(p))),[],'colormap',jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('Elitist','fontname','times','fontsize',15)
    end
    
    
    subaxis(5,length(ems_sel),4*length(ems_sel)+p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
    imshow(A_fractional_full_im(:,:,EMs_cons(ems_sel(p))),[],'colormap',jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('Fractional','fontname','times','fontsize',15)
    end
    set(gcf,'colormap', jet)
end
set(gcf,'color', 'white')


%% Red roofs

EMs_cons = find(groups == 2); % find corresponding group

if extract_ems == 0
    ems_sel = [1,3,5,9,4];  % select instances to display
elseif sum(groups == 4) > 4
    ems_sel = randi(sum(groups == 2),5);
else
    disp('less than 5 elements in bundle -- results are not displayed')
    return
end

figure,
for p = 1:length(ems_sel)
    
    
    subaxis(5,length(ems_sel),p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
    imshow(A_FCLSU_bundle_full_im(:,:,EMs_cons(ems_sel(p))),[],'colormap',jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('FCLSU','fontname','times','fontsize',15)
    end
    
    subaxis(5,length(ems_sel),length(ems_sel)+p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
    imshow(A_collaborative_full_im(:,:,EMs_cons(ems_sel(p))),[],'colormap',jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('Collaborative','fontname','times','fontsize',15)
    end
    
    subaxis(5,length(ems_sel),2*length(ems_sel)+p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
    imshow(A_group_full_im(:,:,EMs_cons(ems_sel(p))),[],'colormap',jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('Group','fontname','times','fontsize',15)
    end
    
    
    subaxis(5,length(ems_sel),3*length(ems_sel)+p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
    imshow(A_elitist_full_im(:,:,EMs_cons(ems_sel(p))),[],'colormap',jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('Elitist','fontname','times','fontsize',15)
    end
    
    
    subaxis(5,length(ems_sel),4*length(ems_sel)+p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
    imshow(A_fractional_full_im(:,:,EMs_cons(ems_sel(p))),[],'colormap',jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('Fractional','fontname','times','fontsize',15)
    end
    set(gcf,'colormap', jet)
end
set(gcf,'color', 'white')

%% Concrete

EMs_cons = find(groups == 1); % find corresponding group

if extract_ems == 0
    ems_sel = [1,4,6,10,11]; % select instances to display
elseif sum(groups == 4) > 4
    ems_sel = randi(sum(groups == 1),5);
else
    disp('less than 5 elements in bundle -- results are not displayed')
    return
end

figure,
for p = 1:length(ems_sel)
    
    
    subaxis(5,length(ems_sel),p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
    imshow(A_FCLSU_bundle_full_im(:,:,EMs_cons(ems_sel(p))),[],'colormap',jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('FCLSU','fontname','times','fontsize',15)
    end
    
    subaxis(5,length(ems_sel),length(ems_sel)+p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
    imshow(A_collaborative_full_im(:,:,EMs_cons(ems_sel(p))),[],'colormap',jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('Collaborative','fontname','times','fontsize',15)
    end
    
    subaxis(5,length(ems_sel),2*length(ems_sel)+p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
    imshow(A_group_full_im(:,:,EMs_cons(ems_sel(p))),[],'colormap',jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('Group','fontname','times','fontsize',15)
    end
    
    
    subaxis(5,length(ems_sel),3*length(ems_sel)+p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
    imshow(A_elitist_full_im(:,:,EMs_cons(ems_sel(p))),[],'colormap',jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('Elitist','fontname','times','fontsize',15)
    end
    
    
    subaxis(5,length(ems_sel),4*length(ems_sel)+p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
    imshow(A_fractional_full_im(:,:,EMs_cons(ems_sel(p))),[],'colormap',jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('Fractional','fontname','times','fontsize',15)
    end
    set(gcf,'colormap', jet)
end
set(gcf,'color', 'white')



%% Show endmember distributions in first three PCA axes

pca_viz_global(X,reshape(S_FCLSU_bundle,L,P*N),reshape(S_collaborative,L,P*N),reshape(S_group,L,P*N),reshape(S_elitist,L,P*N),reshape(S_fractional,L,P*N),bundle)


 

