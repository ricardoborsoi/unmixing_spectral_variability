% This script is an example of use of the algorithm of
%
%   L. Drumetz, M. A. Veganzones, S. Henrot, R. Phlypo, J. Chanussot and 
%   C. Jutten, "Blind Hyperspectral Unmixing Using an Extended Linear
%   Mixing Model to Address Spectral Variability," in IEEE Transactions on 
%   Image Processing, vol. 25, no. 8, pp. 3890-3905, Aug. 2016.
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
% We provide the data cube and reference endmembers. The unmixing is
% performed using the Extended Linear Mixing Model (ELMM), which
% incorporates spectral variability. Several illustrative results are then
% displayed.
%
% Latest Revision: 17-November-2016
% Revision: 1.0
%% DEMO ELMM algorithm

close all
clear

%% load data

        
load real_data_1 

v = [57,30,20]; % RGB bands

[m,n,L] = size(data);

rgb(:,:,1) = imadjust(rescale(data(:,:,v(1)),1));
rgb(:,:,2) = imadjust(rescale(data(:,:,v(2)),1));
rgb(:,:,3) = imadjust(rescale(data(:,:,v(3)),1));
figure, imshow(rgb); % display RGB image

data_r = reshape(data,m*n,L);

load endmembers_houston % load initial reference endmembers
P = size(S0,2);
N = m*n;

materials{1} = 'vegetation'; % identify materials
materials{2} = 'red roofs';
materials{3} = 'concrete';
materials{4} = 'asphalt';



%% Fully Constrained Least Squares Unmixing (FCLSU)

disp('FCLSU...')
tic
A_FCLSU = FCLSU(data_r',S0)';
toc

%% Scaled version of the (partially) Constrained Least Squares (SCLSU)

disp('S-CLSU...')

tic
[A_SCLSU, psis_SCLSU] = SCLSU(data, S0);
toc

S_SCLSU = zeros(L,P,N);

parfor i = 1:m*n
   S_SCLSU(:,:,i) = S0*diag(psis_SCLSU(:,i)); 
end

%% Full Extended Linear Mixing Model

disp('ELMM')

% initialization with S-CLSU

A_init = A_SCLSU;
psis_init = ones(size(A_init));

% tune regularization parameters
% For the regularizations on the abundances and scaling factors, the
% regularization parameters can be either the same for all materials 
% (scalar value), or different for each of them (on the form of a vector).

        lambda_s = 0.5; 
        lambda_a = 0.015;
        lambda_psi = [0.06,0.05,0.06,0.001]; 

% optional parameters

norm = '1,1'; % Use a Total Variation on the abundances
verbose = true; % display

maxiter_anls = 100;
maxiter_admm = 100;
epsilon_s = 10^(-3);
epsilon_a = 10^(-3);
epsilon_psi = 10^(-3);
epsilon_admm_abs = 10^(-2);
epsilon_admm_rel = 10^(-2);

% run

% simple version with default parameters for the stopping criteria

% tic
% [A_ELMM, psis_ELMM, S_ELMM, optim_struct] = ELMM_ADMM(data, A_init, psis_init, S0,lambda_s,lambda_a,lambda_psi,norm,verbose);
% toc

% complete syntax

tic
[A_ELMM, psis_ELMM, S_ELMM, optim_struct] = ELMM_ADMM(data, A_init, psis_init, S0,lambda_s,lambda_a,lambda_psi,norm,verbose,maxiter_anls,maxiter_admm,epsilon_s,epsilon_a,epsilon_psi,epsilon_admm_abs,epsilon_admm_rel);
toc

%% display some results

% reconstruction errors (Root Mean Squared Errors on the reconstruction)

H_FCLSU = S0*A_FCLSU; % reconstruction for S-FCLSU
H_SCLSU= zeros(L,m*n); % reconstruction for S-CLSU
H_ELMM= zeros(L,m*n); % reconstruction for ELMM

parfor i=1:m*n
   H_ELMM(:,i) = squeeze(S_ELMM(:,:,i))*A_ELMM(:,i); 
   H_SCLSU(:,i) = squeeze(S_SCLSU(:,:,i))*A_SCLSU(:,i);
end

  RMSE_FCLSU = sqrt(1/L*sum((H_FCLSU'-data_r).^2,2));
  RMSE_SCLSU = sqrt(1/L*sum((H_SCLSU'-data_r).^2,2));
  RMSE_ELMM = sqrt(1/L*sum((H_ELMM'-data_r).^2,2));
            
  RMSE_FCLSU_im = reshape(RMSE_FCLSU,m,n);
  RMSE_SCLSU_im = reshape(RMSE_SCLSU,m,n);
  RMSE_ELMM_im = reshape(RMSE_ELMM,m,n);
  
  maxRMSE = max([max(RMSE_SCLSU(:)),max(RMSE_ELMM(:))]);
  minRMSE = min([min(RMSE_SCLSU(:)),min(RMSE_ELMM(:))]);
  
  figure, 
  imshow(RMSE_FCLSU_im,[])
  title(['RMSE FCLSU, mean: ',num2str(mean(RMSE_FCLSU(:)))])
  colormap jet
  colorbar
  
  figure
  subplot(1,2,1)
  imshow(RMSE_SCLSU_im,[])
  title(['RMSE SCLSU, mean: ',num2str(mean(RMSE_SCLSU(:)))])
  set(gca,'clim',[minRMSE,maxRMSE])
  colormap jet
  colorbar
  
  subplot(1,2,2)
  imshow(RMSE_ELMM_im,[])
  title(['RMSE ELMM, mean: ',num2str(mean(RMSE_ELMM(:)))])
  set(gca,'clim',[minRMSE,maxRMSE])
  colormap jet
  colorbar
  
  set(gcf,'color', 'white')

   

% abundance maps

A_FCLSU_im = reshape(A_FCLSU',m,n,P);
A_SCLSU_im = reshape(A_SCLSU',m,n,P);
A_ELMM_im = reshape(A_ELMM',m,n,P);

figure,
for p = 1:P
   subplot(3,P,p)
   imshow(A_FCLSU_im(:,:,p),[])
   set(gca,'clim',[0,1])
   if p == 1
   ylabel(' abundances FCLSU')
   end
   colormap jet
   
   subplot(3,P,P+p)
   imshow(A_SCLSU_im(:,:,p),[])
   set(gca,'clim',[0,1])
   if p == 1
   ylabel('abundances SCLSU')
   end
   colormap jet

   subplot(3,P,2*P+p)
   imshow(A_ELMM_im(:,:,p),[])
   set(gca,'clim',[0,1])
   if p == 1
   ylabel('abundances ELMM')
   end
   colormap jet
   xlabel(materials{p})
   
end

set(gcf,'color','white')

% scaling factor maps

psis_SCLSU_im = reshape(psis_SCLSU',m,n,P);
psis_ELMM_im = reshape(psis_ELMM',m,n,P);

figure,
   imshow(psis_SCLSU_im(:,:,1),[])
   title('SCLSU')
   colormap jet
   colorbar  
 ylabel('SCLSU')
 
figure
for p = 1:P
   curr_psis = psis_ELMM_im(:,:,p);
   subplot(2,2,p)
   imshow(psis_ELMM_im(:,:,p),[])
%    set(gca,'clim',[prctile(curr_psis(:),1),prctile(curr_psis(:),99)])
   set(gca,'clim',[quantile(curr_psis(:),0),quantile(curr_psis(:),1)])
 
   colormap jet
   colorbar
   title(materials{p})
   ylabel('scaling factors of the ELMM')
end

set(gcf,'color','white')

% sources (local endmembers)

figure,
for p = 1:P
   subplot(2,2,p)
   plot(squeeze(S_SCLSU(:,p,:)))
   title([materials{p}, ' for SCLSU'])
end
set(gcf,'color','white')

figure,
for p = 1:P
   subplot(2,2,p)
   plot(squeeze(S_ELMM(:,p,:)))
   title(materials{p})
     title([materials{p}, ' for full ELMM'])
end
set(gcf,'color','white')



% scatterplot in 3D for the full ELMM

pca_viz(data_r',reshape(S_ELMM,L,P*N));
