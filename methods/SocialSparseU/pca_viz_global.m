function [] = pca_viz_global(data_r,S_FCLSU,S_collaborative,S_group,S_elitist,S_fractional,allEMsorted_rs)

% PCA_VIZ Display 3D scatterplot of the data, along with the endmembers
% for all the algorithms used in the paper
% Given the data matrix and the endmembers, we display the information
% visually by projecting the data and endmembers onto the 3D space formed 
% by the first three principal components of the data.
%
% The data are displayed in blue, and the endmembers in red.
%
% data: L*N data matrix
% S_*****: L*(N*nbems) matrix containing all the endmembers to display for a
% given algorithm
% - allEMsorted_rs: matrix containing the extracted bundles
%
% Author: Lucas Drumetz
% Latest Revision: 22-July-2019
% Revision: 1.2
%

sigma = cov(data_r'); % covariance matrix
[U,eigenvalues] = eigs(sigma); % U contains the first eigenvectors and D the largest eigenvalues

projection_matrix = U(:,1:3)*U(:,1:3)'; % projection matrix on the first 3 eigen axes


projected_data = projection_matrix*data_r; %projection of the data on the first 3 components
projected_FCLSU = projection_matrix*S_FCLSU;
projected_collaborative = projection_matrix*S_collaborative;
projected_group = projection_matrix*S_group;
projected_elitist = projection_matrix*S_elitist;
projected_fractional = projection_matrix*S_fractional;


projected_bundle_var = projection_matrix*allEMsorted_rs;

% projected_EM_bundle = projection_matrix*EMs_bundle;

PCs = U(:,1:3)'*projected_data; % same with the data
PCs_FCLSU = U(:,1:3)'*projected_FCLSU;
PCs_collaborative = U(:,1:3)'*projected_collaborative;
PCs_group = U(:,1:3)'*projected_group ;
PCs_elitist = U(:,1:3)'*projected_elitist;
PCs_fractional = U(:,1:3)'*projected_fractional ;

PCs_bundle_var = U(:,1:3)'*projected_bundle_var ;

%  scatterplot_NCM = figure;

figure,

subaxis(1,5,1,'spacinghorizontal',0.05,'spacingvertical',0.05);
 scatter3(PCs(1,:),PCs(2,:),PCs(3,:),'.')   % scatter plot of the data
   hold on
 scatter3(PCs_FCLSU(1,:),PCs_FCLSU(2,:),PCs_FCLSU(3,:),200,'.','linewidth',2)
hold on
 scatter3(PCs_bundle_var(1,:),PCs_bundle_var(2,:),PCs_bundle_var(3,:),1000,'k.','linewidth',2)
 set(gca,'fontname','times','fontsize',10)
 axis equal
 
xlabel('PC1','fontname','times','fontsize',10)
ylabel('PC2','fontname','times','fontsize',10)
zlabel('PC3','fontname','times','fontsize',10)
title('FCLSU','fontname','times','fontsize',10)
 subaxis(1,5,2,'spacinghorizontal',0.05,'spacingvertical',0.05);
 scatter3(PCs(1,:),PCs(2,:),PCs(3,:),'.')   % scatter plot of the data
   hold on
 scatter3(PCs_collaborative(1,:),PCs_collaborative(2,:),PCs_collaborative(3,:),200,'.','linewidth',2)
hold on
 scatter3(PCs_bundle_var(1,:),PCs_bundle_var(2,:),PCs_bundle_var(3,:),1000,'k.','linewidth',2)
 set(gca,'fontname','times','fontsize',10)
 axis equal
  
xlabel('PC1','fontname','times','fontsize',10)
ylabel('PC2','fontname','times','fontsize',10)
zlabel('PC3','fontname','times','fontsize',10)

title('Collaborative','fontname','times','fontsize',10)

 subaxis(1,5,3,'spacinghorizontal',0.05,'spacingvertical',0.05);
 scatter3(PCs(1,:),PCs(2,:),PCs(3,:),'.')   % scatter plot of the data
   hold on
 scatter3(PCs_group(1,:),PCs_group(2,:),PCs_group(3,:),200,'.','linewidth',2)
hold on
 scatter3(PCs_bundle_var(1,:),PCs_bundle_var(2,:),PCs_bundle_var(3,:),1000,'k.','linewidth',2)
 set(gca,'fontname','times','fontsize',10)
 axis equal
  
xlabel('PC1','fontname','times','fontsize',10)
ylabel('PC2','fontname','times','fontsize',10)
zlabel('PC3','fontname','times','fontsize',10)

title('Group','fontname','times','fontsize',10)

 set(gcf,'color','white')

 subaxis(1,5,4,'spacinghorizontal',0.05,'spacingvertical',0.05);
 scatter3(PCs(1,:),PCs(2,:),PCs(3,:),'.')   % scatter plot of the data
   hold on
 scatter3(PCs_elitist(1,:),PCs_elitist(2,:),PCs_elitist(3,:),200,'.','linewidth',2)
hold on
 scatter3(PCs_bundle_var(1,:),PCs_bundle_var(2,:),PCs_bundle_var(3,:),1000,'k.','linewidth',2)
 set(gca,'fontname','times','fontsize',10)
 axis equal

  
xlabel('PC1','fontname','times','fontsize',10)
ylabel('PC2','fontname','times','fontsize',10)
zlabel('PC3','fontname','times','fontsize',10)

title('Elitist','fontname','times','fontsize',10)

 subaxis(1,5,5,'spacinghorizontal',0.05,'spacingvertical',0.05);
 scatter3(PCs(1,:),PCs(2,:),PCs(3,:),'.')   % scatter plot of the data
   hold on
 scatter3(PCs_fractional(1,:),PCs_fractional(2,:),PCs_fractional(3,:),200,'.','linewidth',2)
hold on
 scatter3(PCs_bundle_var(1,:),PCs_bundle_var(2,:),PCs_bundle_var(3,:),1000,'k.','linewidth',2)
 set(gca,'fontname','times','fontsize',10)
 axis equal
  
xlabel('PC1','fontname','times','fontsize',10)
ylabel('PC2','fontname','times','fontsize',10)
zlabel('PC3','fontname','times','fontsize',10)

title('Fractional','fontname','times','fontsize',10)
set(gcf,'color','white')

end


