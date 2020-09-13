function [] = pca_viz(data,EMs)
% PCA_VIZ Display 3D scatterplot of the data, along with the endmembers
%
% Given the data matrix and the endmembers, we display the information
% visually by projecting the data and endmembers onto the 3D space formed 
% by the first three principal components of the data.
%
% The data are displayed in blue, and the endmembers in red.
%
% data: L*N data matrix
% EMs: L*nbems matrix containing all the endmembers to display
%
% Author: Lucas Drumetz
% Latest Revision: 17-November-2016
% Revision: 1.3

sigma = cov(data'); % covariance matrix
[U,eigenvalues] = eigs(sigma); % U contains the first eigenvectors and D the largest eigenvalues

projection_matrix = U(:,1:3)*U(:,1:3)'; % projection matrix on the first 3 eigen axes
EM =EMs;


projected_data = projection_matrix*data; %projection of the data on the first 3 components
projected_EM = projection_matrix*EM;

PCs = U(:,1:3)'*projected_data; % same with the data
PCs_EM = U(:,1:3)'*projected_EM;

 scatterplot = figure;
 scatter3(PCs(1,:),PCs(2,:),PCs(3,:),'.')   % scatter plot of the data
   hold on
 scatter3(PCs_EM(1,:),PCs_EM(2,:),PCs_EM(3,:),'r+','linewidth',1)

 set(gca,'fontname','times','fontsize',25)
%  axis equal

xlabel('PC1','fontname','times','fontsize',25)
ylabel('PC2','fontname','times','fontsize',25)
zlabel('PC3','fontname','times','fontsize',25)
set(gcf,'color','white')
