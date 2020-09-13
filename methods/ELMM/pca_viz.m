function [] = pca_viz(data,EM)
% PCA_VIZ Display 3D scatterplot of the data, along with the endmembers
%
% Given the data matrix and the endmembers, we display the information
% visually by projecting the data and endmembers onto the 3D space formed 
% by the first three principal components of the data.
%
% The data are displayed in blue, and the endmembers in red.
%
% data: L*N data matrix
% EM: L*nbems matrix containing all the endmembers to display
%
% Latest Revision: 17-November-2016
% Revision: 1.3

sigma = cov(data'); % covariance matrix
[U,~] = eigs(sigma); % U contains the first eigenvectors

PCs = U(:,1:3)'*data; % orthogonal projection of the data on the 3D space 
% (the columns of U(:,1:3) form an orthonormal basis of the space spanned by the first 3 PCs)
PCs_EM = U(:,1:3)'*EM; % orthogonal projection of the endmembers on the 3D space

figure
scatter3(PCs(1,:),PCs(2,:),PCs(3,:),'+')   % scatter plot of the data
hold on
scatter3(PCs_EM(1,:),PCs_EM(2,:),PCs_EM(3,:),'r.','linewidth',1)

 set(gca,'fontname','times','fontsize',25)
 
%   axis equal % use this for orthonormal axes

xlabel('PC1','fontname','times','fontsize',25)
ylabel('PC2','fontname','times','fontsize',25)
zlabel('PC3','fontname','times','fontsize',25)
set(gcf,'color','white')
