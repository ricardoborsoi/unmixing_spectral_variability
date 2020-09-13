function [ groups, components] = batchvca( pixels, p, bundles, percent)

% Function performing VCA several times on randomly sampled subsets of a
% dataset, and clustering the obtained signatures into bundles using
% kmeans with the spectral angle as a similarity measure.

% Reference:

% Somers, B., Zortea, M., Plaza, A., & Asner, G. P. (2012). Automated 
% extraction of image-based endmember bundles for improved spectral 
% unmixing. IEEE Journal of Selected Topics in Applied Earth Observations
% and Remote Sensing, 5(2), 396-408.

% Authors: Lucas Drumetz & Travis R. Meyer
% Latest Revision: 22-July-2019
% Revision: 1.2

% pixels: data pixels, L*N
% p: number of endmember classes
% bundles: number of subsets of the data
% percent: percentage of the data taken each time
% clustering: either 'kmeans' or 'spectral_clustering'
% sampling is without replacement within one bundle the same pixel can be
% sampled in several bundles.

    B = [];

    m = percent/100 * size(pixels,2);
    runs = 1;
    pixels_update = pixels;
    
    for b = 1:bundles
        [C,I] = datasample(pixels_update,floor(m),2,'Replace',false);
        B = [B, EIA_VCA(C, runs, p,false)];
        pixels_update(:,I) = [];
    end

    components = [];
    for i = 1:bundles
        components = [components, B(i).E];
    end
    
    % clustering part
   
        
    groups = kmeans(components',p,'distance','cosine');
        
  
end
