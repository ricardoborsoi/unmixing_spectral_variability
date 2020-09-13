function [mapped_dataAEC, mappingAEC, synthSignaturesMan, Zref]=train_DeepGen_EM_model(pure_pxs, M0, dimAut, AECtype, actFunStr)
% -------------------------------------------------------------------------
% Trains thedeep generative models for the endmembers using the pure pixels 
% extracted from the image
% 
% INPUTS:
%   pure_pxs  : Cell array with P cells. Each cell is an L * K matrix with K 
%               pure pixels to train the generative model (AEC)
%   M0        : Reference/average endmember signatures
%   dimAut    : Dimension of the latent space
%   AECtype   : String, either 'deterministic' or 'VAE'
%   actFunStr : Activation function for the middle network layers (only valid for 'VAE')
% 
% OUTPUTS:
%   mapped_dataAEC     : Cell array with training pure pixels mapped to the latent space
%   mappingAEC         : Cell array with the learned generative/discriminative EM models
%   synthSignaturesMan : Signatures of each class randomly synthesized using the generative model
%   Zref               : Reference EM signatures mapped to the latent space
% 
% -------------------------------------------------------------------------
% This code is part of the DeepGUn algorithm, referring to the following publication:
% 
%   "Deep Generative Endmember Modeling: An Application to Unsupervised Spectral Unmixing"
%   Ricardo Augusto Borsoi, Tales Imbiriba, José Carlos Moreira Bermudez
%   IEEE Transactions on Computational Imaging
% 
% The DeepGUn algorithm performs spectral unmixing with spectral variability
% modeling the endmembers using deep generative models (variational
% autoencoders). 
% -------------------------------------------------------------------------
% Author: Ricardo Borsoi
% last revision: 02/09/2019
%
% -------------------------------------------------------------------------



    % $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
if strcmp(AECtype,'deterministic')
    % -------------------------------------------------------------------------
    % How do I project low-dimensional data back into the data space? 
    % Back-projection can only be implemented for linear techniques, for 
    %     autoencoders, and for the GPLVM. For some of these models, the 
    %     toolbox implements back-projection via the reconstruct_data.m function.
    addpath(genpath('../utils/drtoolbox'))

    % number of additional samples to interpolate
    numSamp = 5;
    spreadFact = 3;
    % dimAut = 2;

    P = length(pure_pxs);
    
    mapped_dataAEC = cell(1,P);
    mappingAEC     = cell(1,P);
    synthSignaturesMan = cell(1,P);

    % learn manifold and interpolate ------------------------------
    for i=1:P
        % [mapped_dataAEC{i}, mappingAEC{i}] = compute_mapping(trainingSmps{i}', 'PCA', dimAut);
        [mapped_dataAEC{i}, mappingAEC{i}] = compute_mapping(pure_pxs{i}', 'Autoencoder', dimAut);

        % generate samples in the latent space
        muu = mean(mapped_dataAEC{i});
        cxx = cov(mapped_dataAEC{i});
        mappedX = spreadFact * randn(numSamp, size(mapped_dataAEC{i},2)) * sqrtm(cxx) + repmat(muu,[numSamp 1]);
        synthSignaturesMan{i} = real(reconstruct_data(mappedX, mappingAEC{i}));
    end

    
    % $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
elseif strcmp(AECtype,'VAE')
    % train a VAE in keras (python)
    P = length(pure_pxs);
    
    for em_idx=1:P
        trainingData = pure_pxs{em_idx}';
        latent_dim = dimAut;
        % batchSize = 20; % batch size for network training
        batchSize = ceil(size(trainingData,1)/3) + 2;
        beta_loss = 1; % scale the KL divergence term in the AEC cost function
        m_idx = M0(:,em_idx);
        save('python/training_EM_data.mat','trainingData','em_idx','latent_dim','batchSize','beta_loss','actFunStr','m_idx')

        % train VAE and interpolate some samples to test
        [~,~]=system('python python/vae_keras_ems_train.py');
    end
    
    % TODO: get mapped training data, and synth signatures
    mapped_dataAEC = [];
    mappingAEC = loadNetworks('vae_EM_idx',P);
    synthSignaturesMan = [];
    
    % Get average EM matrix
    Zref = zeros(dimAut,P);
    for i=1:P
        load(['latent_ref_EM_idx' num2str(i) '.mat'])
        Zref(:,i) = z_mean;
    end
    
    % $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
else
    error('Unknown autoencoder type!') 
end


