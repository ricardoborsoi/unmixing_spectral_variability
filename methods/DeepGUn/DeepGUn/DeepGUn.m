function [A_deepGen,M_DeepGen]=DeepGUn(r_cube, A_init, M0, dimAut, lambda_zref, lambda_a, flag_Npx, vec_Npx, flag_useparfor)
% -------------------------------------------------------------------------
% DeepGUn algorithm, referring to the following publication:
% 
%   "Deep Generative Endmember Modeling: An Application to Unsupervised Spectral Unmixing"
%   Ricardo Augusto Borsoi, Tales Imbiriba, José Carlos Moreira Bermudez
%   IEEE Transactions on Computational Imaging
% 
% The DeepGUn algorithm performs spectral unmixing with spectral variability
% modeling the endmembers using deep generative models (variational
% autoencoders). 
% 
% 
% INPUTS:
%   data           : nr * nc * L  hyperspectral image cube
%   A_init         : nr * nc * P  abundance maps initialization
%   M0             : Reference endmember matrix
%   dimAut         : Dimension of the VAEs latent space (recomended = 2)
%   lambda_zref    : Regularization parameter
%   lambda_a       : Regularization parameter
%   flag_Npx       : Flag - if "true" extracts as EM bundles the number of
%                            pure pixels specified in "vec_Npx" that are
%                            closest to M0 
%                         - if "false" extracts as EM bundles all pixels 
%                            whose spectral angles to M0 are smaller than 
%                            those specified in M0
%   vec_Npx        : 1 * P vector containing the specification (i.e. the number
%                     of pure pixels if flag_Npx=true or the spectral angles
%                     if flag_Npx=false) of the pure pixels to be extracted
%                     to construct the bundles, for each material
%   flag_useparfor : Allows us to disable parfor if necesary
% 
% OUTPUTS:
%   A_deepGen : Estimated abundances
%   M_DeepGen : Estimated endmember tensor
% 
% -------------------------------------------------------------------------
% Author: Ricardo Borsoi
% last revision: 02/09/2019
%
% -------------------------------------------------------------------------


%% Extract bundles from the image
% flag_Npx = true;
% vec_Npx = [100 100 100];
[bundleLibs,avg_M,PPidx,EM_pix,IDX_comp] = extract_bundles_by_angle(r_cube,M0,vec_Npx,flag_Npx);


save('extracted_bundles.mat','bundleLibs')


%% Train VAEs with the extracted pure pixels
actFunStr = 'relu'; % network activation function
[mapped_dataAEC, mappingAEC, synthSignaturesMan, Zref]=train_DeepGen_EM_model(bundleLibs, M0, dimAut, 'VAE', actFunStr);


%% Solve the unmixing problem with parametrized endmembers

[A_deepGen,M_DeepGen] = DeepGenU_als_vae(r_cube, A_init, mappingAEC, dimAut, Zref, lambda_zref, lambda_a, flag_useparfor);



