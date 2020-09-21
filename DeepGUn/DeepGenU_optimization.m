function [EMs_tensor,Z_latent]=DeepGenU_optimization(Z_init,data,A,mappingAEC,Zref,lambda_zref,flag_useparfor)
% -------------------------------------------------------------------------
% Find the endmember representations in the VAE latent space which 
% best represent the HS data, given an abundance matrix
% 
% INPUTS:
%   Z_init         : Initialization of endmember variables in latent space (Lz * P * nr * nc)
%   data           : nr * nc * L  hyeprspectral data cube
%   A              : nr * nc * P  Estimated abundance maps
%   mappingAEC     : Cell array with generative model of each endmember class
%   Zref           : Lz * P  reference endmembers in latent space
%   lambda_zref    : Regularization parameter
%   flag_useparfor : Allows us to disable parfor if necesary
% 
% OUTPUTS:
%   EMs_tensor     : Cell array with estimated endmember tensors 
%   Z_latent       : Optimal EM latent representations
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



% optimization parameters
% % 'OptimalityTolerance',1e-2 --> size of the gradient
% options = optimoptions(@fminunc,'Display','off','Algorithm','quasi-newton',...
%     'OptimalityTolerance',1e-3, 'StepTolerance',1e-3);

% optimOpts1 = struct('GradObj','off','Display','off','LargeScale','off',...
%     'HessUpdate','bfgs','InitialHessType','identity','GoalsExactAchieve',1,...
%     'StoreN',10,'TolX',1e-3,'TolFun',1e-3);

optimOpts2 = struct('Display','off','Method','lbfgs','numDiff',2, ...
    'progTol', 1e-3,'optTol',1e-3,'Corr',25);


% get constants
[nr,nc,L] = size(data);
N = nr*nc;
P = size(A,1);

% Get autoencoder latent space dimension
dimAut = size(mappingAEC{1}{1}.W,2);


% order data as a matrix
Y = reshape(data,nr*nc,L)';
Z = zeros(dimAut,P,nr*nc);


disp('Estimating latent EM representations...'), fprintf('-')
parfor (n=1:N, flag_useparfor)
% for n=1:N
    
    % working...
    % if mod(n,4)==0, fprintf('\b\\'), elseif  mod(n,4)==1, fprintf('\b|'), elseif mod(n,4)==2, fprintf('\b/'), else, fprintf('\b-'), end
    
%     % Define cost function
%     func = @(z)cf_optim_err(z, Y(:,n), A(:,n), mappingAEC, Zref, lambda_zref);
    
%     % Solve opt problem: 950s ------------------------------
%     z0 = Z_init(:,:,n);
%     zout = fminunc(func,z0,options);
%     Z(:,:,n) = zout;

%     % test other solver: 360s ------------------------------
%     z0 = Z_init(:,:,n);
%     zout = fminlbfgs(func, z0, optimOpts1);
%     Z(:,:,n) = zout;
    
    % test other solver: 108s ------------------------------
	func = @(z)cf_optim_err(reshape(z,dimAut,P), Y(:,n), A(:,n), mappingAEC, Zref, lambda_zref);
	z0 = Z_init(:,:,n);
    z0 = z0(:);
    zout = minFunc(func, z0, optimOpts2)';
    zout = reshape(zout, dimAut, P);
    Z(:,:,n) = zout;
end


% Construct EM tensor from the optimized latent representations
EMs_tensor = cell(1,P);
for i=1:P
    EMs_tensor{i} = zeros(L,N);
    for n=1:N
        EMs_tensor{i}(:,n) = real(perform_net_pass(mappingAEC{i}, Z(:,i,n)));
    end
    EMs_tensor{i} = reshape(EMs_tensor{i}',nr,nc,L);
end

Z_latent = Z;
end


% -------------------------------------------------------------------------
function [nerr]=cf_optim_err(z, y,a,mappingAEC, Zref, lambda_zref)

% get constants
P = size(a,1);
L = size(y,1);
M = zeros(L,P);

% compute EM spectra from latent representations
for i=1:P
    M(:,i) = real(perform_net_pass(mappingAEC{i}, z(:,i)));
end

% Compute the cost
nerr = norm(y-M*a)^2 + lambda_zref * norm(z-Zref,'fro')^2;
end




