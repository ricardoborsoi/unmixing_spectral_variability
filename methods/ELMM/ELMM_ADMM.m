function [A,psi_maps,S,optim_struct] = ELMM_ADMM(data, A_init, psis_init, S0,lambda_s,lambda_a,lambda_psi,varargin)
%ELMM_ADMM Unmix hyperspectral data using the Extended Linear Mixing Model
%(ELMM).
%
%   We find a stationary point of the following functional:
%   J(S,A,PSI) = 1/2 * sum_{k=1}^{N} (||x_k - S_k*a_k||_{2}^{2} + 
%   ||S_k - S0*psi_{k}||_{F}^{2}) + lambda_A R(A) + lambda_PSI R(PSI)
%
%   with S a collection of endmember matrices for each pixel, A the
%   abundances in each pixel and for each endmember, and PSI the scaling
%   factors advocated by the ELMM.
%
%   The abundances are subject to the usual nonnegativity and sum to one
%   constraints. The scaling factors and endmember spectral are nonnegative
%   as well.
%
%   R(A) is a spatial regularization term on the abundances. It can be 
%   either an anisotropic total variation term TV(A) applied on each 
%   material or a Tikhonov like regularization on the spatial gradients of 
%   the abundance maps. R(PSI) is a differentiable Tikhonov regularization 
%   on the spatial gradients of the scaling factor maps.
%
%   Mandatory inputs:
%   -data: m*n*L image cube, where m is the number of rows, n the number of
%   columns, and L the number of spectral bands.
%   -A_init: P*N initial abundance matrix, with P the number of endmembers
%   to consider, and N the number of pixels (N=m*n)
%   -psis_init: P*N initial scaling factor matrix
%   -S0: L*P reference endmember matrix
%   -lambda_s: regularization parameter on the ELMM tightness
%   -lambda_a: regularization parameter for the spatial regularization on
%   the abundances.
%   -lambda_psi: regularization parameter for the spatial regularization on
%   the scaling factors
%   The spatial regularization parameters can be scalars, in which case 
%   they will apply in the same way for all the terms of the concerned
%   regularizations. If they are vectors, then each term of the sum
%   (corresponding to each material) will be differently weighted by the
%   different entries of the vector.
%
%   Optional inputs (arguments are to be provided in the same order as in 
%   the following list):
%   -norm_sr: choose norm to use for the spatial regularization on the
%   abundances. Can be '2,1' (Tikhonov like penalty on the gradient) or
%   '1,1' (Total Variation) (default: '1,1')
%   -verbose: flag for display in console. Display if true, no display
%   otherwise (default: true)
%   -maxiter_anls: maximum number of iterations for the ANLS loop (default:
%   100)
%   -maxiter_admm: maximum number of iterations for the ADMM loop (default:
%   100)
%   -epsilon_s: tolerance on the relative variation of S between two
%   consecutive iterations (default: 10^(-3))
%   -epsilon_a: tolerance on the relative variation of A between two
%   consecutive iterations (default: 10^(-3))
%   -epsilon_psi: tolerance on the relative variation of psi between two
%   consecutive iterations (default: 10^(-3))
%   -epsilon_admm_abs: tolerance on the absolute part of the primal and
%   dual residuals (default: 10^(-2))
%   -epsilon_admm_rel: tolerance on the relative part of the primal and
%   dual residuals (default: 10^(-2))
%
%   Outputs:
%   -A: P*N abundance matrix
%   -psi_maps: P*N scaling factor matrix
%   -S: L*P*N tensor constaining all the endmember matrices for each pixel
%   -optim_struct: structure containing the values of the objective
%   function and its different terms at each iteration
%
%   An example of the use of this function on a real dataset is shown in 
%   the file demo.m.
%
%   The algorithm is presented in detail in:
%
%   L. Drumetz, M. A. Veganzones, S. Henrot, R. Phlypo, J. Chanussot and 
%   C. Jutten, "Blind Hyperspectral Unmixing Using an Extended Linear
%   Mixing Model to Address Spectral Variability," in IEEE Transactions on 
%   Image Processing, vol. 25, no. 8, pp. 3890-3905, Aug. 2016.
%
%   Last Revision: 17-November-2016.
%   Revision: 1.0
%


minnargin = 7; % check number of inputs
maxnargin = 16;
narginchk(minnargin,maxnargin);

P = size(A_init,1); % number of endmembers

% find out if regularization parameters are scalar or vectors


% abundances

if length(lambda_a) == 1
    scalar_lambda_a = 1;
elseif length(lambda_a) == P
    scalar_lambda_a = 0;
    if size(lambda_a,1) == 1
        lambda_a = lambda_a';
    end
else
    error('ELMM_ADMM:wrong_size_regularization','lambda_a must be a scalar or a P-dimensional vector')
end

% scaling factors

if length(lambda_psi) == 1
    scalar_lambda_psi = 1;
elseif length(lambda_psi) == P
    scalar_lambda_psi = 0;
    if size(lambda_psi,1) == 1
        lambda_psi = lambda_psi';
    end
else
    error('ELMM_ADMM:wrong_size_regularization','lambda_psi must be a scalar or a P-dimensional vector')
end

% initialize optional parameters

maxiter_anls = 100;
maxiter_admm = 100;
epsilon_s = 10^(-3);
epsilon_a = 10^(-3);
epsilon_psi = 10^(-3);
epsilon_admm_abs = 10^(-2);
epsilon_admm_rel = 10^(-2);
norm_sr = '1,1';
verbose = true;

% change optional parameters if provided

for i = 1:length(varargin)
    switch i
        case 1
            if strcmp(varargin{i},'2,1') || strcmp(varargin{i},'1,1')
                norm_sr = varargin{i};
            else
                error('ELMM_ADMM:unknown_norm','norm string should be either 2,1 or 1,1');
            end
        case 2
            verbose = varargin{i};
        case 3
            maxiter_anls = varargin{i};
        case 4
            maxiter_admm = varargin{i};
        case 5
            epsilon_s = varargin{i};
        case 6
            epsilon_a = varargin{i};
        case 7
            epsilon_psi = varargin{i};
        case 8
            epsilon_admm_abs = varargin{i};
        case 9
            epsilon_admm_rel = varargin{i};
    end
end

%% initialize variables


[m,n,L] = size(data); % dimensions of the data cube
N = m*n; % number of pixels

data_r = reshape(data,N,L)'; % store image as a matrix

% relative variations of the variables

rs = zeros(maxiter_anls,1);
ra = zeros(maxiter_anls,1);
rpsi =  zeros(maxiter_anls,1);

% initialize variables

A = A_init; % initialize the abundance matrix
S = repmat(S0,[1,1,N]); % initialize pixel-dependent endmember matrix
psi_maps = psis_init; % initialize scaling factors

S0ptS0 = diag(S0'*S0); % precomputing

% optimization output

objective = zeros(maxiter_anls,1);
norm_fitting = zeros(maxiter_anls,1);
source_model = zeros(maxiter_anls,1);

if scalar_lambda_a
    TV_a = zeros(maxiter_anls,1);
else
    TV_a = zeros(maxiter_anls,P);
end

if scalar_lambda_psi
    smooth_psi = zeros(maxiter_anls,1);
else
    smooth_psi = zeros(maxiter_anls,P);
end

%% ADMM parameters

% build handlers and necessary stuff

% forward first order horizontal difference operator
FDh = zeros(m,n);
FDh(end,1) = -1;
FDh(end,end) = 1;
FDh = fft2(FDh);
FDhC = conj(FDh);

% forward first order vertical difference operator
FDv = zeros(m,n);
FDv(1,end) = -1;
FDv(end,end) = 1;
FDv = fft2(FDv);
FDvC = conj(FDv);

% barrier parameter of ADMM and related stuff

rho = zeros(maxiter_admm,1);
rho(1) = 10;

tau_incr = 2;
tau_decr = 2;
nu = 10;

% auxiliary functions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define a circular convolution (the same for all bands) accepting a
% matrix  and returning a matrix

ConvC = @(X,FK)  reshape(real(ifft2(fft2(reshape(X', m,n,P)).*repmat(FK,[1,1,P]))), m*n,P)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% convert matrix to image
conv2im  = @(A)  reshape(A',m,n,P);
% convert image to matrix
conv2mat = @(A)  reshape(A,m*n,P)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% optimization

for i = 1:maxiter_anls
    
    % store previous values of the variables
    
    S_old = S;
    psi_maps_old = psi_maps;
    A_old_anls = A;
    
    %% S_update
    
    if verbose
        fprintf('updating S...\n')
    end
        
%         parfor k =1:N
        for k =1:N
            S(:,:,k) = (data_r(:,k)*A(:,k)' + lambda_s*S0*diag(psi_maps(:,k))) ...   %
                /(A(:,k)*A(:,k)' + lambda_s*eye(P));
            S(:,:,k) = max(1e-6,S(:,:,k));
        end
        
    if verbose
        fprintf('Done!\n')
    end
    
    %% A_update
    
    if verbose
        fprintf('updating A...\n')
    end
    
    if any(lambda_a)  % ADMM (with spatial regularization)
        
        % initialize split variables
        
    
        
        v1 = A;
        v1_im = conv2im(v1); % v2 needs to be an image for its update
        v2 = ConvC(A,FDh);
        v3 = ConvC(A,FDv);
        v4 = A;
        
        % initialize Lagrange multipliers
        
        d1 = zeros(P,N);
        d2 = zeros(size(v2));
        d3 = zeros(size(v3));
        d4 = zeros(size(psi_maps));
        
        mu = zeros(1,N); % this vector is for the ASC
        
        % initialize primal and dual variables
        
        primal = zeros(maxiter_admm,1);
        dual = zeros(maxiter_admm,1);
        
        % precomputing
        
        Hvv1 = ConvC(v1,FDv);
        Hhv1 = ConvC(v1,FDh);
        
        
        for iter = 1:maxiter_admm
            
            % store values from the previous iterations
            
            A_old = A;
            
            v1_old = v1;       
            p_res2_old = v2 - Hhv1;
            p_res3_old = v3 - Hvv1;
            v4_old = v4;
            d1_old = d1;
            d4_old = d4;
            
            % min w.r.t. A and mu
            
%             parfor k = 1:N
            for k = 1:N
                ALPHA = (S(:,:,k)'*S(:,:,k)+2*rho(iter)*eye(P));
                ALPHA_INVERTED = inv(ALPHA);
                BETA = ones(P,1);
                s = sum(ALPHA_INVERTED(:));
                SEC_MEMBER = [ ((S(:,:,k)'*(data_r(:,k))) + rho(iter)*(v1(:,k) + d1(:,k) + v4(:,k)+ d4(:,k))); 1];
                OMEGA_INV = [ALPHA_INVERTED*(eye(P)-1/s*ones(P,P)*ALPHA_INVERTED), 1/s * ALPHA_INVERTED*BETA ; 1/s * BETA'*ALPHA_INVERTED, -1/s];
                X = OMEGA_INV * SEC_MEMBER;
                
                A(:,k) = X(1:end-1);
                mu(k) = X(end);
            end
            

            % min w.r.t. v1
            
            % reshape necessary variables into images
            
            A_im = conv2im(A);
            d1_im = conv2im(d1);
            d2_im = conv2im(d2);
            d3_im = conv2im(d3);
            v2_im = conv2im(v2);
            v3_im = conv2im(v3);
            
            
            % update in the Fourier domain
            
%             parfor p = 1:P
            for p = 1:P
                second_term_in_spectral_domain = fft2(squeeze(A_im(:,:,p)) - squeeze(d1_im(:,:,p))) + ...
                    fft2(squeeze((v2_im(:,:,p)) + squeeze(d2_im(:,:,p)))).*FDhC + fft2(squeeze((v3_im(:,:,p)) + squeeze(d3_im(:,:,p)))).*FDvC;
                v1_im(:,:,p) = real(ifft2((second_term_in_spectral_domain)./(ones(m,n) + abs(FDh).^2+ abs(FDv).^2))); % this is the convolution performed in the fourier domain, done band by band.
            end
            
            % convert back necessary variables into matrices
            
            v1 = conv2mat(v1_im);
            Hvv1 = ConvC(v1,FDv);
            Hhv1 = ConvC(v1,FDh);
            
            
            % min w.r.t. v2 and v3;
            
            if scalar_lambda_a
                
                if strcmp(norm_sr,'2,1')
                    v2 = vector_soft_col(-(d2 - Hhv1),lambda_a/rho(iter)); % l21 norm
                    v3 = vector_soft_col(-(d3 - Hvv1),lambda_a/rho(iter));
                elseif strcmp(norm_sr,'1,1')
                    v2 = soft(-(d2 - Hhv1),lambda_a/rho(iter)); % l11 norm
                    v3 = soft(-(d3 - Hvv1),lambda_a/rho(iter));
                end
                
            else
                if strcmp(norm_sr,'2,1')
                    for p =1:P
                        v2(p,:) = vector_soft_col(-(d2(p,:) - Hhv1(p,:)),lambda_a(p)/rho(iter)); % l21 norm
                        v3(p,:)= vector_soft_col(-(d3(p,:) - Hvv1(p,:)),lambda_a(p)/rho(iter));
                    end
                elseif strcmp(norm_sr,'1,1')
                    for p = 1:P
                        v2(p,:) = soft(-(d2(p,:) - Hhv1(p,:)),lambda_a(p)/rho(iter)); % l11 norm
                        v3(p,:) = soft(-(d3(p,:) - Hvv1(p,:)),lambda_a(p)/rho(iter));
                    end
                end
            end
            
            % min wrt. v4
            
            v4 = max(A - d4,zeros(size(A)));
            
            
            % dual update
            
            % compute necessary variables for the residuals, and update the Lagrange multipliers

            
            p_res1 = v1 - A;
            p_res2 = v2 - Hhv1;
            p_res3 = v3 - Hvv1;
            p_res4 = v4 - A;
            
            d1 = d1 + p_res1;
            d2 = d2 + p_res2;
            d3 = d3 + p_res3;
            d4 = d4 + p_res4;
            
            % primal and dual residuals
            
            primal(iter) = sqrt( norm(p_res1,'fro')^2 + norm(p_res2,'fro')^2 + norm(p_res3,'fro')^2 + norm(p_res4,'fro')^2);
            dual(iter) = rho(iter) * sqrt( norm(v1_old-v1,'fro')^2 + norm(v4_old-v4,'fro')^2);
            
            % compute termination values
            
            epsilon_primal = sqrt(4*P*N) *  epsilon_admm_abs + epsilon_admm_rel*max(sqrt(2*norm(A,'fro')^2),...
                sqrt(norm(v1_old,'fro')^2 + norm(p_res2_old,'fro')^2 + norm(p_res3_old,'fro')^2 + norm(v4_old,'fro')^2));
            epsilon_dual = sqrt(P*N)*epsilon_admm_abs + rho(iter) * epsilon_admm_rel * sqrt(norm(d1_old,'fro') + norm(d4_old,'fro')^2);
            
            rel_A = abs(norm(A,'fro')-norm(A_old,'fro'))/norm(A_old,'fro');
            
            % display of ADMM results
            
            if verbose
                fprintf('iter %d of %d, rel_A = %f, primal = %f, eps_p = %d dual = %f, eps_d = %f, rho = %f \n',iter,maxiter_admm,rel_A,primal(iter),epsilon_primal,dual(iter),epsilon_dual,rho(iter));
            end
            
            if iter > 1 && ((primal(iter) < epsilon_primal && dual(iter) < epsilon_dual));
                break
            end
            
            % rho update
            
            if iter < maxiter_admm
                if norm(primal(iter)) > nu * norm(dual(iter))
                    rho(iter+1) = tau_incr*rho(iter);
                    A = A/tau_incr;
                elseif norm(dual(iter)) > nu * norm(primal(iter))
                    rho(iter+1) = rho(iter)/tau_decr;
                    A = tau_decr*A;
                else
                    rho(iter+1) = rho(iter);
                end
            end
            
        end
        
    else % FCLSU (without spatial regularization)
        
%         parfor k =1:N
        for k =1:N
            A(:,k) =  FCLSU(data_r(:,k),S(:,:,k));
        end
        
    end
    
    if verbose
        fprintf('Done!\n')
    end
    
    %% psi_update
    
    if verbose
        fprintf('updating psi...\n')
    end
    
    if any(lambda_psi) % with spatial regularization
        
        if scalar_lambda_psi %% update with scalar or vector lambda_s (done in the Fourier domain)
%                 parfor p = 1:P
                for p = 1:P
                    numerator = reshape(lambda_s*squeeze(S(:,p,:))'*S0(:,p),m,n);
                    psi_maps_im =  real(ifft2(fft2(numerator)./((lambda_psi*(abs(FDh).^2+ abs(FDv).^2) + lambda_s * S0ptS0(p)))));
                    psi_maps(p,:) = psi_maps_im(:);
                end
        else
%                 parfor p = 1:P
                for p = 1:P
                    numerator = reshape(lambda_s*squeeze(S(:,p,:))'*S0(:,p),m,n);
                    psi_maps_im =  real(ifft2(fft2(numerator)./((lambda_psi(p)*(abs(FDh).^2+ abs(FDv).^2) + lambda_s * S0ptS0(p)))));
                    psi_maps(p,:) = psi_maps_im(:);
                end

        end
    else % without spatial regularization
        for p=1:P
            psi_maps_temp = zeros(N,1);
%             parfor k = 1:N
            for k = 1:N
                psi_maps_temp(k) =  (S0(:,p)'*S(:,p,k))/S0ptS0(p);
            end
            psi_maps(p,:) = psi_maps_temp;
        end
    end
    
    if verbose
        fprintf('Done!\n')
    end
    
    % residuals of the ANLS loops
    
    rs_vect = zeros(N,1);
    
%     parfor k=1:N
    for k=1:N
        rs_vect(k) = norm(squeeze(S(:,:,k))-squeeze(S_old(:,:,k)),'fro')/norm(squeeze(S_old(:,:,k)),'fro');
    end
    
    rs(i) = mean(rs_vect);
    ra(i) = norm(A(:)-A_old_anls(:),2)/norm(A_old_anls(:),2);
    rpsi(i) = norm(psi_maps-psi_maps_old,'fro')/(norm(psi_maps_old,'fro'));
    
    % compute objective function value
    
    SkAk = zeros(L,N);
    
%     parfor k = 1:N
    for k = 1:N
        SkAk(:,k) = squeeze(S(:,:,k)*A(:,k));
        S0_psi(:,:,k) = S0 *diag(psi_maps(:,k));
    end
    
    norm_fitting(i) = 1/2*norm(data_r(:)-SkAk(:))^2;
    
    source_model(i) = 1/2 * norm(S(:)-S0_psi(:))^2;

    
    if any(lambda_psi) && any(lambda_a) % different objective functions depending on the chosen regularizations
        
        if scalar_lambda_psi
            smooth_psi(i) = 1/2 * (sum(sum((ConvC(psi_maps,FDh).^2))) +  sum(sum((ConvC(psi_maps,FDv).^2))));
        else
            
            CvCpsih = ConvC(psi_maps,FDh);
            CvCpsiv = ConvC(psi_maps,FDv);
            
%             parfor p =1:P
            for p =1:P
                smooth_psi(i,p) = 1/2 * (sum(sum((CvCpsih(p,:).^2))) +  sum(sum((CvCpsiv(p,:).^2))));
            end
        end
        
        if scalar_lambda_a
            if strcmp(norm_sr,'2,1')
                TV_a(i) = sum(sum(sqrt(ConvC(A,FDh).^2 + ConvC(A,FDv).^2)));
            elseif  strcmp(norm_sr,'1,1')
                TV_a(i) = sum(sum(abs(ConvC(A,FDh)) + abs(ConvC(A,FDv))));
            end
        else
            
            CvCAh = ConvC(A,FDh);
            CvCAv = ConvC(A,FDv);
            
            if strcmp(norm_sr,'2,1')
%                 parfor p =1:P
                for p =1:P
                    TV_a(i,p) = sum(sum(sqrt(CvCAh(p,:).^2 + CvCAv(p,:).^2)));
                end
            elseif  strcmp(norm_sr,'1,1')
%                 parfor p =1:P
                for p =1:P
                    TV_a(i,p) = sum(sum(abs(CvCAh(p,:)) + abs(CvCAv(p,:))));
                end
            end
        end
        
        objective(i) = norm_fitting(i) + lambda_s * source_model(i) + lambda_a' * TV_a(i,:)' + lambda_psi' * smooth_psi(i,:)';
        
        % display
        
        if verbose
            fprintf('iteration %d of %d, rs = %f, ra = %f, rpsi= %f\n',i,maxiter_anls,rs(i),ra(i),rpsi(i));
            fprintf('objective = %f, norm_fitting = %f, source_model = %f, smooth_psi = %f, TV_A = %f\n',...
                objective(i), norm_fitting(i), lambda_s * source_model(i), lambda_psi' * smooth_psi(i,:)',  lambda_a' * TV_a(i,:)');
        end
        
    elseif ~any(lambda_psi) && any(lambda_a)
        
        
        if scalar_lambda_a
            if strcmp(norm_sr,'2,1')
                TV_a(i) = sum(sum(sqrt(ConvC(A,FDh).^2 + ConvC(A,FDv).^2)));
            elseif  strcmp(norm_sr,'1,1')
                TV_a(i) = sum(sum(abs(ConvC(A,FDh)) + abs(ConvC(A,FDv))));
            end
        else
            
            CvCAh = ConvC(A,FDh);
            CvCAv = ConvC(A,FDv);
            
            if strcmp(norm_sr,'2,1')
%                 parfor p =1:P
                for p =1:P
                    TV_a(i,p) = sum(sum(sqrt(CvCAh(p,:).^2 + CvCAv(p,:).^2)));
                end
            elseif  strcmp(norm_sr,'1,1')
%                 parfor p =1:P
                for p =1:P
                    TV_a(i,p) = sum(sum(abs(CvCAh(p,:)) + abs(CvCAv(p,:))));
                end
            end
        end
        
        objective(i) = norm_fitting(i) + lambda_s * source_model(i) + lambda_a' * TV_a(i,:)';
        % display
        
        if verbose
            fprintf('iteration %d of %d, rs = %f, ra = %f, rpsi= %f\n',i,maxiter_anls,rs(i),ra(i),rpsi(i));
            fprintf('objective = %f, norm_fitting = %f, source_model = %f, TV_A = %f\n',objective(i), norm_fitting(i), lambda_s * source_model(i), lambda_a' * TV_a(i,:)');
        end
        
    elseif any(lambda_psi) && ~any(lambda_a);
        
        if scalar_lambda_psi
            smooth_psi(i) = 1/2 * (sum(sum((ConvC(psi_maps,FDh).^2))) +  sum(sum((ConvC(psi_maps,FDv).^2))));
        else
            
            CvCpsih = ConvC(psi_maps,FDh);
            CvCpsiv = ConvC(psi_maps,FDv);
            
%             parfor p =1:P
            for p =1:P
                smooth_psi(i,p) = 1/2 * (sum(sum((CvCpsih(p,:).^2))) +  sum(sum((CvCpsiv(p,:).^2))));
            end
        end
        
        objective(i) = norm_fitting(i) +  lambda_s * source_model(i) + lambda_psi' * smooth_psi(i,:)';
        
        % display
        
        if verbose
            fprintf('iteration %d of %d, rs = %f, ra = %f, rpsi= %f\n',i,maxiter_anls,rs(i),ra(i),rpsi(i));
            fprintf('objective = %f, norm_fitting = %f, source_model = %f, smooth_psi = %f\n',objective(i), norm_fitting(i), lambda_s * source_model(i), lambda_psi' * smooth_psi(i,:)');
        end
        
    else
        objective(i) = norm_fitting(i) + lambda_s * source_model(i);
        
        % display
        
        if verbose
            fprintf('iteration %d of %d, rs = %f, ra = %f, rpsi= %f\n',i,maxiter_anls,rs(i),ra(i),rpsi(i));
            fprintf('objective = %f, norm_fitting = %f, source_model = %f\n',objective(i), norm_fitting(i), lambda_s * source_model(i));
        end
        
    end
    
    % termination test
    
    if ((rs(i) < epsilon_s) && (ra(i) < epsilon_a)) && (rpsi(i)< epsilon_psi)
        break;
    end
    
end

% recover optimization results

optim_struct.obj = objective;
optim_struct.fit = norm_fitting;
optim_struct.regul = source_model;
optim_struct.TVa = TV_a;
optim_struct.smoothpsi = smooth_psi;

end

