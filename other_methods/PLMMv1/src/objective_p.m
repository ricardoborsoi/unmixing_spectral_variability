function [ f ] = objective_p(Y,M,A,dM,H,W,alpha,gamma,varargin)
% Computation of the objective function value.
%%
% Code : Pierre-Antoine Thouvenin, February 16th 2015.
%%
%-------------------------------------------------------------------------%
%   Y   =  M    A   + [dM1*a1|...|dMn*an] +   B     (LMM)
% (L|N)  (L|K)(K|N)           (dMA)         (L|N)
%
% 0.5*|Y - MA - dMA|² + 0.5*alpha\phi(A) + 0.5*beta\psi(M) + 0.5*gamma|dM|²
%  
% s.t.            A >=  0         with N : number of pixels;    
%                 M >=  0              K : number of endmembers;
% forall n  M + dMn >=  0              L : number of spectral bands.
%             A'   1 =  1
%           (N|K)(K|1)(N|1)
%-------------------------------------------------------------------------%
% Input: 
% > Y   hypespectral image (L|N);
% > M   initial endmembers (L|K);
% > A   initial abundances (L|N);
% > dM  initial variability terms [cell(1|N)];
% > H   image height;
% > W   image width.
% > gamma  penalty parameter (variability);
% > alpha  penalty parameter (spatial smoothness);
%
% Optional parameters
% > delta  lower bound of the Lipschitz's coefficients;
% > beta   penalty parameter ('NONE','DISTANCE','MUTUAL DISTANCE','VOLUME')
%
% >> DISTANCE
% --> M0   reference signatures
%
% >> VOLUME
% --> Y_proj  projected data(PCA)(K-1|N)
% --> Y_bar   data mean (mean(Y,2))(L|1)
% --> U       inverse projector (PCA variables -> M)(L|K-1);
% --> V       projector (M -> PCA variables)(K-1|L);
%
% Output:
% < f   objective function value.
%-------------------------------------------------------------------------%
%%
%--------------------------------------------------------------
% Local variables
%--------------------------------------------------------------
% Size
[K,N] = size(A);
L = size(M,1);

%--------------------------------------------------------------
% Optional parameters
%--------------------------------------------------------------
if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'PENALTY'
                type = varargin{i+1}{1};
                beta = varargin{i+1}{2};
                switch upper(type)
                    case 'DISTANCE'
                        M0 = varargin{i+1}{3};    
                    case 'VOLUME'
                        Y_proj = varargin{i+1}{3};
                        Y_bar = varargin{i+1}{4};
                        U = varargin{i+1}{5};
                        V = varargin{i+1}{6};
                end
            otherwise
                type = 'NONE';
        end
    end
end

%--------------------------------------------------------------
% Abundance penalization term
%--------------------------------------------------------------
if alpha == 0
        termA = 0;
else
        %-Spatial smoothness constraint
        UpA = zeros(K,N);
        for k = 1:N-W
            UpA(:,W+k) = A(:,W+k) - A(:,k);
        end
        DownA  = [-UpA(:,W+1:end),zeros(K,W)];
        LeftA  = zeros(K,N);
        RightA = zeros(K,N);
        for k = 0:H-1
            LeftA(:,1+k*W:(k+1)*W)  = [zeros(K,1), diff(A(:,1+k*W:(k+1)*W),1,2)];
            RightA(:,1+k*W:(k+1)*W) = [-LeftA(:,2+k*W:(k+1)*W), zeros(K,1)];
        end
        AH = [LeftA,RightA,UpA,DownA];
        termA = 0.5*alpha*(norm(AH,'fro')^2);
end

%--------------------------------------------------------------
% Endmember penalization term
%--------------------------------------------------------------
switch upper(type)
    case 'NONE'
        termM = 0;
    case 'DISTANCE'
        termM = 0.5*beta*(norm(M-M0,'fro')^2); 
    case 'MUTUAL DISTANCE'
        %-Distance between endmembers
        R = zeros(K,K*K);
        for l = 0:K-1
            el = zeros(K,1);
            el(l+1) = 1;
            R(:,1+l*K:(l+1)*K) = -eye(K) + el*ones(1,K);
        end
        termM = 0.5*beta*(norm(M*R,'fro')^2);
    case 'VOLUME'
        %-Volume
        T = V*(M - Y_bar(:,ones(K,1))); 
        Volume = det([T;ones(1,K)])/factorial(K-1);
        termM = 0.5*beta*(Volume^2);
    otherwise
        error(['Unrecognized penalty type: ''' type '''']);
end

%--------------------------------------------------------------
% Objective function
%--------------------------------------------------------------
Norm_dM = zeros(1,N);
dMA     = zeros(L,N);
for n = 1:N
    Norm_dM(n) = norm(dM{n},'fro')^2; 
    dMA(:,n)   = dM{n}*A(:,n);
end

% Objective
f = 0.5*(norm(Y-M*A-dMA,'fro')^2) + termA + termM + ...
    0.5*gamma*sum(Norm_dM);

end
