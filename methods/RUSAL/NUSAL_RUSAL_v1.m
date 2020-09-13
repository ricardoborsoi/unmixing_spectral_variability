
function [zt,stop_flag,i] = NUSAL_RUSAL_v1(ymat,M,tau1,tau2,AL_iters,tol,Model,Ab_init)

%% NUSAL_RUSAL_v1(ymat,M,tau1,tau2,AL_iters,tol,Model)
%  Paper : A. Halimi, J. M. Bioucas-Dias, N. Dobigeon, G. S. Buller
%          and S. McLaughlin, "Fast Hyperspectral Unmixing in Presence of
%          Nonlinearity or Mismodelling Effects", IEEE Trans. Comput. Imaging., 2017.
%  Code  :  NUSAL + RUSAL
%  Version (Nov., 2016) by Abderrahim Halimi (a.halimi@hw.ac.uk)
%  For any comments contact the author
%
%% --------------- Inputs --------------------------------------------
%%% ymat      : pixels of size (row x col x L)
%%% MPlus     : endmembers of size (L x R)
%%% tau1      : regularization parameter for the l_1 norm
%%% tau2      : regularization parameter for the l_{21} norm
%%% AL_iters  : maximum number of iterations  (> 10)
%%% tol       : stopping value for the residuals
%%% Model     : choice between NL or ME models
%%%           : ='NUSAL2'  for the NL model of order 2
%%%           : ='NUSAL3'  for the NL model of order 3
%%%           : ='RUSAL'  for the ME model
%%%
%% --------------- Outputs --------------------------------------------
%%% zt         : Chain of the estimated parameters   (R+D x N x Iter/10)
%%% stop_flag  : Function ended because
%%%              stop_flag=0: the maximum number of iterations was attained
%%%              stop_flag=1: the convergence criterion was satisfied
%%% i          : Number of iterations until convergence
%% --------------------------------------------------------------------
%% -------------------------------------------------------------------------
%
% Copyright (Nov., 2016):        Abderrahim Halimi (a.halimi@hw.ac.uk)
%
% CDA-NL is distributed under the terms of
% the GNU General Public License 3.0.
%
% Permission to use, copy, modify, and distribute this software for
% any purpose without fee is hereby granted, provided that this entire
% notice is included in all copies of any software which is or includes
% a copy or modification of this software and in all copies of the
% supporting documentation for such software.
% This software is being provided "as is", without any express or
% implied warranty.  In particular, the authors do not make any
% representation or warranty of any kind concerning the merchantability
% of this software or its fitness for any particular purpose."
% ---------------------------------------------------------------------

%---------------------------------------------
%   Constants
%---------------------------------------------
[LM R]        = size(M);
[row col L]  = size(ymat);
N            = row*col;
y            = reshape(ymat,N,L)'; % L x N
F            = M'*M;

if(LM ~= L)
    error('mixing matrix M and data set y are inconsistent');
end
if (AL_iters <= 0 )
    error('AL_iters must be a positive integer');
end
if (AL_iters <= 10 )
    error('Please choose AL_iters higher than 10');
end
if (tau1 <= 0 | tau2 <= 0 )
    error('Regularization parameters tau1 and tau2 should be strictly positive');
end

%%%%%%%%%%%
if(strcmp(Model,'RUSAL'))  %%%%% Robust unmixing
    D    = 20;             % size of the DCT coefficients
    dL   = idct(eye(L)); %dct(eye(L));
    Matn = dL(:,1:D);
    DtD  = eye(D);
    InvKtK   = 1/3;
    cstePSTO = eye(R)-1/R*ones(R);
    %%%%%%%%%%%
elseif(strcmp(Model,'NUSAL2'))  %%%%% Nonlinear unmixing
    D = R*(R-1)/2+R;      % size of the nonlinear coefficients
    M_NL    = [];
    for i=1:R
        for j=i+1:R
            M_NL   = [M_NL  M(:,i).*M(:,j)];
        end
    end
    M_NL   = [sqrt(2)*M_NL M.^2]; % LxD
    Matn   = M_NL;
    D      = size(M_NL,2);
    DtD    = Matn'*Matn;
    InvKtK   = repmat([1/3*ones(R,1);1/4*ones(D,1)],1,N);%
    cstePSTO = eye(R)-1/R*ones(R);
    %%%%%%%%%%%
elseif(strcmp(Model,'NUSAL3'))  %%%%% Nonlinear unmixing
    M_NL2    = [];M_NL3    = [];
    for i=1:R
        for j=i+1:R
            M_NL2   = [M_NL2  M(:,i).*M(:,j)];
        end
    end
    M_NL2   = [sqrt(2)*M_NL2 M.^2]; % LxD
    
    for i=1:R
        for j=1:R
            if(i~=j)
                M_NL3   = [M_NL3  sqrt(3)*M(:,i).^2.*M(:,j)];
            end
        end
    end
    for i=2:R
        for j=2:i-1
            for k=1:j-1
                M_NL3   = [M_NL3  sqrt(6)*M(:,i).*M(:,j).*M(:,k)];
            end
        end
    end
    M_NL3   = [M_NL3 M.^3]; % LxD
    M_NL    = [M_NL2 M_NL3];
    Matn   = M_NL;
    D      = size(M_NL,2);
    DtD    = Matn'*Matn;
    InvKtK   = repmat([1/3*ones(R,1);1/4*ones(D,1)],1,N);%
    cstePSTO = eye(R)-1/R*ones(R);
else
    error('Model should take the values: NUSAL2, NUSAL3 or RUSAL')
end

mu       = 0.01;%1; %
Md       = [M, Matn]; % LX  D+R
idM      = Matn'*M;% DxR;
[UF,SF]  = svd([F  idM';idM DtD] ) ;
sF       = diag(SF); % L+R x 1
Igam     = UF*diag(1./(sF+mu))*UF';% L+R x L+R
yCst     = Igam* Md'*y; % L+R x N



%---------------------------------------------
%   Initialization
%---------------------------------------------

%%%%
%%%% least squares constrained (sum(x) = 1)
%%%%
SMALL = 1e-12;
B     = ones(1,R);
a     = ones(1,N);
% test if F is invertible
if rcond(F) > SMALL
    % compute the solution explicitly
    IF = inv(F);
    Ab_init = IF*M'*y-IF*B'*inv(B*IF*B')*(B*IF*M'*y-a);
end

%%%%
%%%% Auxiliar variables
%%%%
u1   = [Ab_init;10^(-4)*ones(D,N)]; % L+R x N
u2   = 10^(-4)*ones(D,N); %L x1
u3   = 10^(-4)*ones(D,N); %L x1

if(strcmp(Model,'RUSAL'))  % Sparsity only
    u4   = Ab_init;  % R xN
elseif(strcmp(Model,'NUSAL2') | strcmp(Model,'NUSAL3'))
    u4   = [Ab_init;10^(-4)*ones(D,N)];  % R xN
end

u5   = Ab_init;  % R xN
d1   = zeros(size(u1));
d2   = zeros(size(u2));
d3   = zeros(size(u3));
d4   = zeros(size(u4));
d5   = zeros(size(u5));
z    = [Ab_init;zeros(D,N)]; % L+R x N

%---------------------------------------------
%  AL iterations - main body
%---------------------------------------------
i             = 1;
stop_flag     = 0;
res_p         = 10^6;  
mu_changed    = 0;
tol1          = sqrt(N*R)*tol;
tol2          = tol1;

while (i <= AL_iters) && (stop_flag==0)
    if mod(i,10) == 1
        u10 = u1;
        u20 = u2;
        u30 = u3;
        u40 = u4;
        u50 = u5;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%     First part: linear system   %%%%%%%
    
    psi1  = u1 + d1; % L+R x N
    psi2  = u2 + d2; % D x N
    psi3  = u3 + d3; % D x N
    psi4  = u4 + d4; % R x N
    
    switch Model
        case 'RUSAL'
            %         if(strcmp(Model,'RHUSpLa'))  % Sparsity only
            psi5    = u5 + d5; % R x N
            Hppsi   = psi3;
            z       = InvKtK*(psi1 + [psi4+psi5;psi2+Hppsi]);  % L+R x N
            
        case 'NUSAL2'
            %         elseif(strcmp(Model,'NLUSpLa')  ) % Sparsity and positivity
            psi5    = u5 + d5; % R x N
            %%% coeff sparse et non dct
            Hppsi = psi3;
            z       = InvKtK.*(psi1 + psi4 + [psi5;psi2+Hppsi]);  % L+R x N
        case 'NUSAL3'
            %         elseif(strcmp(Model,'NLUSpLa')  ) % Sparsity and positivity
            psi5    = u5 + d5; % R x N
            %%% coeff sparse et non dct
            Hppsi = psi3;
            z       = InvKtK.*(psi1 + psi4 + [psi5;psi2+Hppsi]);  % L+R x N
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%  Second part: MPO  %%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%     G1      %%%%%%%%%%%%%%%%
    v1            = z - d1;  % D+R x N
    u1            = yCst + mu*Igam*v1 ; % D+R x N
    d1            = -v1 + u1;
    
    
    %%%%%%%%%%%%%%%%%     G2      %%%%%%%%%%%%%%%%
    ztemp   = z(R+1:D+R,:) ;
    ztemp2  = ztemp;
    v2      = ztemp2  -  d2; % D x N
    u2      = soft(v2,tau1/mu); % L x N  check
    d2      = -v2 + u2;
    
    
    %%%%%%%%%%%%%%%%%     G3      %%%%%%%%%%%%%%%%
    ztemp3  = ztemp;
    v3      = ztemp   -  d3; % D x N
    % L12 constraint
    u3      = vector_soft_row(v3',tau2/mu); % L x N
    u3      = u3';
    d3      = -v3 + u3;
    
    
    %%%%%%%%%%%%%%%%%     G4      %%%%%%%%%%%%%%%%
    if(strcmp(Model,'RUSAL'))
        v4    = z(1:R,:) -  d4;  % R x N
        u4    = max(v4,0);       % R x N
        d4    = -v4 + u4;
    elseif(strcmp(Model,'NUSAL2') |  strcmp(Model,'NUSAL3'))
        v4    = z  -  d4;        % R x N
        u4    = max(v4,0);       % R x N
        d4    = -v4 + u4;
    end
    
    %%%%%%%%%%%%%%%%%     G5      %%%%%%%%%%%%%%%%
    v5    = z(1:R,:) -  d5;  % R x N
    u5    = cstePSTO*v5+1/R*ones(R,N); % R x N
    d5    = -v5 + u5;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%     Convergence   %%%%%%%%%%%%%%%%
    
    if mod(i,10) == 1
        
        if(strcmp(Model,'RUSAL'))  % Sparsity only
            res_p        = sqrt(norm(z-u1,'fro')^2 + norm(ztemp2-u2,'fro')^2+ norm(ztemp3-u3,'fro')^2 + norm(z(1:R,:)-u4,'fro')^2+ norm(z(1:R,:)-u5,'fro')^2);
            res_d        = mu*norm((u1-u10) + [(u4-u40)+(u5-u50);(u2-u20)+(u3-u30)],'fro');
        elseif(strcmp(Model,'NUSAL2') | strcmp(Model,'NUSAL3') ) % Sparsity and positivity
            res_p        = sqrt(norm(z-u1,'fro')^2 + norm(ztemp2-u2,'fro')^2+ norm(ztemp3-u3,'fro')^2 + norm(z-u4,'fro')^2+ norm(z(1:R,:)-u5,'fro')^2);
            res_d        = mu*norm((u1-u10) +(u4-u40)+ [(u5-u50);(u2-u20)+(u3-u30)],'fro');
        end
        
        % update mu
        if res_p > 10*res_d
            mu = mu*2;
            d1 = d1/2;
            d2 = d2/2;
            d3 = d3/2;
            d4 = d4/2;
            d5 = d5/2;
            
            mu_changed = 1;
        elseif res_d > 10*res_p
            mu = mu/2;
            d1 = d1*2;
            d2 = d2*2;
            d3 = d3*2;
            d4 = d4*2;
            d5 = d5*2;
            mu_changed = 1;
        end
        if  mu_changed
            Igam     = UF*diag(1./(sF+mu))*UF';% L+R x L+R
            yCst     = Igam* Md'*y; % L+R x N
            mu_changed = 0;
        end  
        stop_flag      = ((abs (res_p) < tol1) || (abs (res_d) < tol2));
        zt(:,:,(i-1)/10+1) = z;
    end
    
    i=i+1;
    
end






