function [A,r_err,bestIdx] = MESMA(y, Lib )
% =========================================================================
% Runs Multiple Endmember Spectral Mixture Analysis (MESMA), unmixing
% each pixel of y with the combination of endmembers in the library
% Lib that results in the smallest reconstruction error.
% 
% INPUT
%   y: observed image, L-by-N (bands*pixels)
%   Lib: spectral library, P-by-1 cell array of L-by-N_lib(p) endmembers
% 
% OUTPUT
%   A: estimated abundances, P-by-N
%   r_err: reconstruction error for each pixel (N-by-1)
%   bestIdx: index of the signatures from the library which were used to
%            unmix each pixel (N-by-P)
% 
% Ricardo Borsoi, 2020
% =========================================================================


P = length(Lib);
L = size(Lib{1},1);
N = size(y,2);

% get library sizes
N_lib = zeros(P,1);
for i=1:P
    N_lib(i) = size(Lib{i},2);
end


% initializations
N_tot = 1;
for i=1:P
    N_tot = N_tot * N_lib(i); % total number of iterations
end


% --------------------------------
% indexing

linidxs = 1:N_tot; % linear indices
% convert to vector indices
sizess = [];
for i=1:P
    sizess = [sizess N_lib(i)];
end
allidxvec = ind2sub2(sizess,linidxs); % allidxvec = [id_N(1) id_N(2) ... id_N(P)]


rec_err = inf*ones(N_tot,N);
% M = zeros(L,P);
% parfor ii=1:N_tot
for ii=1:N_tot
    % disp(ii)
    MM = zeros(L,P);
    for p=1:P
        MM(:,p) = Lib{p}(:,allidxvec(ii,p));
    end
    
    a = qpas2(MM,y);
    rec_err(ii,:) = sum((y-MM*a).^2,1);

end




% Find combinations with the least RE
[minval,minidx] = min(rec_err);

bestIdx = zeros(N,P);
A = zeros(P,N);
r_err = zeros(1,N);

% M = zeros(0);
% parfor i=1:N
for i=1:N
    
    bestIdx(i,:) = allidxvec(minidx(i),:)';
    M = zeros(L,P);
    for p=1:P
        M(:,p) = Lib{p}(:,allidxvec(minidx(i),p));
    end
    A(:,i) = qpas2(M,y(:,i));
    r_err(i) = norm(y(:,i)-M*A(:,i),2)^2;
end

end


function [a] = qpas2(M,y)

option_alg = 5;

switch option_alg
    case 1
        P = size(M,2);
        N = size(y,2);
        H = M'*M;
        Aeq = ones(1,P);
        beq = 1;
        LB = zeros(P,1);
        UB = inf*ones(P,1);
        options = optimoptions('quadprog', 'Algorithm','interior-point-convex','TolX',1e-6,'Display','off'); 
        a = zeros(P,N);
        for i=1:N
            f = -M'*y(:,i);
            a(:,i) = quadprog(H,f,[],[],Aeq,beq,LB,UB,[],options);
        end
    
    case 2
        % other option
        P = size(M,2);
        N = size(y,2);
        C = M;
        Aeq = ones(1,P);
        beq = 1;
        lb = zeros(P,1);
        ub = inf*ones(P,1);
        options = optimoptions('TolX',1e-6,'Display','off'); 
        a = zeros(P,N);
        for i=1:N
            d = y(:,i);
            a(:,i) = lsqlin(C,d,[],[],Aeq,beq,lb,ub,[],options);
        end

    case 3
        L = size(y,1);
        P = size(M,2); 
        N = size(y,2); 
        ddelta = 2e3; 2e3; 1e4; 1e4; 1e5;
        C = [M;  ddelta * ones(1,P)];
        d = [y;  ddelta * ones(1,N)];
%         options = optimoptions('TolX',1e-6,'Display','off'); 
        options = optimset('Display','off','TolX',1e-16);
        % options = optimset('Display','off'); 
        a = lsqnonnegvect(C,d, options);
        % a = lsqnonnegvect(C,d); 

        
    case 4
        P = size(M,2);
        a = qpas(M'*M, -M'*y,[],[],ones(1,P),1,zeros(P,1),1e15*ones(P,1));

    case 5
        % sunsal's code for FCLS
        a = sunsal(M,y,'POSITIVITY','yes','ADDONE','yes'); 
        % a = sunsal(M,y,'POSITIVITY','yes','ADDONE','yes','TOL',1e-4);
        
        
        
    otherwise
        error('Unknown NNLS/FLCS algorithm selected for MESMA!')
end

end



