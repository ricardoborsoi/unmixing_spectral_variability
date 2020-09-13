function gam_NL = sunsal_gam_c2(Y_NL2,M_NL2,EVe,EVa,s,R,Posi,C_EV);


% data set size
[L,N]    = size(Y_NL2);
D        = R*(R+1)/2;

AL_iters = 100;%500;
tol      = 1e-4;
tol1     = sqrt(D*N)*tol;
tol2     = sqrt(D*N)*tol;


MY       = (M_NL2'* (Y_NL2)).*(repmat(C_EV(:).^2,1,D))'; %DxN
mu       = 100;%0.01;%


small    = 10^(-6);
EVa2     = EVa*(C_EV(:).^4)' +ones(D,1)*(1./s(:)'); % DxN
IF0      = (repmat(EVe,N,1)./kron(EVa2'+small, ones(D,1)))*EVe'; % NDxD
MY0      = kron((MY)' , ones(D,1)); % DN x D
x        = reshape(sum(IF0.*MY0,2), D,N);

if(Posi == 0)
    gam_NL = x;
    return
end

IF   = (repmat(EVe,N,1)./kron(EVa2'+mu, ones(D,1)))*EVe'; % NDxD
x    = max(0,x);
z    = x;
% scaled Lagrange Multipliers
d     = 0*z;
z0    = x;
i     = 1;
res_p = inf;
res_d = inf;
mu_changed = 0;
muu   = [mu];

while (i <= AL_iters) && ((abs (res_p) > tol1) || (abs (res_d) > tol2))
    %     i
    % save z to be used later
    if mod(i,10) == 1
        z0 = z;
    end
    z    = max(0,x);
    
    %%%%%%%%%%%%
    %%%%%%%%%%%%
    MY2  = kron((MY + mu*(z +d))' , ones(D,1)); % DN x D
    x    = reshape(sum(IF.*MY2,2), D,N);% D x N
    
    % Lagrange multipliers update
    d = d -(x-z);
    
    % update mu so to keep primal and dual residuals whithin a factor of 10
    if mod(i,10) == 1
        %         i
        % primal residue
        res_p = norm(x-z,'fro');
        % dual residue
        res_d = mu*norm(z-z0,'fro');
        % update mu
        if res_p > 10*res_d
            mu = mu*2;%muu = [muu mu];
            d = d/2;
            mu_changed = 1;
        elseif res_d > 10*res_p
            mu = mu/2;%muu = [muu mu];
            d = d*2;
            mu_changed = 1;
        end
        muu = [muu mu];
        
        if  mu_changed
            % update IF and IF1
            
            IF  = (repmat(EVe,N,1)./kron(EVa2'+mu, ones(D,1)))*EVe'; % NDxD
            
        end
    end
    i=i+1;
end
gam_NL =  (x);



