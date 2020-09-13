function [alphaPlus_out, rho,sig] = sample_alphaPlus(alphaPlus,y,MUr,sigma2r,m_compt,Tx,sig)

[L R] = size(MUr);
rho = zeros(1,R); % rho = 0 if proposed sample rejected, rho = 1 else

% selecting the first abundance coefficient to be sampled
 temp = randperm(R);
 r = temp(1);
 comp_r = setdiff(temp,r);

 for k = comp_r
    ind_k = find(comp_r == k);
    Sk = sigma2r(comp_r); Sk(ind_k) = [];
    Muk = MUr(:,comp_r); Muk(:,ind_k) = [];
    alphak = alphaPlus(comp_r); alphak(ind_k) = [];
    if isempty(alphak) & isempty(Sk) & isempty(Muk) % case R = 2
        alphak = 0;
        Sk = 0;
        Muk = 0;
    end    

    % random walk
    
    % determining the Gaussian variance in function of the acceptance rate
    % every 100 iterations
    
       if Tx(k) > 0.4 && rem(m_compt-1,100) == 0 
             sig(k) = sig(k)*5;
             
        elseif Tx(k) < 0.3 && rem(m_compt-1,100) == 0 
             sig(k) = sig(k)/5;
       end    
            %proposed sample
            alpha = alphaPlus(k) + sqrt(sig(k)).*randn ;
             
     if (alpha > 0 & alpha < (1 - sum(alphak))) %constraints

         alpha_star = alpha;

         mu_alpha = MUr(:,1:R-1)*alphaPlus(1:R-1) + MUr(:,R)*(1-sum(alphaPlus(1:R-1)));
         C_alpha = sigma2r(1:R-1)*(alphaPlus(1:R-1).^2) + sigma2r(R)*(1-sum(alphaPlus(1:R-1)))^2;
         mu_alpha_star =MUr(:,k)*alpha_star + Muk*alphak + MUr(:,r)*(1-(alpha_star + sum(alphak)));
         C_alpha_star = sigma2r(k)*(alpha_star.^2) + Sk*(alphak.^2) + sigma2r(r)*(1-(alpha_star + sum(alphak)))^2;

         % difference between the logarithms of the distributions
         d = 0.5*(((norm(y-mu_alpha_star))^2)/C_alpha_star - ((norm(y-mu_alpha))^2)/C_alpha ) + (L/2)*log(C_alpha_star/C_alpha);

     else
         d = +Inf;
     end    
         
         % Acceptance-reject procedure to determine alphaPlus(t+1)

         if (d<0.0 | exp(-d)>rand)

             alphaPlus(k) = alpha_star;
             rho(k) = 1;
         else     
             rho(k) = 0;
         end  
 end
 
alphaPlus_out(comp_r) = alphaPlus(comp_r);
alphaPlus_out(r) = 1 - sum(alphaPlus(comp_r));
alphaPlus_out = alphaPlus_out';
