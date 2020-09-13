function [sig2r_out]= sample_sigma2r(alpha_r, MUr, y, delta)

% Endmember variance

[L R] = size(MUr);

B = (norm(y - MUr*alpha_r)^2 + 2*sum(alpha_r.^2)*delta)./(2*sum(alpha_r.^2));

sig2inv = gamrnd(L/2+1,inv(B),1,1);
sig2r_out = 1/sig2inv ; 