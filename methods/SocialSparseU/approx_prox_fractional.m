function Y = approx_prox_fractional(X,tau,p)
   Y = max(abs(X)-tau^(2-p)*abs(X).^(p-1),0).*sign(X);   
end
