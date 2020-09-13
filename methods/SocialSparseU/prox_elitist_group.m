function Y = prox_elitist_group(X,groups,tau)

nbg = max(groups);

[P,N] = size(X);

Y = zeros(P,N);

for p = 1:nbg
%     for k = 1:N
%    tau_p = tau/(1+tau) * sum(abs(X(groups == p,k))); 
%    Y(groups == p,k) = soft(X(groups == p,k),tau_p);
%     end
   curr_X = X(groups == p,:);
   tau_p = tau/(1+tau) * sum(abs(curr_X)); 
   Y(groups == p,:) = max(abs(curr_X)-repmat(tau_p,sum(groups == p),1),0).*sign(curr_X ); 
    

end

