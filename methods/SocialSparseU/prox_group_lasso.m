function Y = prox_group_lasso(X,groups,tau)

nbg = max(groups);

[P,N] = size(X);

Y = zeros(P,N);

for p = 1:nbg
   Y(groups == p,:) = vector_soft_col(X(groups == p,:),tau);
end