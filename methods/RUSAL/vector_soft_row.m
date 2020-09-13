function Y = vector_soft_row(X,tau)
%
%  computes the vector soft columnwise


NU = sqrt(sum(X.^2,2));
A = max(0, NU-tau);
Y = repmat((A./(A+tau)),1,size(X,2)).* X;
