function w = ProjectOntoSimplex(v, b)
% PROJECTONTOSIMPLEX Projects point onto simplex of specified radius.
%
% w = ProjectOntoSimplex(v, b) returns the vector w which is the solution
%   to the following constrained minimization problem:
%
%    min   ||w - v||_2
%    s.t.  sum(w) = b, w >= 0.
%
%   That is, performs Euclidean projection of v to the positive simplex of
%   radius b.
%
% Author: John Duchi (jduchi@cs.berkeley.edu)
% Correction: Pierre-Antoine Thouvenin, February 12th 2015.

if (b < 0)
  error('Radius of simplex is negative: %2.3f\n', b);
end
% v = (v > 0) .* v; % faux et inutile
u = sort(v,'descend');
sv = cumsum(u);
rho = find(u > (sv - b) ./ (1:length(u))', 1, 'last');
% theta = max(0, (sv(rho) - b) / rho); % idem !
theta = (sv(rho) - b) / rho;
w = max(v - theta, 0);
