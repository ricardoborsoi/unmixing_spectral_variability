% // ====================================================================
% // This file is part of the Endmember Induction Algorithms Toolbox for MATLAB 
% // Copyright (C) Grupo de Inteligencia Computacional, Universidad del 
% // País Vasco (UPV/EHU), Spain, released under the terms of the GNU 
% // General Public License.
% //
% // Endmember Induction Algorithms Toolbox is free software: you can redistribute 
% // it and/or modify it under the terms of the GNU General Public License 
% // as published by the Free Software Foundation, either version 3 of the 
% // License, or (at your option) any later version.
% //
% // Endmember Induction Algorithms Toolbox is distributed in the hope that it will
% // be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
% // of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
% // General Public License for more details.
% //
% // You should have received a copy of the GNU General Public License
% // along with Endmember Induction Algorithms Toolbox. 
% // If not, see <http://www.gnu.org/licenses/>.
% // ====================================================================
%
%% results = EIA_VCA(data,runs,p,verbosity)
%
% Miguel Angel Veganzones <miguel-angel.veganzones[AT]gipsa-lab.fr>
% GIPSA-lab, Grenoble-INP
% http://www.gipsa-lab.grenoble-inp.fr
%
% Manuel Grana <manuel.grana[AT]ehu.es>
% Grupo de Inteligencia Computacional (GIC), Universidad del Pais Vasco /
% Euskal Herriko Unibertsitatea (EHU/UPV)
% http://www.ehu.es/computationalintelligence
% 
% Copyright (2011,2012) Grupo de Inteligencia Computacional @ Universidad del Pais Vasco, Spain.
% Copyright (2005)  Instituto de Telecomunicações @ Instituto Superior Técnico Lisboa, Portugal.
%
% Vertex Component Analysis endmembers induction algorithm.
% ------------------------------------------------------------------------------
% Input:   data      : column data matrix [nvariables x nsamples]
%          runs      : number of times the EIA is run over the data
%          p         : number of endmembers to be induced.
%          verbosity : true (default) | false
%
% Output:  results   : an struct with the variables described above. Some variables are
%                      represented as a cell with the a size equal to the
%                      number of times the EIA algorith is run ('runs'
%                      parameter)
%                      * E   : set of induced endmembers [nvariables x p]
%                      * C   : induced endmembers indexes vector [nsamples] with {0,1} 
%                              values, where '1' indicates that the corresponding sample
%                              has been identified as an endmember. Some of the
%                              algorithms do not select pixels as the endmembers
%                              and, in such cases C is empty.
%                      * t   : the time in seconds that took the EIA to finish
%
% Bibliographical references:
% [1] Jose M. P. Nascimento and Jose M. B. Dias "Vertex Component Analysis: A Fast Algorithm to Unmix Hyperspectral Data", IEEE Trans. Geosci. Remote Sensing, vol. 43, no. 4, pp. 898-910, 2005. 
function results = vca(data,runs,p,verbosity)

%% Parameters
if nargin < 3
    error('Insufficient parameters');
end
if ~isnumeric(runs) || runs <= 0
    error('runs must be a positive integer');
end
if ~isnumeric(p) || p <= 0
    error('p must be a positive integer');
end
if nargin < 4
    verbosity = true;
else
    if ~islogical(verbosity)
        error('verbosity must be a logical value');
    end
end
if verbosity
    disp('Running VCA algorithm ...');
end

%% data size
[~,nsamples] = size(data);

%% SNR Estimates
if verbosity
    disp('... estimating SNR');
end
r_m = mean(data,2);      
R_m = repmat(r_m,[1 nsamples]); % mean of each band
R_o = data - R_m;           % data with zero-mean 
[Ud,~,~] = svd(R_o*R_o'/nsamples);  % computes the p-projection matrix 
Ud = Ud(:,[1:p]);
x_p = Ud'*R_o;                 % project the zero-mean data onto p-subspace
SNR = estimate_snr(data,r_m,x_p);
SNR_th = 15 + 10*log10(p);
         
%% Choosing Projective Projection or projection to p-1 subspace
if SNR < SNR_th,   
    if verbosity
        disp('... projective projection');
    end
    d = p-1;
    Ud= Ud(:,1:d);    
    Rp =  Ud * x_p(1:d,:) + repmat(r_m,[1 nsamples]);      % again in dimension L
    x = x_p(1:d,:);             %  x_p =  Ud' * R_o; is on a p-dim subspace
    c = max(sum(x.^2,1))^0.5;
    y = [x ; c*ones(1,nsamples)] ;
else
    if verbosity
        disp('... projection to p-1 subspace');
    end
    d = p;
    [Ud,~,~] = svd(data*data'/nsamples);         % computes the p-projection matrix 
    Ud = Ud(:,[1:p]);
    x_p = Ud'*data;
    Rp =  Ud * x_p(1:d,:);      % again in dimension L (note that x_p has no null mean)
    x =  Ud' * data;
    u = mean(x,2);        %equivalent to  u = Ud' * r_m
    y =  x./ repmat( sum( x .* repmat(u,[1 nsamples]) ) ,[d 1]);
end
 
%% Initialization
results = struct('E',cell(runs,1),'C',cell(runs,1),'t',cell(runs,1));

%% runs
for r = 1:runs
    if verbosity
        fprintf('... run %g\n',r);
    end
    tic;
    C = zeros(1,p);
    A = zeros(p,p);
    A(p,1) = 1;
    for i=1:p
          w = rand(p,1);   
          f = w - A*pinv(A)*w;
          f = f / sqrt(sum(f.^2));
          v = f'*y;
          [~,C(i)] = max(abs(v));
          A(:,i) = y(:,C(i));        % same as x(:,C(i))
    end
    E = Rp(:,C);
    t = toc;
    if verbosity
        fprintf('... finished in %g seconds\n',t);
    end
    results(r,1).E = E;
    results(r,1).C = C;
    results(r,1).t = t;
end
    
%% Internal functions
function snr_est = estimate_snr(R,r_m,x)

[L,~]=size(R);           % L number of bands (channels)
                         % N number of pixels (Lines x Columns) 
[p,N]=size(x);           % p number of endmembers (reduced dimension)
P_y = sum(R(:).^2)/N;
P_x = sum(x(:).^2)/N + r_m'*r_m;
snr_est = 10*log10( (P_x - p/L*P_y)/(P_y- P_x) );
return;