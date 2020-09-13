function [sampledSigs]=libGeneratorHapke(inputSignature,N_samples)
% ==================================================
% Generate synthetic spectral signatures for different 
% observation angles, considering a Lambertian approximation.
% 
% 
% ==================================================


% First, regress the input spectra to obtain the albedo
% https://crustal.usgs.gov/speclab/data/HTMLmetadata/README/splib07_DataSeries.htm
mu = cosd(30);
mu0 = cosd(30);

% afun = @(w)( norm( x - w./ ((1+2*mu*sqrt(1-w)).*(1+2*mu0*sqrt(1-w))) )^2);

% regress to obtain the albedo under Lambretian assumption
x = inputSignature;
temp = ( ((mu0+mu)^2*(x.^2) + (1+4*mu0*mu*x).*(1-x)).^0.5 - (mu0+mu)*x ) ./ (1+4*mu0*mu*x);
w = -(temp.^2-1);


% generate signatures for different angles
sampledSigs = [];
rangevars1 = linspace(0,80,N_samples);
rangevars2 = linspace(0,80,N_samples);

% for k=0:15:60
% for k=0:15:80
for k1=rangevars1
for k2=rangevars2
    % mu_b = cosd(0);
    % mu0_b = cosd(0);
    mu_b = cosd(k1);
    mu0_b = cosd(k2);

    x2 = w./ ((1+2*mu_b*sqrt(1-w)).*(1+2*mu0_b*sqrt(1-w)));
    sampledSigs = [sampledSigs x2];
end
end


