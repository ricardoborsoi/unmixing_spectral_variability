function [M_all] = libGeneratorSimpleAtmospheric(m0_sample,N_samples,desiredWavelenghths)
% generate spectral variability according to a simple illumination model, for different
% viewing angles
%
% INPUT: m0_sample - sample of the endmember spectra
%        N_samples - number of additional samples to generate (at different angles)
%        desiredWavelenghths - wavelength values for the bands of 'm0_sample'
%
% Ricardo Borsoi, 06/2020

% load spectral profile of direct and diffuse illumination sunlilght
load atmosIlluminationModel.mat

if exist('desiredWavelenghths','var')
    % Interpolate it to aviris bands
    Esun_direct_light  = interp1(wavelength_nm, Esun_direct_light,  desiredWavelenghths,'linear','extrap');
    Esky_diffuse_light = interp1(wavelength_nm, Esky_diffuse_light, desiredWavelenghths,'linear','extrap');
end

%N_samples = 2000;
%L = size(M_vegetation,1);
L = length(desiredWavelenghths);
M_all = zeros(L,N_samples);
aalpha0 = 0;
for i=1:N_samples
    %sample random angle
    aalpha_obs =  0.9*(pi/2) * rand; %0.6*rand; % 0.05 *(pi/2) * (rand-0.5) + 0;
    M_all(:,i) = m0_sample .* ((Esun_direct_light*cos(aalpha_obs)+Esky_diffuse_light)./(0.00005+Esun_direct_light*cos(aalpha0)+Esky_diffuse_light));
end 

