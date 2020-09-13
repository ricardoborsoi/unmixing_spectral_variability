function [Y,Yim,A,A_cube,Mth,M_avg]=generate_image(SNR)
% ------------------------------------
% Generates a hyperspectral image following the linear mixing model
% with spectral variability, for a given SNR. The spectral variability
% is geenrated accorging to simplifications of radiative transfer models.
% The mixed images contain vegetation, soil and water as materials.
%
% INPUT: SNR (in decibels)
%
% OUTPUT: Y - mixed image ordered as a matrix (bands * pixels)
%         Yim - mixed image ordered as a tensor (rows * cols * bands)
%         A - abundances ordered as a matrix (#ems * pixels)
%         A_cube - abundances ordered as a tensor (rows * cols * #ems)
%         Mth - tensor containign the true endmember matrix for each pixel (bands * #ems * pixels)
%         M_avg - average endmember matrix
%
% Ricardo Borsoi, 06/2020
% ------------------------------------


% load synthetic abundance maps (generated according to a Gaussian random fields model)
load('abundance_spatial_3em.mat')
[nr,nc,P] = size(A_cube);
N = nr*nc;


% vegetation, soil, water
% load reference endmembers
load('endmembers_b.mat')
load('EM_b_SlectBands.mat')
load('AVIRIS_wavelengths.mat')

M0_water = M(:,2);
M0_dirt  = M(:,3);


reference_wavelenths = AVIRIS_wavlen(SlectBands);


% generate spectral signatures of vegetation using the PROSPECT-D model
% N_samples = 5;
% paramOption = 'all';
% [M_vegetation,~]=libGeneratorPROSPECTD(paramOption, N_samples, reference_wavelenths);


N_samples = 10;
M_vegetation = [];
paramOption = 'structure';
[resv,~]=libGeneratorPROSPECTD(paramOption, N_samples, reference_wavelenths);
M_vegetation = [M_vegetation,resv];
paramOption = 'chlorophyll';
[resv,~]=libGeneratorPROSPECTD(paramOption, N_samples, reference_wavelenths);
M_vegetation = [M_vegetation,resv];
paramOption = 'waterthickness';
[resv,~]=libGeneratorPROSPECTD(paramOption, N_samples, reference_wavelenths);
M_vegetation = [M_vegetation,resv];
paramOption = 'drymatter';
[resv,~]=libGeneratorPROSPECTD(paramOption, N_samples, reference_wavelenths);
M_vegetation = [M_vegetation,resv];



% generate spectral signatures of Dirt following a simplified version of Hapke's model
N_samples = 10;
[M_dirt]=libGeneratorHapke(M0_dirt,N_samples);



% generate spectral signatures of Water using a simplified illumination model

% load spectral profile of direct and diffuse illumination sunlilght
load atmosIlluminationModel.mat

% Interpolate it to aviris bands
Esun_direct_light  = interp1(wavelength_nm, Esun_direct_light,  AVIRIS_wavlen(SlectBands),'linear','extrap');
Esky_diffuse_light = interp1(wavelength_nm, Esky_diffuse_light, AVIRIS_wavlen(SlectBands),'linear','extrap');

N_samples = 200;
L = size(M_vegetation,1);
M_water = zeros(L,N_samples);
aalpha0 = 0;
for i=1:N_samples
    %sample random angle
    aalpha_obs =  0.9*(pi/2) * rand; %0.6*rand; % 0.05 *(pi/2) * (rand-0.5) + 0;
    M_water(:,i) = M0_water .* ((Esun_direct_light*cos(aalpha_obs)+Esky_diffuse_light)./(0.00005+Esun_direct_light*cos(aalpha0)+Esky_diffuse_light));
end 


% -------------------------------------------------
%{
figure
% [ha, pos] = tight_subplot(1, 3, [0.1 0.1], 0.1, 0.1);
% axes(ha(1));
subplot(1,3,1)
plot(AVIRIS_wavlen(SlectBands)/1000, M_vegetation)
xlim([0.4 2.457]), %ylim([0 0.9])
xlabel('Wavelength [$\mu$m]','interpreter','latex','fontsize',14)
ylabel('Reflectance','interpreter','latex','fontsize',14)
% axes(ha(2));
subplot(1,3,2)
plot(AVIRIS_wavlen(SlectBands)/1000, M_dirt)
xlim([0.4 2.457]), %ylim([0 0.9])
xlabel('Wavelength [$\mu$m]','interpreter','latex','fontsize',14)
% axes(ha(3));
subplot(1,3,3)
plot(AVIRIS_wavlen(SlectBands)/1000, M_water)
xlim([0.4 2.457]), %ylim([0 0.9])
xlabel('Wavelength [$\mu$m]','interpreter','latex','fontsize',14)
%}

% -------------------------------------------------
A = reshape(A_cube,N,P)';
Mth = zeros(L,P,N);
for i=1:N

    Mth(:,1,i) = M_vegetation(:, randi(size(M_vegetation,2),1,1));
    Mth(:,2,i) = M_dirt(:, randi(size(M_dirt,2),1,1));
    Mth(:,3,i) = M_water(:, randi(size(M_water,2),1,1));
    
    Y(:,i) = Mth(:,:,i) * A(:,i);
end

M_avg = [mean(M_vegetation,2),...
         mean(M_dirt,2),...
         mean(M_water,2)];

%----------------------
% add noise
pw_signal = norm(Y,'fro')^2/(prod(size(Y)));
npower    = pw_signal/(10^(SNR/10));
Y = Y + sqrt(npower) * randn(size(Y));
%----------------------

Yim = reshape(Y',[nr nc L]);






