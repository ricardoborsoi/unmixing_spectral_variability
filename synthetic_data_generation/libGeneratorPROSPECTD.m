function [resv,wavelengths]=libGeneratorPROSPECTD(paramOption, N_samples,desiredWavelenghths)
% ============================================================
% This code generates a set of vegetation spectral signatures 
% containing spectral variability ugint eh PROSPECT-D model
% by varying biophysicla parameters
% 
% INPUT
%   paramOption: type of parameter to sample from, can be
%                'structure', 'chlorophyll', 'waterthickness', or 'drymatter'
%   N_samples: number of samples
%   desiredWavelenghths: optional parameter, desired wavelengths of the signatures
%
% OUTPUT
%   resv: sampled spectral signatures
%   wavelengths: corresponding spectral wavelengths in nanometers
% ============================================================

if ~exist('N_samples','var')
    N_samples = 10;
end

% ********************************************************************************
% This code is directly based on the routine to 
% generate leaf optical properties using PROSPECT-D ********************************************************************************
% for any question or request about PROSPECT-D, please contact: 
%
% Jean-Baptiste FERET
% UMR-TETIS, IRSTEA Montpellier
% Maison de la Télédétection
% 500 rue Jean-Fracois Breton
% 34093 Montpellier cedex 5
% E-mail: jb.feret@teledetection.fr
%
% Stéphane JACQUEMOUD
% Université Paris Diderot / Institut de Physique du Globe de Paris
% 35 rue Hélène Brion
% 75013 Paris, France
% E-mail: jacquemoud@ipgp.fr
%
% http://teledetection.ipgp.fr/prosail/
% ********************************************************************************
% ********************************************************************************
% Main script running with PROSPECT-DB version 6.0 (16 January 2017)
% ********************************************************************************

% load biochemistry and structure variables
% N      = 1.2		% structure coefficient
% Cab    = 30.		% chlorophyll content (µg.cm-2) 
% Car    = 10.		% carotenoid content (µg.cm-2)
% Anth   = 1.		% Anthocyanin content (µg.cm-2)
% Cbrown = 0.0		% brown pigment content (arbitrary units)
% Cw     = 0.015	% EWT (cm)
% Cm     = 0.009	% LMA (g.cm-2)

% load standard lead parameters
% load('PROSPECT-D_Matlab/leaf_parameter.txt');
load('leaf_parameter.txt');
N     = leaf_parameter(1); 
Cab   = leaf_parameter(2);
Car   = leaf_parameter(3); 
Anth  = leaf_parameter(4); 
Cbrown= leaf_parameter(5);
Cw    = leaf_parameter(6); 
Cm    = leaf_parameter(7);



switch paramOption

case 'structure'
    % N      = 1.2		% structure coefficient
    % rangevar1 = 1:0.25:3.5; % N   = leaf structure parameter
    rangevar1 = linspace(1,2,N_samples);
    for iidx = 1:length(rangevar1)
        N = rangevar1(iidx);
        LRT=prospect_DB(N,Cab,Car,Anth,Cbrown,Cw,Cm);
        if iidx == 1
            resv	= LRT(:,2);
        else
            resv	= [resv LRT(:,2)];
        end
    end
    wavelengths = LRT(:,1); % in nm


case 'chlorophyll'
    % Cab    = 30.		% chlorophyll content (µg.cm-2) 
    % rangevar1 = [10:5:80]; % chlorophyll content (µg.cm-2) 
    rangevar1 = linspace(10,80,N_samples);
    for iidx = 1:length(rangevar1)
        Cab = rangevar1(iidx);
        LRT=prospect_DB(N,Cab,Car,Anth,Cbrown,Cw,Cm);
        if iidx == 1
            resv	= LRT(:,2);
        else
            resv	= [resv LRT(:,2)];
        end
    end
    wavelengths = LRT(:,1); % in nm


case 'waterthickness'
    % Cw     = 0.015	% EWT (cm)
    % 0.000063 .. 0.04
    % rangevar1 = [0.000005,0.00001,0.00005, 0.0001, 0.0005, 0.001, 0.005, 0.01, 0.03]; % Cw  = equivalent water thickness in g/cm² or cm
    rangevar1 = logspace(log10(0.005),log10(0.03),N_samples);
    for iidx = 1:length(rangevar1)
        Cw = rangevar1(iidx);
        LRT=prospect_DB(N,Cab,Car,Anth,Cbrown,Cw,Cm);
        if iidx == 1
            resv	= LRT(:,2);
        else
            resv	= [resv LRT(:,2)];
        end
    end
    wavelengths = LRT(:,1); % in nm


case 'drymatter'
    % Cm     = 0.009	% LMA (g.cm-2)
    % 0..1
    % rangevar1 = 0.001:0.002:0.022; %  Cm  = dry matter content in g/cm²
    rangevar1 = linspace(0.001,0.022,N_samples);
    for iidx = 1:length(rangevar1)
        Cm = rangevar1(iidx);
        LRT=prospect_DB(N,Cab,Car,Anth,Cbrown,Cw,Cm);
        if iidx == 1
            resv	= LRT(:,2);
        else
            resv	= [resv LRT(:,2)];
        end
    end
    wavelengths = LRT(:,1); % in nm


case 'all'
    if N_samples^4 > 50000, error('too many samples, reduce N_samples.'); end
    rangevar1 = linspace(1,2,N_samples); % N   = leaf structure parameter
    rangevar2 = linspace(10,80,N_samples); % chlorophyll content (µg.cm-2) 
    rangevar3 = logspace(log10(0.005),log10(0.03),N_samples); % Cw  = equivalent water thickness in g/cm² or cm
    rangevar4 = linspace(0.001,0.022,N_samples); %  Cm  = dry matter content in g/cm²
    resv = [];
    for iidx = 1:length(rangevar1)
    for jidx = 1:length(rangevar2)
    for kidx = 1:length(rangevar3)
    for lidx = 1:length(rangevar4)
        N   = rangevar1(iidx);
        Cab = rangevar2(jidx);
        Cw  = rangevar3(kidx);
        Cm  = rangevar4(lidx);
        LRT = prospect_DB(N,Cab,Car,Anth,Cbrown,Cw,Cm);
        resv = [resv LRT(:,2)];
        wavelengths = LRT(:,1); % in nm
    end
    end
    end
    end


otherwise 
    error('unknown option for PROSPECT-D library generation.')

end


% if desired wavelengths are provided, interpolate the signatures
tot_samples = size(resv,2);
if exist('desiredWavelenghths','var')
    resvTmp = resv;
    resv = zeros(length(desiredWavelenghths),tot_samples);
    for iidx=1:tot_samples
        resv(:,iidx) = interp1(wavelengths,resvTmp(:,iidx),desiredWavelenghths,'linear','extrap');
    end
    wavelengths = desiredWavelenghths;
end



