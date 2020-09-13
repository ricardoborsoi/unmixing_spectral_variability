% ********************************************************************************
% main routine to generate leaf optical properties using PROSPECT-D
% ********************************************************************************
% _______________________________________________________________________
% for any question or request, please contact: 
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
% _______________________________________________________________________
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

load('leaf_parameter.txt');
N     = leaf_parameter(1); Cab   = leaf_parameter(2);
Car   = leaf_parameter(3); Anth  = leaf_parameter(4); Cbrown= leaf_parameter(5);
Cw    = leaf_parameter(6); Cm    = leaf_parameter(7);

LRT=prospect_DB(N,Cab,Car,Anth,Cbrown,Cw,Cm);

dlmwrite('leaf_spectrum.txt',LRT,'delimiter','\t','precision',5)
