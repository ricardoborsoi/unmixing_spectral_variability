% _______________________________________________________________________
% prospect_DB.m
% PROSPECT-Dynamic with brown pigments 
% PROSPECT version 6.0 (January, 16th 2017)
% subroutines required: calctav.m, dataSpec_PDB.m
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
% Plant leaf reflectance and transmittance are calculated from 400 nm to
% 2500 nm (1 nm step) with the following parameters:
%
%       - N   = leaf structure parameter
%       - Cab = chlorophyll a+b content in µg/cm²
%       - Car = carotenoids content in µg/cm²
%       - Anth = Anthocyanin content in µg/cm²
%       - Cbrown= brown pigments content in arbitrary units
%       - Cw  = equivalent water thickness in g/cm² or cm
%       - Cm  = dry matter content in g/cm²
%
% Here are some examples observed during the LOPEX'93 experiment on
% fresh (F) and dry (D) leaves :
%
% ---------------------------------------------
%                N     Cab     Cw        Cm    
% ---------------------------------------------
% min          1.000    0.0  0.004000  0.001900
% max          3.000  100.0  0.040000  0.016500
% corn (F)     1.518   58.0  0.013100  0.003662
% rice (F)     2.275   23.7  0.007500  0.005811
% clover (F)   1.875   46.7  0.010000  0.003014
% laurel (F)   2.660   74.1  0.019900  0.013520
% ---------------------------------------------
% min          1.500    0.0  0.000063  0.0019
% max          3.600  100.0  0.000900  0.0165
% bamboo (D)   2.698   70.8  0.000117  0.009327
% lettuce (D)  2.107   35.2  0.000244  0.002250
% walnut (D)   2.656   62.8  0.000263  0.006573
% chestnut (D) 1.826   47.7  0.000307  0.004305
% ---------------------------------------------
% _______________________________________________________________________
% _______________________________________________________________________%
% if no information about Anth or Cbrown and work on green / mature leaves
% set it to 0 													 		 %
% if no information about Car and work on green / mature leaves 		 %
% set the Chl / Car ratio between 4 and 5. this is not appropriate for 	 %
% senescent leaves 														 %
% _______________________________________________________________________%
% this code includes numerical optimizations proposed in the FLUSPECT code
% Authors: Wout Verhoef, Christiaan van der Tol (c.vandertol@utwente.nl) & 
% Joris Timmermans
% Date: 2007
% Update from PROSPECT to FLUSPECT: January 2011 (CvdT)
% for more info about FLUSPECT, see publication: 
% Vilfan, N., van der Tol, C., Muller, O., Rascher, U., Verhoef, W., 2016. 
% Fluspect-B: A model for leaf fluorescence, reflectance and transmittance 
% spectra. Remote Sens. Environ. 186, 596–615. doi:10.1016/j.rse.2016.09.017

function LRT=prospect_DB(N,Cab,Car,Ant,Brown,Cw,Cm)
% ***********************************************************************
% Jacquemoud S., Baret F. (1990), PROSPECT: a model of leaf optical
% properties spectra, Remote Sens. Environ., 34:75-91.
% Reference: 
% Féret, Gitelson, Noble & Jacqumoud (2017). PROSPECT-D: Towards modeling 
% leaf optical properties through a complete lifecycle
% Remote Sensing of Environment, 193:204–215
% DOI: http://doi.org/10.1016/j.rse.2017.03.004
% The specific absorption coefficient corresponding to brown pigment is
% provided by Frederic Baret (EMMAH, INRA Avignon, baret@avignon.inra.fr)
% and used with his autorization.
% ***********************************************************************

data    = dataSpec_PDB;
lambda  = data(:,1);    nr      = data(:,2);
Kab     = data(:,3);    Kcar    = data(:,4);
Kant    = data(:,5);    KBrown  = data(:,6);
Kw      = data(:,7);    Km      = data(:,8);
Kall    = (Cab*Kab+Car*Kcar+Ant*Kant+Brown*KBrown+Cw*Kw+Cm*Km)/N;
j       = find(Kall>0);               % Non-conservative scattering (normal case)
t1      = (1-Kall).*exp(-Kall);
t2      = Kall.^2.*expint(Kall);
tau     = ones(size(t1));
tau(j)  = t1(j)+t2(j);

% ***********************************************************************
% reflectance and transmittance of one layer
% ***********************************************************************
% Allen W.A., Gausman H.W., Richardson A.J., Thomas J.R. (1969),
% Interaction of isotropic ligth with a compact plant leaf, J. Opt.
% Soc. Am., 59(10):1376-1379.
% ***********************************************************************
% reflectivity and transmissivity at the interface
%-------------------------------------------------
talf    = calctav(40,nr);
ralf    = 1-talf;
t12     = calctav(90,nr);
r12     = 1-t12;
t21     = t12./(nr.^2);
r21     = 1-t21;

% top surface side
denom   = 1-r21.*r21.*tau.^2;
Ta      = talf.*tau.*t21./denom;
Ra      = ralf+r21.*tau.*Ta;

% bottom surface side
t       = t12.*tau.*t21./denom;
r       = r12+r21.*tau.*t;

% ***********************************************************************
% reflectance and transmittance of N layers
% Stokes equations to compute properties of next N-1 layers (N real)
% Normal case
% ***********************************************************************
% Stokes G.G. (1862), On the intensity of the light reflected from
% or transmitted through a pile of plates, Proc. Roy. Soc. Lond.,
% 11:545-556.
% ***********************************************************************
D       = sqrt((1+r+t).*(1+r-t).*(1-r+t).*(1-r-t));
rq      = r.^2;
tq      = t.^2;
a       = (1+rq-tq+D)./(2*r);
b       = (1-rq+tq+D)./(2*t);

bNm1    = b.^(N-1);                  %
bN2     = bNm1.^2;
a2      = a.^2;
denom   = a2.*bN2-1;
Rsub    = a.*(bN2-1)./denom;
Tsub    = bNm1.*(a2-1)./denom;

% Case of zero absorption
j       = find(r+t >= 1);
Tsub(j) = t(j)./(t(j)+(1-t(j))*(N-1));
Rsub(j)	= 1-Tsub(j);

% Reflectance and transmittance of the leaf: combine top layer with next N-1 layers
denom   = 1-Rsub.*r;
tran    = Ta.*Tsub./denom;
refl    = Ra+Ta.*Rsub.*t./denom;

LRT     = [lambda refl tran];
