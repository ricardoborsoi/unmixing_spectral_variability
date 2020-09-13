function [M_var DM epsi] = generation_variabilite( M,nk )
%% generation_variabilite( M,nk ) 
%  Paper :  A. Halimi, P. Honeine, J. Bioucas-Dias, "Hyperspectral Unmixing
%          in Presence of Endmember Variability, Nonlinearity or Mismodelling
%          Effects", IEEE Trans. Image Process., 2016.
%  Code  :  Generation of EV
%  Version (April, 2016) by Abderrahim Halimi (a.halimi@hw.ac.uk)
%  For any comments contact the author
%
%% --------------- Inputs --------------------------------------------
%%% M      : endmembers of size (L x R)
%%% nk     : Number of pixels to be generated

%%%
%% --------------- Outputs --------------------------------------------
%%% M_var  : Pixel dependent endmembers (L x R x N)
%%% DM     : Pixel dependent variability  (L x R x N)
%%% epsi   : Variance coefficient of the EV  (1 x 1)
%% --------------------------------------------------------------------


[L R]= size(M); 
epsi = 0.001; 
 
Hdir     = toeplitz(exp(-[1:L].^2/100^2));
Kinv     = inv(Hdir);
Hdir     = chol(Kinv);
Kdir     = inv(Kinv); 

%%%% Fixed m_r in each region, eps_r for each element
for r=1:R 
    size(sqrt(epsi)* chol(Kdir)'*randn(L,nk))
        DM(:,:,r) = sqrt(epsi)* chol(Kdir)'*randn(L,nk);% L N R
end
DM     = permute(DM ,[1 3 2]);

Mrep       = repmat(M,[1 1 nk]);
M_var      = Mrep + DM; % L R N
ind        = find(M_var<0);
DM(ind)    = -Mrep(ind)+0.00001;
M_var      = Mrep + DM; % L R N

 
