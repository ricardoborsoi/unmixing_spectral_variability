function [A,M,alg_time,Yhat] = adaptor_BCM(Yim,M0,Lib,A_init,opt)


if ~exist('opt','var')
    % tune regularization parameters
else
end


[nr,nc,L] = size(Yim);
Yim(Yim<0) = 0;
Yim(Yim>1) = 1;

Y = reshape(Yim,nr*nc,L)';

for ii=1:size(M0,2)
    Lib{ii}(Lib{ii}<1e-6) = (1e-6)*rand;
    Lib{ii}(Lib{ii}>1-1e-6) = 1-(1e-6)*rand;
end
for ii=1:size(M0,2)
    for jj=1:L
        if var(Lib{ii}(jj,:)) < 1e-6 
            Lib{ii}(jj,:) = Lib{ii}(jj,:) + (1e-7)*rand(1,size(Lib{ii},2));
        end
    end
end


cd methods/BCM

Xim = Yim;


% convert lib from Px1 to a 1xP cell:
Lib2 = cell(1,size(M0,2));
for  ii=1:size(M0,2)
    Lib2{1,ii} = Lib{ii};
end


tic

[Parameters] = BCMParameters(Lib2);


%% 	BCM Unmixing
%%% Approximate running time: 70 seconds per BCM approach*. 
%%% *Measured on a desktop PC with Intel i7 3.20 GHz processor and 12 GB RAM.

%disp('BCM-Spectral-QP Unmixing...');
[Prop] = BCM(Xim, Parameters, 1);%BCM-Spectral-QP
%P1r = reshape(P1,[13 19 4]);

%disp('BCM-Spectral-MH Unmixing...');
%[Prop] = BCM(Xim, Parameters, 2);%BCM-Spectral-MH
%P2r = reshape(P2,[13 19 4]);

%disp('BCM-Spatial-QP Unmixing...');
%[Prop] = BCM(Xim, Parameters, 3);%BCM-Spatial-QP
%P3r = reshape(P3,[13 19 4]);

%disp('BCM-Spatial-MH Unmixing...');
%[Prop] = BCM(Xim, Parameters, 4);%BCM-Spatial-MH
%P4r = reshape(P4,[13 19 4]);

%A_BCM = reshape(Prop,[nr nc size(M0,2)]);
A_BCM = Prop';

alg_time = toc;

A = A_BCM;
M = [];

Yhat = M0 * A;

cd ../..




