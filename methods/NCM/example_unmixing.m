clc
clear all
close all

% loading: endmember spectra
spectra = load ('spectra.txt');
spectra = spectra-min(spectra(:));
[L tmp] = size(spectra);
N_spectra = tmp-1;

% wavelength
wave = spectra(:,1);

% matrix of all the spectra
M_all = spectra(:,2:end);

% step for the number of bands
p=1;

% vector of the endmember chosen for the mixture
endm_chosen = [1 2];
M = M_all(1:p:end,endm_chosen);

% L = number of bands
% R = number of endmembers 
[L R] = size(M);

figure
subplot(2,1,1)
for i=1:R
    hold on
    switch i 
        case 1
            motif = '-';
        case 2
            motif = '--';
        case 3
            motif = ':';
        case 4
            motif = '-:';
    end
    set(gca,'box','on')
    set(gca,'fontsize',20)
    plot(wave,M(:,i),motif,'linewidth',1.6)    
    hold off
end
ylabel('Reflectance')
xlabel('Wavelength')
axis([0.3 2.6 0 1.1*max(M(:,2))])

% number of pixels
N_pixels = 1;

alphaPlus = [0.3 0.7];

y = zeros(L,N_pixels);
Mplus = zeros(L,R);

% endmember variance
s = 0.005; 


% resulting signal
% endmember vector construction
for r=1:R
    Mplus(:,r)=multrandn(M(:,r), s*eye(L),1);    
end   
y = Mplus*alphaPlus';

SNR = 10*log10((norm(M*alphaPlus')^2)/(L*s*(sum(alphaPlus.^2))))

%%

% Markov chain length
Nmc = 10000;

% number of burn-in period
Nbi = 1000;

% Unmixing procedure
tic
[TalphaPlus,Tsigma2r] = unmixing(y,M,Nmc,Nbi);
toc
%%
% exploitation

% abundances
 result = exploitation(TalphaPlus,Nbi,Nmc); 
 %variance
 result = exploitation_sigma(Tsigma2r,Nbi,Nmc);
