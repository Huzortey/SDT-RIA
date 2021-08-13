
%%  SECOND DERIVATIVE TECHNIQUE & RANGE INDEPENDENT ALGORITHM
%   By Andrew A. Huzortey, Benjamin Anderson & Alfred Owusu
%   8/13/2021
close all; clear; clc
%%  Import the Data into MatLab 
%   should be arranged with wavelength and Intesity 
%   as 1st and 2nd columns respectively
ay=importdata('Tourmaline_785_RAW_2.xls'); x=ay(:,1)'; y=ay(:,2)';
figure;plot(x,y); 
F1=y;  

%% OPTIONAL: select spectral range of interest (ROI) 
%  R0I1 is first & and ROI2 is last
%  ROI1=150; ROI2=1959;
%  F1=F1(ROI1:ROI2);  x=x(ROI1:ROI2);  figure; plot(x,F1);%return
%  Save Range of Interest (ROI) of Raw Spectrum

%% Compute Second Derivative 
%  derivxy function from from https://terpconnect.umd.edu isrecommended
SF1=derivxy(x,F1);     % First derivative
SF2=derivxy(x,SF1);    % Second derivative
df2=SF2;
%% OPTIONAL:
%  Smoothen to remove noise using 1st order savisty golay filter
%  Number of Iteration can be varied
i=0;
for i=1:3
    df2=sgolayfilt(df2,1,11);
    i=i+1;
end
figure; plot(x,df2)

%% Find Initial values and Input to Fit the Derivative Spectra
%  Initial values 
df3=df2(10:length(df2)-10); % Avoids first and last points  
[pks, locs]=max(-df3);      % Invert spectrum to identify minimum point easily
locs=locs+10;               % locs is position of maximum point
vv=[12 x(locs)];            %% 12/cm is the general width of a raman peak

% Second derivative Lorentzian function
Ao3=df2./((((8*(x-vv(2)).^2)./(((vv(1).^4).*((1+(((x-vv(2))./vv(1))... 
    .^2)).^3)))))-(2./((vv(1).^2).*((1+(((x-vv(2))./vv(1)).^2)).^2))));

v=[ Ao3(locs) vv(1) x(locs)];
F=@(v,x)v(1).*((((8*(x-v(3)).^2)./(((v(2).^4).*((1+(((x-v(3))./v(2)).^2)).^3)))))...
    -(2./((v(2).^2).*((1+(((x-v(3))./v(2)).^2)).^2))));
lsqOpts = optimoptions('lsqcurvefit', 'MaxFunEvals', 1000000);
r1=lsqcurvefit(F,v,x,df2, [], [], lsqOpts);
figure; plot(x,F(r1,x),'r'); hold on; plot(x,df2)

r=lsqcurvefit(F,r1,x,df2, [], [], lsqOpts);
figure; plot(x,F(r,x),'r'); hold on; plot(x,df2)
% 
%% modified iterative fitting using a savitzy golay filter of order zero (lowpass filter)
%%
G6=F1;          % The original data
CC2=F1(:,locs); % Amplitude value at the maximum peak position in original data
CC3=CC2-r(1);   % Amplitude of baseline estimated  
hi_diff=0;      % Iterarions end when hi_diff is equal to r(1)
tic
while hi_diff< r(1)
 G6=sgolayfilt(G6,1,11); 
 G6(G6>F1)=F1(F1<G6);
 Q1=G6(:,locs);
 hi_diff=CC2-Q1; % counts the number of iterations
 r(1);
end
figure
subplot(2,1,1);plot(x,F1,'b',x,G6,'--r','linewidth',1.5)
legend('Original Specta','Estimated baseline')
%% Subtract Baseline from Original Spectrum to obtain Raman Spectra
G7=F1-G6;
subplot(2,1,2);plot(x,G7,'b','linewidth',1.5)
legend('Recovered Raman Spectra')
