%{
addpath('./osc_decomp')
load('fNIRS_sample.mat')
J = size(Y,1);
T = size(Y,2);
MAX_OSC = 7;
MAX_VAR = 10;
[osc_param,osc_AIC,osc_mean,osc_cov,osc_phase] = osc_decomp(Y,fs,MAX_OSC,MAX_VAR);
[minAIC,K] = min(osc_AIC);
osc_a = osc_param(K,1:K);
osc_f = osc_param(K,K+1:2*K);
osc_sigma2 = osc_param(K,2*K+1:3*K);
osc_c = osc_param(K,3*K+1:(2*J+1)*K);
osc_tau2 = osc_param(K,(2*J+1)*K+1);

periodogram1 = zeros(1,T/2);
spectrum1 = zeros(K,T/2);
periodogram2 = zeros(1,T/2);
spectrum2 = zeros(K,T/2);
for j=1:T/2
    periodogram1(j) = abs(Y(1,:)*exp(-1i*2*pi*j/T*[1:T])')^2/2/pi/T;
    periodogram2(j) = abs(Y(2,:)*exp(-1i*2*pi*j/T*[1:T])')^2/2/pi/T;
    for k=1:K
        a = osc_a(k);
        theta = 2*pi*osc_f(k)/fs;
        A = (1-2*a^2*cos(theta)^2+a^4*cos(2*theta))/a/(a^2-1)/cos(theta);
        b = (A-2*a*cos(theta)+sqrt((A-2*a*cos(theta))^2-4))/2;
        spectrum1(k,j) = -osc_sigma2(k)*a*cos(theta)/b*abs(1+b*exp(-1i*2*pi*j/T))^2/abs(1-2*a*cos(theta)*exp(-1i*2*pi*j/T)+a^2*exp(-1i*4*pi*j/T))^2/2/pi;
        spectrum2(k,j) = -norm(osc_c(2*k-1:2*k))^2*osc_sigma2(k)*a*cos(theta)/b*abs(1+b*exp(-1i*2*pi*j/T))^2/abs(1-2*a*cos(theta)*exp(-1i*2*pi*j/T)+a^2*exp(-1i*4*pi*j/T))^2/2/pi;
    end
end
%}
load('fNIRS_sample');
load('fNIRS_sample_result');
J = size(Y,1);
T = size(Y,2);

figure
spectrogram(Y(1,:),32,16,[0:64]/128*fs,fs,'yaxis'),set(gca,'FontSize',20);

figure
spectrogram(Y(2,:),32,16,[0:64]/128*fs,fs,'yaxis'),set(gca,'FontSize',20);

figure,hold on
plot([1:T/2]/T*fs,log10(periodogram1),'k--')
plot([1:T/2]/T*fs,log10(sum(spectrum1)),'b*-')
for k=1:K
    plot([1:T/2]/T*fs,log10(spectrum1(k,:)),'r-')
end
xlabel('frequency (Hz)','FontSize',20)
ylabel('log10 power','FontSize',20)
set(gca,'FontSize',16);

figure,hold on
plot([1:T/2]/T*fs,log10(periodogram2),'k--')
plot([1:T/2]/T*fs,log10(sum(spectrum2)),'b*-')
for k=1:K
    plot([1:T/2]/T*fs,log10(spectrum2(k,:)),'r-')
end
xlabel('frequency (Hz)','FontSize',20)
ylabel('log10 power','FontSize',20)
set(gca,'FontSize',16);
