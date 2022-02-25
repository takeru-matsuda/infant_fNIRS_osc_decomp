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
osc_mean = osc_mean(1:2*K,:,K);
osc_cov = osc_cov(1:2*K,1:2*K,:,K);
osc_phase = osc_phase(1:K,:,K);
%}
load('fNIRS_sample');
load('fNIRS_sample_result');
J = size(Y,1);
T = size(Y,2);

figure
subplot(2,1,1),plot((1:200)/fs,Y(1,1:200),'r','LineWidth',2),set(gca,'fontsize',16)
subplot(2,1,2),plot((1:200)/fs,Y(2,1:200),'b','LineWidth',2),set(gca,'fontsize',16)

figure
for k=1:K
    subplot(K,1,k)
    plot((1:200)/fs,osc_mean(2*k-1,1:200),'k','LineWidth',2)
    set(gca,'FontSize',16);
end

figure
for k=1:K
    subplot(K,1,k)
    plot((1:200)/fs,osc_c(2*k-1:2*k)*osc_mean(2*k-1:2*k,1:200),'k','LineWidth',2)
    set(gca,'FontSize',16);
end

figure
for k=1:K
    subplot(K,1,k)
    plot((1:200)/fs,osc_phase(k,1:200)/pi*180,'k','LineWidth',2)
    ylim([-180 180])
    yticks([-180 0 180])
    set(gca,'FontSize',16);
end
