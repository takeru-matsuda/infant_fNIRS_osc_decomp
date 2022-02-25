addpath('./osc_decomp')
J = 2;
T = 200;
K = 3;
fs = 200;
a = [0.9 0.85 0.8];
f = [20 50 75];
sigma2 = [.1 .1 .1];
tau2 = 0.1;
x = zeros(2*K,T);
for k=1:K
    x(2*k-1:2*k,1) = sqrt(sigma2(k))/sqrt(1-a(k)^2)*randn(2,1);
end
for t=2:T
    for k=1:K
        x(2*k-1:2*k,t) = a(k)*[cos(2*pi*f(k)/fs) -sin(2*pi*f(k)/fs); sin(2*pi*f(k)/fs) cos(2*pi*f(k)/fs)]*x(2*k-1:2*k,t-1)+sqrt(sigma2(k))*randn(2,1);
    end
end
c = zeros(J,2*K);
c(1,1:2:2*K) = 1;
c(2,:) = [1.8 1.2 0 1 -0.1 0.4];
Y = c*x+sqrt(tau2)*randn(J,T);

MAX_OSC = 7;
[osc_param,osc_AIC,osc_mean,osc_cov,osc_phase] = osc_decomp(Y,fs,MAX_OSC);
[minAIC,K] = min(osc_AIC);
osc_a = osc_param(K,1:K);
osc_f = osc_param(K,K+1:2*K);
osc_sigma2 = osc_param(K,2*K+1:3*K);
osc_c = osc_param(K,3*K+1:(2*J+1)*K);
osc_tau2 = osc_param(K,(2*J+1)*K+1);

figure
subplot(2,1,1)
plot((1:T)/fs,Y(1,:),'k','LineWidth',2)
set(gca,'FontSize',16);
subplot(2,1,2)
plot((1:T)/fs,Y(2,:),'k','LineWidth',2)
set(gca,'FontSize',16);

figure
for k=1:K
    subplot(K,1,k)
    plot((1:T)/fs,osc_mean(2*k-1,:,K),'k','LineWidth',2)
    set(gca,'FontSize',16);
end

figure
for k=1:K
    subplot(K,1,k)
    plot((1:T)/fs,osc_c(2*k-1:2*k)*osc_mean(2*k-1:2*k,:,K),'k','LineWidth',2)
    set(gca,'FontSize',16);
end

figure
for k=1:K
    subplot(K,1,k)
    plot((1:T)/fs,osc_phase(k,:,K)/pi*180,'k','LineWidth',2)
    set(gca,'FontSize',16);
    ylim([-180 180])
    yticks([-180 0 180])
end
