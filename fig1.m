addpath('./osc_decomp')
T = 200;
K = 4;
fs = 100;
a = [0.9 0.9 0.9 0.9];
f = [5 10 20 30];
sigma2 = [.1 .1 .1 .1];
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
y = sum(x(1:2:end,:))+sqrt(tau2)*randn(1,T);

MAX_OSC = 7;
[osc_param,osc_AIC,osc_mean,osc_cov,osc_phase] = osc_decomp(y,fs,MAX_OSC);
[minAIC,K] = min(osc_AIC);
osc_a = osc_param(K,1:K);
osc_f = osc_param(K,K+1:2*K);
osc_sigma2 = osc_param(K,2*K+1:3*K);
osc_tau2 = osc_param(K,3*K+1);

figure
plot((1:T)/fs,y,'k','LineWidth',2)
xlabel('Time')
ylabel('y1')
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
    plot((1:T)/fs,osc_phase(k,:,K)/pi*180,'k','LineWidth',2)
    ylim([-180 180])
    yticks([-180 0 180])
    set(gca,'FontSize',16);
end
