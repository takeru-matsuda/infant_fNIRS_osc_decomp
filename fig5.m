addpath('./osc_decomp')
J = 2;
T = 1800;
K = 3;
fs = 10;
a = [0.99 0.4 0.99];
f = [.03 2 2];
sigma2 = (1-a.^2)*10^-2;
tau2 = 10^-6;
x = zeros(2*K,T);
phi = zeros(K,T);
for k=1:K
    x(2*k-1:2*k,1) = sqrt(sigma2(k))/sqrt(1-a(k)^2)*randn(2,1);
end
for t=2:T
    for k=1:K
        x(2*k-1:2*k,t) = a(k)*[cos(2*pi*f(k)/fs) -sin(2*pi*f(k)/fs); sin(2*pi*f(k)/fs) cos(2*pi*f(k)/fs)]*x(2*k-1:2*k,t-1)+sqrt(sigma2(k))*randn(2,1);
        phi(k,t) = atan2(x(2*k,t),x(2*k-1,t));
    end
end
c = zeros(J,2*K);
c(1,1:2:2*K) = 1;
c(2,:) = [0.5*cos(230/180*pi) 0.5*sin(230/180*pi) -1 0 0.2 0];
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
plot((1:200)/fs,Y(1,1:200),'r-','LineWidth',2)
set(gca,'FontSize',16);
ylabel('y1','FontSize',20)
subplot(2,1,2)
plot((1:200)/fs,Y(2,1:200),'b-','LineWidth',2)
set(gca,'FontSize',16);
ylabel('y2','FontSize',20)
xlabel('Time (s)','FontSize',20)

figure
for k=1:K
    subplot(K,1,k)
    plot((1:200)/fs,osc_mean(2*k-1,1:200,K),'k-','LineWidth',2)
    set(gca,'FontSize',16);
    ylabel(texlabel(sprintf('x^{(%d)}',k)),'FontSize',20)
end
xlabel('Time (s)','FontSize',20)
