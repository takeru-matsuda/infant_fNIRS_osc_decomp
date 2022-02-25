addpath('./osc_decomp')
J = 2;
T = 1800;
K = 3;
fs = 10;
a = [0.99 0.4 0.99];
f = [.03 2 2];
sigma = (1-a.^2)*10^-2;
x = zeros(2*K,T);
phi = zeros(K,T);
for k=1:K
    x(2*k-1:2*k,1) = sqrt(sigma(k))/sqrt(1-a(k)^2)*randn(2,1);
end
for t=2:T
    for k=1:K
        x(2*k-1:2*k,t) = a(k)*[cos(2*pi*f(k)/fs) -sin(2*pi*f(k)/fs); sin(2*pi*f(k)/fs) cos(2*pi*f(k)/fs)]*x(2*k-1:2*k,t-1)+sqrt(sigma(k))*randn(2,1);
        phi(k,t) = atan2(x(2*k,t),x(2*k-1,t));
    end
end
c = zeros(J,2*K);
c(1,1:2:2*K) = 1;
c(2,:) = [0.5*cos(230/180*pi) 0.5*sin(230/180*pi) -1 0 0.2 0];

sigman = 2.^[-30:-1];
MAX_OSC = 7;
osc_param = zeros(MAX_OSC,(2*J+1)*MAX_OSC+1,length(sigman));
osc_AIC = zeros(MAX_OSC,length(sigman));
osc_mean = zeros(2*MAX_OSC,T,MAX_OSC,length(sigman));
MSE_osc = zeros(1,length(sigman));
for i=1:length(sigman)
    Y = c*x+sqrt(sigman(i))*randn(J,T);
    [osc_param(:,:,i),osc_AIC(:,i),osc_mean(:,:,:,i)] = osc_decomp(Y,fs,MAX_OSC);
    [~,Khat] = min(osc_AIC(:,i));
    osc_f = osc_param(Khat,Khat+1:2*Khat,i);
    [~,k] = min(abs(osc_f-f(1)));
    MSE_osc(i) = norm(osc_mean(2*k-1,:,Khat,i)-x(1,:))^2/T;
end

figure
plot(log10(sigman),MSE_osc,'k-','LineWidth',2)
set(gca,'FontSize',16);
xlabel('log10 \tau^2','FontSize',20);
ylabel('MSE','FontSize',20);
