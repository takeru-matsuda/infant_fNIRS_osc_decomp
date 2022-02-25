load('osc_summary.mat')

figure
histogram(K(:))
xticks(3:7)
set(gca,'FontSize',16)
xlabel('number of oscillators (K)')
ylabel('count');

figure
histogram(freq(freq>0 & freq<.1))
xlim([0 .1])
set(gca,'FontSize',16)
xlabel('frequency (Hz)','FontSize',20)
ylabel('count','FontSize',20);

figure
histogram(freq(freq>0))
xlim([0 5])
set(gca,'FontSize',16)
xlabel('frequency (Hz)','FontSize',20)
ylabel('count','FontSize',20);

figure
plot(freq(:),log10(pow(:)),'k+');
xlim([0 0.1])
ylim([-6 4])
set(gca,'FontSize',16);
xlabel('frequency (Hz)','FontSize',20);
ylabel('log10 power','FontSize',20);

figure
plot(freq(:),log10(pow(:)),'k+');
ylim([-5 4])
set(gca,'FontSize',16);
xlabel('frequency (Hz)','FontSize',20);
ylabel('log10 power','FontSize',20);
