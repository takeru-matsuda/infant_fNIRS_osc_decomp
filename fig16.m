load('osc_summary.mat')

figure

subplot(1,3,1)
histogram(freq_brain(freq_brain>0),.01:.01:.1)
xlim([.01 .1])
ylim([0 120])
set(gca,'FontSize',16)
xlabel('frequency','FontSize',20)
ylabel('count','FontSize',20);

subplot(1,3,2)
histogram(freq_pulse(freq_pulse>0),1.6:.1:2.4)
xlim([1.6 2.4])
ylim([0 120])
set(gca,'FontSize',16)
xlabel('frequency','FontSize',20)

subplot(1,3,3)
histogram(freq_mirror(freq_mirror>0),1.6:.1:2.4)
xlim([1.6 2.4])
ylim([0 120])
set(gca,'FontSize',16)
xlabel('frequency','FontSize',20)
