load('osc_summary.mat')

figure
polarhistogram(phase_diff(freq>.01 & freq <.1),18);

figure
plot((phase_diff(freq>.01 & freq <.1))/pi*180,proj_norm(freq>.01 & freq <.1),'k+')
xlim([0 360])
xticks(0:60:360)
ylim([0 2.5])
set(gca,'FontSize',16);
xlabel('phase diff.','FontSize',20);
ylabel('norm','FontSize',20);

figure
polarhistogram(phase_diff(freq>1.6 & freq <2.4),18);

figure
plot((phase_diff(freq>1.6 & freq <2.4))/pi*180,proj_norm(freq>1.6 & freq <2.4),'k+')
xlim([0 360])
xticks(0:60:360)
ylim([0 2.5])
set(gca,'FontSize',16);
xlabel('phase diff.','FontSize',20);
ylabel('norm','FontSize',20);

figure
histogram(log10(pow(freq>1.6 & freq<2.4 & abs(phase_diff-pi)>pi/2)),20)
xlim([-7 1])
set(gca,'FontSize',16)
xlabel('log10 power','FontSize',20)
ylabel('count','FontSize',20);

figure
histogram(log10(pow(freq>1.6 & freq<2.4 & abs(phase_diff-pi)<pi/2)),20)
xlim([-7 1])
set(gca,'FontSize',16)
xlabel('log10 power','FontSize',20)
ylabel('count','FontSize',20);
