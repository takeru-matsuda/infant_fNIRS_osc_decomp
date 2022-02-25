load('osc_summary.mat')

figure
plot(freq(a>0),a(a>0),'k+');
set(gca,'FontSize',16);
xlabel('frequency (Hz)','FontSize',20);
ylabel('a','FontSize',20);

figure
histogram((a_brain(a_brain>0)),0:.01:1)
set(gca,'FontSize',16);
xlabel('a','FontSize',20);
ylabel('count','FontSize',20);

figure
histogram((a_pulse(a_pulse>0)),0:.01:1)
set(gca,'FontSize',16);
xlabel('a','FontSize',20);
ylabel('count','FontSize',20);

figure
histogram((a_mirror(a_mirror>0)),0:.01:1)
set(gca,'FontSize',16);
xlabel('a','FontSize',20);
ylabel('count','FontSize',20);
