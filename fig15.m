load('AIC_results')

figure
AICdiff = AICtog-AICsep;
histogram(AICdiff(:))
xlim([-4000 0])
set(gca,'FontSize',16)
xlabel('AIC difference','FontSize',20)
ylabel('count','FontSize',20);

figure
AICdiff_control = AICtog_control-AICsep_control;
histogram(AICdiff_control(:))
xlim([-200 800])
set(gca,'FontSize',16)
xlabel('AIC difference','FontSize',20)
ylabel('count','FontSize',20);
