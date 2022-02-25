%% Figure 18
load osc_corr.mat

ShortRange = [1 2; 2 3; 3 4; 4 5; 5 6; 6 7; 7 8; 8 9; 10 11; 11 12; 12 13; 13 14; 14 15; 15 16; 16 17; 17 18];
Contralateral = [1 18; 1 17; 2 18; 2 17; 2 16; 3 17; 3 16; 3 15; 4 16; 4 15; 4 14; 5 15; 5 14; 5 13; 6 14; 6 13; 6 12; 7 13; 7 12; 7 11; 8 12; 8 11; 8 10; 9 11; 9 10];
Ipsilateral = [1 9; 1 8; 2 9; 2 8; 2 7; 3 8; 3 7; 3 6; 4 7; 4 6; 10 17; 10 18; 11 16; 11 17; 11 18; 12 15; 12 16; 12 17; 13 15; 13 16];

brainoscShortRange = NaN(18,18,21);
brainoscContralateral = NaN(18,18,21);
brainoscIpsilateral = NaN(18,18,21);
brainoscControl = NaN(18,18,21);
pulseShortRange = NaN(18,18,21);
pulseContralateral = NaN(18,18,21);
pulseIpsilateral = NaN(18,18,21);
pulseControl = NaN(18,18,21);
mirrorShortRange = NaN(18,18,21);
mirrorContralateral = NaN(18,18,21);
mirrorIpsilateral = NaN(18,18,21);
mirrorControl = NaN(18,18,21);

counter = 0;
counter2 = 0;
for ii = 1:17
    for jj = ii+1:18
        counter = counter + 1;
        Allcombination(counter,:) = [ii,jj];
        if size(find(ShortRange(:,1) == ii & ShortRange(:,2) == jj),1) + size(find(Contralateral(:,1) == ii & Contralateral(:,2) == jj),1) + size(find(Ipsilateral(:,1) == ii & Ipsilateral(:,2) == jj),1) == 0
            counter2 = counter2 + 1;
            Control(counter2,:) = [ii,jj];
            brainoscControl(ii,jj,:) = brainosc_r(ii,jj,:);
            pulseControl(ii,jj,:) = pulse_r(ii,jj,:);            
            mirrorControl(ii,jj,:) = mirror_r(ii,jj,:);
        end
        if size(find(ShortRange(:,1) == ii & ShortRange(:,2) == jj),1) == 1
            brainoscShortRange(ii,jj,:) = brainosc_r(ii,jj,:);
            pulseShortRange(ii,jj,:) = pulse_r(ii,jj,:);            
            mirrorShortRange(ii,jj,:) = mirror_r(ii,jj,:);
        end
        if size(find(Contralateral(:,1) == ii & Contralateral(:,2) == jj),1) == 1
            brainoscContralateral(ii,jj,:) = brainosc_r(ii,jj,:);
            pulseContralateral(ii,jj,:) = pulse_r(ii,jj,:);            
            mirrorContralateral(ii,jj,:) = mirror_r(ii,jj,:);
        end
        if size(find(Ipsilateral(:,1) == ii & Ipsilateral(:,2) == jj),1) == 1
            brainoscIpsilateral(ii,jj,:) = brainosc_r(ii,jj,:);
            pulseIpsilateral(ii,jj,:) = pulse_r(ii,jj,:);            
            mirrorIpsilateral(ii,jj,:) = mirror_r(ii,jj,:);
        end        
    end
end

brainoscShortRangenanmean = squeeze(mean(mean(brainoscShortRange,2,'omitnan'),1,'omitnan'));
brainoscContralateralnanmean = squeeze(mean(mean(brainoscContralateral,2,'omitnan'),1,'omitnan'));
brainoscIpsilateralnanmean = squeeze(mean(mean(brainoscIpsilateral,2,'omitnan'),1,'omitnan'));
brainoscControlnanmean = squeeze(mean(mean(brainoscControl,2,'omitnan'),1,'omitnan'));
pulseShortRangenanmean = squeeze(mean(mean(pulseShortRange,2,'omitnan'),1,'omitnan'));
pulseContralateralnanmean = squeeze(mean(mean(pulseContralateral,2,'omitnan'),1,'omitnan'));
pulseIpsilateralnanmean = squeeze(mean(mean(pulseIpsilateral,2,'omitnan'),1,'omitnan'));
pulseControlnanmean = squeeze(mean(mean(pulseControl,2,'omitnan'),1,'omitnan'));
mirrorShortRangenanmean = squeeze(mean(mean(mirrorShortRange,2,'omitnan'),1,'omitnan'));
mirrorContralateralnanmean = squeeze(mean(mean(mirrorContralateral,2,'omitnan'),1,'omitnan'));
mirrorIpsilateralnanmean = squeeze(mean(mean(mirrorIpsilateral,2,'omitnan'),1,'omitnan'));
mirrorControlnanmean = squeeze(mean(mean(mirrorControl,2,'omitnan'),1,'omitnan'));

% Figure 18a
ch18_x = [0.3827 0.2704 0.1837 0.1224 0.1020 0.1224 0.1837 0.2704 0.3827 0.5689 0.6811 0.7679 0.8291 0.8495 0.8291 0.7679 0.6811 0.5689];
ch18_y = [0.8444 0.7934 0.7117 0.5995 0.4796 0.3622 0.2500 0.1684 0.1173 0.1173 0.1684 0.2500 0.3622 0.4796 0.5995 0.7117 0.7934 0.8444];

figure('Position', [1000 1560 1600 200]);
s = 100;

for kk = 1:4
    subplot(1,4,kk)
    hold on
    axis square
    set(gca,'XTick',[],'YTick',[]);
    set(gca,'Visible','off');
    set(gcf,'Color','w');
    if kk == 1
        targetcategory = ShortRange;
    elseif kk == 2
        targetcategory = Contralateral;
    elseif kk == 3
        targetcategory = Ipsilateral;
    else
        targetcategory = Control;
    end
    
    for ii = 1:18
        scatter(ch18_x(ii),ch18_y(ii),s,'k');
    end
    for ii = 1:17
        for jj = ii+1:18
            if size(find(targetcategory(:,1) == ii & targetcategory(:,2) == jj),1) == 1
                if kk == 1
                    line([ch18_x(ii) ch18_x(jj)],[ch18_y(ii) ch18_y(jj)],'Color','r');
                elseif kk == 2
                    line([ch18_x(ii) ch18_x(jj)],[ch18_y(ii) ch18_y(jj)],'Color','g');
                elseif kk == 3
                    line([ch18_x(ii) ch18_x(jj)],[ch18_y(ii) ch18_y(jj)],'Color','b');
                else
                    line([ch18_x(ii) ch18_x(jj)],[ch18_y(ii) ch18_y(jj)],'Color','k');
                end
            end
        end
    end
end

% Figure 18b
figure('Position',[1450 620 700 840])
hold on
scatter(ones(1,21),brainoscShortRangenanmean,'r','LineWidth',2)
scatter(1.1*ones(1,21),brainoscContralateralnanmean,'g','LineWidth',2)
scatter(1.2*ones(1,21),brainoscIpsilateralnanmean,'b','LineWidth',2)
scatter(1.3*ones(1,21),brainoscControlnanmean,'k','LineWidth',2)
scatter(2*ones(1,21),pulseShortRangenanmean,'r','LineWidth',2)
scatter(2.1*ones(1,21),pulseContralateralnanmean,'g','LineWidth',2)
scatter(2.2*ones(1,21),pulseIpsilateralnanmean,'b','LineWidth',2)
scatter(2.3*ones(1,21),pulseControlnanmean,'k','LineWidth',2)
scatter(3*ones(1,21),mirrorShortRangenanmean,'r','LineWidth',2)
scatter(3.1*ones(1,21),mirrorContralateralnanmean,'g','LineWidth',2)
scatter(3.2*ones(1,21),mirrorIpsilateralnanmean,'b','LineWidth',2)
scatter(3.3*ones(1,21),mirrorControlnanmean,'k','LineWidth',2)
xlim([0.5 3.9])
ax = gca;
ax.TickDir = "out";
ax.XTick = [1.15 2.15 3.15];
ax.XTickLabel = [{'brain oscillator'},{'pulse wave'},{'mirroring noise'}];
legend('Short range','Contralateral','Ipsilateral','Control','Location','southwest','FontSize', 12);
legend('boxoff')
ylabel('correlation coefficient')

clear all
