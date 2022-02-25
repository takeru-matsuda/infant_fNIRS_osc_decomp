%% Figure 19
load oxydeoxy_corr.mat

ch18_x = [0.3827 0.2704 0.1837 0.1224 0.1020 0.1224 0.1837 0.2704 0.3827 0.5689 0.6811 0.7679 0.8291 0.8495 0.8291 0.7679 0.6811 0.5689];
ch18_y = [0.8444 0.7934 0.7117 0.5995 0.4796 0.3622 0.2500 0.1684 0.1173 0.1173 0.1684 0.2500 0.3622 0.4796 0.5995 0.7117 0.7934 0.8444];

figure('Position', [2050 25 1600 500]);
s = 100;
for mm = 1:2
    for kk = 1:7
    subplot(2,7,7*(mm-1)+kk)
    hold on
    axis square
    set(gca,'XTick',[],'YTick',[]);
    set(gca,'Visible','off');
    set(gcf,'Color','w');
    if mm == 1
        targetr = oxy_r;
    elseif mm == 2
        targetr = deoxy_r;
    end
    rthresh = 0.7-0.1*(kk-1);
    
    for ii = 1:18
        scatter(ch18_x(ii),ch18_y(ii),s,'k');
    end
    for ii = 1:17
        for jj = ii+1:18
            if targetr(ii,jj) > rthresh
                line([ch18_x(ii) ch18_x(jj)],[ch18_y(ii) ch18_y(jj)],'Color','r');
            end
        end
    end
    end
end

clear all
