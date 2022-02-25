load('fNIRS_sample.mat')

imf = emd(Y(1,:));
figure
for i=1:size(imf,2)
    subplot(size(imf,2),1,i)
    plot((1:200)/fs,imf(1:200,i),'r','LineWidth',1)
    set(gca,'FontSize',12);
end
xlabel('Time (s)','FontSize',20)

imf = emd(Y(2,:));
figure
for i=1:size(imf,2)
    subplot(size(imf,2),1,i)
    plot((1:200)/fs,imf(1:200,i),'b','LineWidth',1)
    set(gca,'FontSize',12);
end
xlabel('Time (s)','FontSize',20)

