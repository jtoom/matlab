% Run STRTH_FGdata_pulse_stats on time series first so the ACF and TS
% variables are loaded

% Plot spacetime plot
rndtrp = 71;
N = floor(length(TS)/rndtrp);
spaceTime = reshape(TS(1:rndtrp*N),[rndtrp,N]);

% Plot distribution of period between every second pulse

fig9 = figure(9);
set(fig9,'Position',[680 200 500 720],'PaperPositionMode','auto')
imagesc(spaceTime')
ylim([0 1000])
colormap(jet(256))
xlabel(['Round Trip Time (' num2str(rndtrp) 'pts)'])
ylabel('No. of Round Trips')



fig1 = figure(1);
set(fig1,'Position',[680 558 800 420],'PaperPositionMode','auto')
[acfPKS,acfLOCS,acfW,acfP] = findpeaks(ACFun(1:311));
findpeaks(ACFun(1:311),'Annotate','extents');
xlabel('Delay (pts)');
ylabel('ACF');
% xlim([0 300]);
templabels = ['0.5\tau';'1.0\tau';'1.5\tau';'2.0\tau';'2.5\tau';'3.0\tau';'3.5\tau';'4.0\tau'];
text(acfLOCS+2,acfPKS+0.1,templabels)


[acfPKS,acfLOCS,acfW,acfP] = findpeaks(ACFun);
acfODDS = acfW(1:2:end);
acfEVENS = acfW(2:2:end);

figure(2)
histogram(acfW,50);
xlabel('ACF peak width (pts)');
ylabel('Counts');
xlim([16.51 16.63]);
LX = get(gca,'XLim');
LY = get(gca,'YLim');
text(LX(1)+0.03*diff(LX),LY(1)+0.92*diff(LY),['Mean = ' num2str(mean(acfW),'%.4f') ' pts'])
text(LX(1)+0.03*diff(LX),LY(1)+0.85*diff(LY),['StDev = ' num2str(std(acfW),'%.4f') ' pts'])

figure(3)
histogram(acfODDS,50);
xlabel('ACF peak width (pts)');
ylabel('Counts');
xlim([16.51 16.63]);
LX = get(gca,'XLim');
LY = get(gca,'YLim');
text(LX(1)+0.03*diff(LX),LY(1)+0.92*diff(LY),['Mean = ' num2str(mean(acfODDS),'%.4f') ' pts'])
text(LX(1)+0.03*diff(LX),LY(1)+0.85*diff(LY),['StDev = ' num2str(std(acfODDS),'%.4f') ' pts'])

figure(4)
histogram(acfEVENS,50);
xlabel('ACF peak width (pts)');
ylabel('Counts');
xlim([16.51 16.63]);
LX = get(gca,'XLim');
LY = get(gca,'YLim');
text(LX(1)+0.03*diff(LX),LY(1)+0.92*diff(LY),['Mean = ' num2str(mean(acfEVENS),'%.4f') ' pts'])
text(LX(1)+0.03*diff(LX),LY(1)+0.85*diff(LY),['StDev = ' num2str(std(acfEVENS),'%.4f') ' pts'])