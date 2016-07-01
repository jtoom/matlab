% Plot jitter distributions for Strathclyde paper
% Generate publication quality file
% **** Run STRTH_FGdata_pulse_stats.m first ****

fig88 = figure(88);
set(fig88,'Position',[50 50 600 800],'PaperPositionMode','auto')
fsp = subplot(3,1,1);
histogram(jitter_ALL/1e-12,'BinWidth',10);
%     xlabel('Difference from Mean Period (ps)')
set(get(gca,'child'),'FaceColor','b');
ylabel('Counts')
set(gca,'Xlim',[-80 80])
legend('All pulses')
LX = get(gca,'XLim');
LY = get(gca,'YLim');
str_ALL = ['StDev = ' num2str(std(jitter_ALL)/1e-12,'%.1f') '$$ \pm $$' num2str(unc_all/1e-12,'%.1f') 'ps'];
text(LX(1)+0.03*diff(LX),LY(1)+0.85*diff(LY),str_ALL, 'Interpreter','latex','FontSize',13)

subplot(3,1,2)
histogram(jitter_EVEN/1e-12,'BinWidth',10);
%     xlabel('Difference from Mean Period - Even Pulses (ps)')
set(get(gca,'child'),'FaceColor','g');
ylabel('Counts')
set(gca,'Xlim',[-80 80])
legend('Even pulses')
LX = get(gca,'XLim');
LY = get(gca,'YLim');
str_EVEN = ['StDev = ' num2str(std(jitter_EVEN)/1e-12,'%.1f') '$$ \pm $$' num2str(unc_EVEN/1e-12,'%.1f') 'ps'];
text(LX(1)+0.03*diff(LX),LY(1)+0.85*diff(LY),str_EVEN, 'Interpreter','latex','FontSize',13)

subplot(3,1,3)
histogram(jitter_ODD/1e-12,'BinWidth',10);
set(get(gca,'child'),'FaceColor','r');
xlabel('Difference from Mean Period (ps)')
ylabel('Counts')
set(gca,'Xlim',[-80 80])
legend('Odd pulses')
LX = get(gca,'XLim');
LY = get(gca,'YLim');
str_ODD = ['StDev = ' num2str(std(jitter_ODD)/1e-12,'%.1f') '$$ \pm $$' num2str(unc_ODD/1e-12,'%.1f') 'ps'];
text(LX(1)+0.03*diff(LX),LY(1)+0.85*diff(LY),str_ODD, 'Interpreter','latex','FontSize',13)

% print('-dpng','-r300',[analysis_loc 'Pulse Jitter_' folder ' - ' number ' - ' direction ' - ' current ' mA.png'])
% print('-depsc','-r300',[analysis_loc 'Pulse Jitter_' folder ' - ' number ' - ' direction ' - ' current ' mA.eps'])