
delay = h5read('E:\Uni\Post Doc\Strathclyde\iDrive\data\Josh Analysis\140709\5\down\ACF.h5','/delay');
current = h5read('E:\Uni\Post Doc\Strathclyde\iDrive\data\Josh Analysis\140709\5\down\ACF.h5','/current');
ACF = h5read('E:\Uni\Post Doc\Strathclyde\iDrive\data\Josh Analysis\140709\5\down\ACF.h5','/ACF');

delay_crop = h5read('E:\Uni\Post Doc\Strathclyde\iDrive\data\Josh Analysis\140709\5\down\ACF.h5','/delay',[1 1],[1 300]);
current_crop = h5read('E:\Uni\Post Doc\Strathclyde\iDrive\data\Josh Analysis\140709\5\down\ACF.h5','/current',[1 5],[1 225]);
ACF_crop = h5read('E:\Uni\Post Doc\Strathclyde\iDrive\data\Josh Analysis\140709\5\down\ACF.h5','/ACF',[1 5],[300 225]);

figure(71)
imagesc(delay_crop,current_crop,ACF_crop)
set(gca,'Ydir','normal')
xlabel('Delay (ns)')
ylabel('Current (mA)')
colormap(jet(256))
colorbar

% print('-dpng','-r300',['E:\Uni\Post Doc\Strathclyde\iDrive\data\Josh Analysis\' folder '\' number '\' direction '\ACF_crop.png'])

figure(72)
imagesc(delay,current,ACF)
set(gca,'Ydir','normal')
xlabel('Delay (ns)')
ylabel('Current (mA)')
xlim([0 300])
ylim([410 414.4])
colormap(jet(256))
colorbar
caxis([-1 1])

print('-dpng','-r300',['E:\Uni\Post Doc\Strathclyde\iDrive\data\Josh Analysis\' folder '\' number '\' direction '\ACF_crop.png'])