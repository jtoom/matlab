% Run STRTH_FGdata_pulse_stats first to get variables loaded

% This script should highlight which pulses have widths that make up the 2
% peaks seen in some distributions of peak width

cut_off = 300e-12;

log_ind = w < cut_off;

test_loc_big = loc;
test_loc_small = loc;
test_pks_big = pks;
test_pks_small = pks;

test_pks_big(log_ind) = [];
test_loc_big(log_ind) = [];

test_pks_small(~log_ind) = [];
test_loc_small(~log_ind) = [];

figure(1)
plot((0:length(TS)-1)*ts,TS,test_loc_big,test_pks_big,'ro',test_loc_small,test_pks_small,'go')
% xlim([6.75e-6 6.81e-6])
xlabel('Time (s)')
ylabel('Signal (V)')

figure(2)
findpeaks(TS,1/ts,'MinPeakHeight',min(TS)+0.5*p2pTS,'MinPeakDistance',mpd,'Annotate','extents');
% xlim([6.75e-6 6.81e-6])
% ylim([-0.02 0.04])
xlabel('Time (s)')
ylabel('Signal (V)')

% print('-dpng','-r300',['E:\Uni\Post Doc\Strathclyde\iDrive\data\Josh Analysis\' folder '\' number '\' direction '\Pulse Stats ' current 'mA_incorrect_width.png'])