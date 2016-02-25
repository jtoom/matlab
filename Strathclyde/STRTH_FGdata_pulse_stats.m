clear all

%==========================================================================
% Settings
folder = '140709';
number = '5';
direction = 'down';
current = '413.40';
testSR = 0;         % set to 1 to drop sample rate by half
addRAND = 1;

file_loc = ['\\10.48.16.125\Strathclyde\iDrive\data\' folder '\' number '\' direction '\'];    % directory where the files are
file_list = ls([file_loc 'dpo' current 'mA.dat']);    % you can enter the whole file name here if you just want to look at a single time series
analysis_loc = ['E:\Uni\Post Doc\Strathclyde\iDrive\data\Josh Analysis\' folder '\' number '\' direction '\'];
%==========================================================================

sz = size(file_list);

% Get injection range
for z = 1:sz(1)
    inj(z) = str2double(file_list(z,4:9));
end

% Initialise arrays
P2Pamp = zeros(sz(1),1);
VarTS = zeros(sz(1),1);

for a = 1:sz(1)  
    disp(inj(a))
    
    % Load time series (ASCII .dat)
    TS = load([file_loc file_list(a,:)]);
    TS = -1*TS;     % invert
    
    %=====================================
    % Add small random variation
    if addRAND == 1
        deltaTS = 0.000001;
        r = -deltaTS + (deltaTS+deltaTS).*rand(length(TS),1);
        TS = TS + r;
    end
    
    %=====================================
    
    % Test dropping the sampling rate by half (every second point)
    if testSR == 1
        TS = TS(1:2:end);
    end
    
%     % Load time series (binary .bin)
%     fileID = fopen([file_loc file_list(a,:)]);
%     TS = fread(fileID,'float32');
%     fclose(fileID);

    % Load OSA data
    fileID = fopen([file_loc 'fp_' current 'mA.dat']);
    C = textscan(fileID,'%s %s %f', 'HeaderLines',4);
    fclose(fileID);
    OSA = C{1,3};
    
    % Compute FFT
    xs = pow2(nextpow2(length(TS)));        % length of FFT
    if testSR == 1
        ts = 20e-12;
    else
        ts = 10e-12;                            % time step (100 GSa/s = 10 ps/pt)
    end
    f = (1:xs/2)./(ts*xs);                  % frequency scale
    [~,ind] = min(abs(f-25E9));             % find index of 25GHz
    f(ind:end) = [];
    [~,ind_L] = min(abs(f - 5E6));
    [~,ind_H] = min(abs(f - 9.5e9));
    fftx = fft(TS,xs);                  % computes FFT
    fftxr = sqrt(fftx(1:xs/2).*conj(fftx(1:xs/2)))/(xs/2);
    fftxr(ind:end) = [];                % removes frequency > 25GHz to save file size (scope BW = 23GHz)
    dBmx = 20*log10(fftxr/(.316));      % converts to dBm (assuming 50 ohm scope input)
    % Find dominant peak in FFT
    [fft_max,ind_max] = max(dBmx(ind_L:ind_H));
    ind_max = ind_max+(ind_L-1);
    % set the minimum distance between peaks (for peak detection) based on half the period of the peak FFT frequency
    mpd = (1/f(ind_max))/2;
    
    % Find Peaks
    avgTS = mean(TS);
    p2pTS = max(TS)-min(TS);
    [pks,loc,w,~] = findpeaks(TS,1/ts,'MinPeakHeight',min(TS)+0.5*p2pTS,'MinPeakDistance',mpd);

    % Statistics on peaks
    d = diff(loc);
    meanT = mean(d);    % Average period between pulses
    meanT_1stpair = mean(d(1:2:end));
    meanT_2ndpair = mean(d(2:2:end));
    stdT_1stpair = std(d(1:2:end));
    stdT_2ndpair = std(d(2:2:end));
    meanT_ODDpulse = mean(diff(loc(1:2:end)));    % Average period between every second pulse
    meanT_EVENpulse = mean(diff(loc(2:2:end)));    % Average period between every second pulse
    meanT_THIRDpulse = mean(diff(loc(1:3:end)));
    meanT_FOURTHpulse = mean(diff(loc(1:4:end)));
    stdT = std(diff(loc));      % standard deviation of period between pulses
    stdT_ODDpulse = std(diff(loc(1:2:end)));      % standard deviation of period between every 2nd pulse
    stdT_EVENpulse = std(diff(loc(2:2:end)));      % standard deviation of period between every 2nd pulse
    stdT_THIRDpulse = std(diff(loc(1:3:end)));
    stdT_FOURTHpulse = std(diff(loc(1:4:end)));
    meanW = mean(w);            % Average pulse width at half prominence
    meanW_ODDpulse = mean(w(1:2:end));
    meanW_EVENpulse = mean(w(2:2:end));
    stdW = std(w);              % Standard deviation of pulse width
    stdW_ODDpulse = std(w(1:2:end));
    stdW_EVENpulse = std(w(2:2:end));
    
    % Permutation Entropy (needs to be calculated prior to this)
    PE = load([analysis_loc 'PE_m3_dpo' current 'mA.csv']); 
    
    % Load ACF
    I = h5read([analysis_loc 'ACF.h5'],'/current');
    ind = find(I == inj);
    ACFun = h5read([analysis_loc 'ACF.h5'],'/ACF',[ind 1],[1 20000]);
    delay = h5read([analysis_loc 'ACF.h5'],'/delay');
end

if sz(1) > 1

else
    
    fig99 = figure(99);
    set(fig99,'Position',[30 30 1200 950],'PaperPositionMode','auto')
    subplot(3,3,1)
    plot((0:length(TS)-1)*ts,TS,loc,pks,'ro')
    xlabel('Time (s)')
    ylabel('Signal (V)')
    xlim([0 1e-8])
    
    subplot(3,3,2)
    plot(f./1e9,dBmx,f(ind_max)/1e9,fft_max,'ro')
    xlabel('Frequency (GHz)')
    ylabel('Power (dBm)')
    xlim([0 10])
    ylim([-110 -20])
    
    subplot(3,3,3)
    plot(f./1e9,dBmx)
    xlabel('Frequency (GHz)')
    ylabel('Power (dBm)')
    xlim([f(ind_max)/1e9-0.08 f(ind_max)/1e9+0.08])
    ylim([-110 -20])
    
    subplot(3,3,4)
    plot(OSA)
    xlabel('Optical Frequency (a.u.)')
    ylabel('Power (a.u.)')
    
    subplot(3,3,5)
    h_period = histogram(diff(loc)/1e-12,'BinWidth',10);
    xlabel('Pulse Period (ps)')
    ylabel('Counts')
    set(gca,'Xlim',[100 350])
    LX = get(gca,'XLim');
    LY = get(gca,'YLim');
    text(LX(1)+0.03*diff(LX),LY(1)+0.92*diff(LY),['Mean = ' num2str(meanT/1e-12,'%.1f') ' ps'])
    text(LX(1)+0.03*diff(LX),LY(1)+0.85*diff(LY),['StDev = ' num2str(stdT/1e-12,'%.1f') ' ps'])
    
    subplot(3,3,6)
    h_width = histogram(w/1e-12);
    xlabel('Pulse Width (ps)')
    ylabel('Counts')
%     set(gca,'Xlim',[70 160])
    LX = get(gca,'XLim');
    LY = get(gca,'YLim');
    text(LX(1)+0.03*diff(LX),LY(1)+0.92*diff(LY),['Mean = ' num2str(meanW/1e-12,'%.1f') ' ps'])
    text(LX(1)+0.03*diff(LX),LY(1)+0.85*diff(LY),['StDev = ' num2str(stdW/1e-12,'%.1f') ' ps'])
    
    subplot(3,3,7)
    plot(1:length(PE),PE)
    xlabel('Delay (pts)')
    ylabel('Permutation Entropy')
    
    subplot(3,3,8)
    plot(delay(1:200),ACFun(1:200))
    xlabel('Delay (pts)')
    ylabel('ACF')
    
    subplot(3,3,9)
    plot(delay(1:20000),ACFun(1:20000))
    xlabel('Delay (pts)')
    ylabel('ACF')
    
    if addRAND == 1
        print('-dpng','-r300',['E:\Uni\Post Doc\Strathclyde\iDrive\data\Josh Analysis\' folder '\' number '\' direction '\Pulse Stats ' current 'mA_with_random_pert.png'])
    else
        print('-dpng','-r300',['E:\Uni\Post Doc\Strathclyde\iDrive\data\Josh Analysis\' folder '\' number '\' direction '\Pulse Stats ' current 'mA.png'])
    end
    
%     figure(2)
%     plot(f./1e9,dBmx,f(ind_max)/1e9,fft_max,'ro')
%     xlabel('Frequency (GHz)')
%     ylabel('Power (dBm)')
%     xlim([0 10])
%     ylim([-110 -20])
%     
%     axes('Position',[0.6 0.65 0.28 0.25])
%     box on
%     plot(f./1e9,dBmx)
%     xlim([f(ind_max)/1e9-0.08 f(ind_max)/1e9+0.08])
%     ylim([-110 -20])
%     set(gca, 'YTick', []);
    
end
