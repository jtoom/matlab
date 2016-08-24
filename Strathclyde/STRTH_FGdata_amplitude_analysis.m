% Strathclyde data 2014/07

clearvars
tic
%==========================================================================
% Settings

folder = '140708';
number = '4';
direction = 'down';
file_loc = ['\\10.48.24.77\Strathclyde\iDrive\data\' folder '\' number '\' direction '\'];    % directory where the files are
% file_loc = ['E:\Uni\Post Doc\Strathclyde\iDrive\data\Josh Analysis\' folder '\' number '\' direction '\with Random Perturbation\'];
file_list = ls([file_loc 'dpo*mA.dat']);
outputFolder = ['E:\Uni\Post Doc\Strathclyde\iDrive\data\Josh Analysis\' folder '\' number '\' direction '\'];
%==========================================================================

sz = size(file_list);

% Get injection range
for z = 1:sz(1)
    inj(z) = str2double(file_list(z,4:9));
end

% Initialise array
Arms = zeros(sz(1),1);
App = zeros(sz(1),1);
Amin = zeros(sz(1),1);
R = zeros(sz(1),1);

for a = 1:sz(1)    
    disp(['Current = ' num2str(inj(a),'%.2f') 'mA  ->  ' num2str((a/sz(1))*100,'%.2f') '% complete'])
    
    % Load time series    
    TS = load([file_loc file_list(a,:)]);
    TS = -1*TS;
    
    % Calculate the ratio between peak to peak and rms amplitude
    Arms(a) = std(TS);
    App(a) = max(TS) - min(TS);
    Amin(a) = min(TS);
    R(a) = App(a)/Arms(a);
end
toc
csvwrite([outputFolder 'Amp_analysis.csv'],[Arms App R]);

figure(1)
subplot(3,1,1)
plot(inj,R)
xlabel('Injection')
ylabel('pp / rms')
subplot(3,1,2)
plot(inj,Arms)
xlabel('Injection')
ylabel('rms')
subplot(3,1,3)
plot(inj,App)
xlabel('Injection')
ylabel('pp')

% print('-dpng','-r300',[outputFolder 'PE_map.png']);