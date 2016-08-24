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
m = 3;                  % Oridnal pattern length
tau_start = 1;          % Start delay value
tau_finish = 200;       % Finish delay value
NP = 20000;             % Number of time series points to use
%==========================================================================

sz = size(file_list);

% Get injection range
for z = 1:sz(1)
    inj(z) = str2double(file_list(z,4:9));
end

% Initialise array
PE = zeros(sz(1),tau_finish);

for a = 1:sz(1)    
    disp(['Current = ' num2str(inj(a),'%.2f') 'mA  ->  ' num2str((a/sz(1))*100,'%.2f') '% complete'])
    
    % Load time series
%     fileID = fopen([file_loc file_list(a,:)]);
%     TS = fread(fileID,'float32');
%     fclose(fileID);
    
    TS = load([file_loc file_list(a,:)]);
    TS = -1*TS(1:NP);
    
    % Compute PE as a function of delay (tau)
    [PET, ~, tau, ~] = PEcalc(TS,m,tau_start,tau_finish,1,0,1,0);
    PE(a,:) = PET';
end
toc
csvwrite([outputFolder 'PE_m' num2str(m) '_d_' num2str(tau_finish) '.csv'],PE);

figure(1)
imagesc(tau,inj,PE)
xlabel('Delay')
ylabel('PE')

% print('-dpng','-r300',[outputFolder 'PE_map.png']);