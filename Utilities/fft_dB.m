function [dBmx,f] = fft_dB(TS,dt)
% Calculate the FFT and convert to dBm based on 50ohm scope load

% Inputs:
%   TS = time series on which to calculate FFT
%   dt = time series sampling interval (e.g. 50 GSample/sec -> dt = 20ps)
% Outputs:
%   dBmx = FFT magnitude in dBm
%   f = frequency values vector

    xs = length(TS);        % length of time series
    f = (1:xs/2)./(dt*xs);  % frequency scale

    fftx = fft(TS,xs);      %computes FFT
    fftxr = sqrt(fftx(1:xs/2).*conj(fftx(1:xs/2)))/(xs/2);
    dBmx = 20*log10(fftxr/(.316));    % converts to dBm