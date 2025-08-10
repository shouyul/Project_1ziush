function [ lineNoiseIn ] = SetLineNoiseRemovalStruct(EEG)

% May be useful to add the following path:
%addpath(genpath('/home/adelaux/Desktop/MATLAB_MT/eeglab-by-marius/eeglab14_1_0b/plugins/PrepPipeline0.5/utilities'));

duration = EEG.times(end)/1000; %duration in s
Nmax_4sec_wdws = floor(duration/4);

lineNoiseIn=struct('Fs', EEG.srate,... Sampling frequency
    'fPassBand', [0, EEG.srate/2],... Frequency band used (default [0, Fs/2] = entire band)
    'fScanBandWidth', 2,...  +/- bandwidth centered on each f0 to scan for significant lines (TM)
    'lineFrequencies', [50 100],... Line frequencies to be removed (default: multiples of 50Hz) + Here there is a huge artifact at 90Hz too, probably due to the VIVE
    'lineNoiseChannels', 1:size(EEG.data,1),... Channels to remove line noise from (default size(data, 1))
    'maximumIterations', 25,... Maximum times to iterate removal (default = 10)
    'p', 0.05,... Significance level cutoff (default = 0.01)
    'pad', 0,... FFT padding factor ( -1 corresponds to no padding, 0 corresponds to padding to next highest power of 2 etc.) (default is 0)
    'taperBandWidth',2,... Taper bandwidth (default 2 Hz)
    'taperWindowSize',duration/Nmax_4sec_wdws,... Taper sliding window length (default 4 sec)
    'taperWindowStep',duration/Nmax_4sec_wdws,... Sliding window step size (default 4 sec = no overlap)
    'tau',100 ... Window overlap smoothing factor (default 100)
    );
end

