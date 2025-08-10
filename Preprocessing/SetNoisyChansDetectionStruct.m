function [ NoisyIn ] = SetNoisyChansDetectionStruct(EEG, globalArchi)

% May be useful to add the following path:
%addpath(genpath('/home/adelaux/Desktop/MATLAB_MT/eeglab-by-marius/eeglab14_1_0b/plugins/PrepPipeline0.5/utilities'));

all_chans = 1:EEG.nbchan;
EEG_chans = setdiff(all_chans, all_chans(strcmp({EEG.chanlocs.type},'EOG')));

NoisyIn=struct(...
    'srate', EEG.srate,... Sampling frequency in Hz
    'samples', size(EEG.data, 2),... number of samples in the data
    'chaninfo', EEG.chaninfo,... standard EEGLAB chaninfo (nose direction is relevant)
    'chanlocs', EEG.chanlocs,... standard EEGLAB chanlocs structure
    'evaluationChannels', EEG_chans... Channels to evaluate
    );

switch globalArchi
    case 'bemobil'
        disp('Using PREP pipeline...')
        
        %NoisyIn.evaluationChannels = 1:size(EEG.data, 1);
        NoisyIn.robustDeviationThreshold = 5; % z score cutoff for robust channel deviation (default = 5)
        NoisyIn.highFrequencyNoiseThreshold = 5; % z score cutoff for SNR (signal above 50 Hz) (default = 5)
        NoisyIn.correlationWindowSeconds = 1; % correlation window size in seconds (default = 1 sec)
        NoisyIn.correlationThreshold = 0.4; % correlation below which window is bad (default = 0.4)
        NoisyIn.badTimeThreshold = 0.01; % cutoff fraction of bad corr windows (default = 0.01)
        NoisyIn.ransacOff = false; % Whether to perform RANSAC or not (default = false)
        NoisyIn.ransacSampleSize = 50; % samples for computing ransac (default = 50)
        NoisyIn.ransacChannelFraction = 0.25; % fraction of channels for robust reconstruction (default = 0.25)
        NoisyIn.ransacCorrelationThreshold = 0.75; % cutoff correlation for abnormal wrt neighbors(default = 0.75)
        NoisyIn.ransacUnbrokenTime = 0.4; % cutoff fraction of time channel can have poor ransac predictability (default = 0.4)
        NoisyIn.ransacWindowSeconds = 5; % correlation window for ransac (default = 5 sec)
        
    case 'APP'
        disp('Using APP pipeline...')
        %NoisyIn.evaluationChannels = EEG_chans;
        NoisyIn.correlationTop = 4; % top best correlation coefficients (excluding self-correlation) that will be averaged
        % for bad correlation rejection
end
end

