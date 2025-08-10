% custom_filter() - Filtering of EEG data
% Adapted from bemobil_filter() from the bemobil pipeline (Marius Klug)
%
% Inputs:
%   EEG                     - current EEGLAB EEG structure
%   lowcutoff               - low cut off frequency for firfilt filering, if [], no filter will be applied
%   highcutoff              - high cut of frequency, if [], no filter will be applied
%
% Outputs:
%   EEG                     - current EEGLAB EEG structure
%
%   .set data file of current EEGLAB EEG structure stored on disk (OPTIONALLY)

function [EEG] = custom_filter(EEG, lowcutoff, highcutoff)

if ~isempty(lowcutoff) || ~isempty(highcutoff)
    EEG.etc.filter.type = 'Hamming windowed sinc FIR filter (zero-phase)';
else
    error('No filter cutoffs specified, what was your plan here?!')
end

% Constants
TRANSWIDTHRATIO = 0.25;
fNyquist = EEG.srate/2;

%% High pass
if ~isempty(lowcutoff)
    if lowcutoff >= 0.5
        filter_order = 1650; % ensures a transition bandwith of 0.5Hz (see Klug&Gramann 2020)
        transition_bandwidth = 0.5; %in Hz, fixed by the filter order
        passband_edge = lowcutoff + transition_bandwidth/2;
        figure;
        [EEG, ~, ~, ~] = pop_eegfiltnew(EEG, passband_edge, 0, filter_order, 0, [], 1);
        EEG = eeg_checkset(EEG);
    else
        figure;
        [EEG, ~, ~, filter_order] = pop_eegfiltnew(EEG, lowcutoff, 0, [], 0, [], 1);
        EEG = eeg_checkset(EEG);
        % Default code from pop_eegfiltnew
        % to compute transition bandwidth when not specifying the order:
        passband_edge = lowcutoff;
        maxDf = passband_edge; % Band-/highpass
        transition_bandwidth = min([max([maxDf * TRANSWIDTHRATIO 2]) maxDf]);
        lowcutoff = passband_edge - transition_bandwidth/2;
    end
    
    disp(['Highpass filtered the data with ' num2str(lowcutoff) 'Hz cutoff, '...
        num2str(transition_bandwidth) 'Hz transition bandwidth, '...
        num2str(passband_edge) 'Hz passband edge, and '...
        num2str(filter_order) ' order.']);
    
    % removing and remaking the field is necessary for the order of the struct fields to be identical
    %if isfield(EEG.etc.filter,'highpass');  EEG.etc.filter = rmfield(EEG.etc.filter, 'highpass'); end
    EEG.etc.filter.highpass.cutoff = lowcutoff;
    EEG.etc.filter.highpass.transition_bandwidth = transition_bandwidth;
    EEG.etc.filter.highpass.passband = passband_edge;
    EEG.etc.filter.highpass.order = filter_order;
    close(gcf)
else
    if ~isfield(EEG.etc.filter,'highpass')
        EEG.etc.filter.highpass = 'not applied';
        %     else
        %         % removing and remaking the filed is necessary for the order of the struct fields to be identical
        %         temp = EEG.etc.filter.highpass;
        %         EEG.etc.filter = rmfield(EEG.etc.filter, 'highpass');
        %         EEG.etc.filter.highpass = temp;
    end
end

%% Low pass
if ~isempty(highcutoff)
    if highcutoff > fNyquist - 1
        disp('Warning: Cannot filter higher than Nyquist frequency.');
        highcutoff = fNyquist - 1;
        disp(['Now continuing with highest possible frequency: ' num2str(highcutoff)]);
    end
    
    %filter_order = 824; % ensures a transition bandwith of ~1Hz
    %transition_bandwidth = 1; %in Hz, fixed by the filter order
    %passband_edge = highcutoff - transition_bandwidth/2;
    
    figure;
    [EEG, ~, ~, filter_order] = pop_eegfiltnew(EEG, 0, highcutoff, [], 0, [], 1);
    EEG = eeg_checkset(EEG);
    
    % Default code from pop_eegfiltnew
    % to compute transition bandwidth when not specifying the order:
    passband_edge = highcutoff;
    maxDf = fNyquist - passband_edge; % Band-/highpass
    transition_bandwidth = min([max([passband_edge * TRANSWIDTHRATIO 2]) maxDf]);
    cutoff_freq = passband_edge + transition_bandwidth/2;
    
    disp(['Lowpass filtered the data with ' num2str(cutoff_freq) 'Hz cutoff, '...
        num2str(transition_bandwidth) 'Hz transition bandwidth, '...
        num2str(passband_edge) 'Hz passband edge, and '...
        num2str(filter_order) ' order.']);
    
    % removing and remaking the filter struct field is necessary for the order of the struct fields to be identical
    %if isfield(EEG.etc.filter,'lowpass'); EEG.etc.filter = rmfield(EEG.etc.filter, 'lowpass'); end
    EEG.etc.filter.lowpass.cutoff = cutoff_freq;
    EEG.etc.filter.lowpass.transition_bandwidth = transition_bandwidth;
    EEG.etc.filter.lowpass.passband = passband_edge;
    EEG.etc.filter.lowpass.order = filter_order;
    close(gcf)
else
    if ~isfield(EEG.etc.filter,'lowpass')
        EEG.etc.filter.lowpass = 'not applied';
        %     else
        %         % removing and remaking the filed is necessary for the order of the struct fields to be identical
        %         temp = EEG.etc.filter.lowpass;
        %         EEG.etc.filter = rmfield(EEG.etc.filter, 'lowpass');
        %         EEG.etc.filter.lowpass = temp;
    end
end
end
