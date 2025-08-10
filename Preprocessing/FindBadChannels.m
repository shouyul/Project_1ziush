function [ noisyOut ] = FindBadChannels(EEG, cfg)
% Search for noisy channels
% Includes some preparatory processing steps to improve the performance of the detection:
% 1. HP filter the data
% 2. Remove Line Noise
% 3. Re-reference the channels in case of APP pipeline

subject = cfg.subjects(cfg.current_subject).id;

% filter
lowcutoff = cfg.filterPreProc.low_cut_off;
highcutoff = cfg.filterPreProc.high_cut_off;
if ~isempty(lowcutoff)
    fprintf('Highpass Filtering (%.1f Hz) for automatic bad channel detection...\n', lowcutoff)
end
if ~isempty(highcutoff)
    fprintf('Lowpass Filtering (%.1f Hz) for automatic bad channel detection...\n', highcutoff)
end
[EEG_HP] = custom_filter(EEG, lowcutoff, highcutoff);

% Add path to prepPipeline subdirectories if not in the list
tmp = which('getPipelineDefaults');
if isempty(tmp)
    myPath = fileparts(which('prepPipeline'));
    addpath(genpath(myPath));
end

% Remove Line Noise with PREP pipeline functions
disp('Removing Line Noise...')
[EEG_HP_noLN, lineNoiseOut] = removeLineNoise_custom(EEG_HP, cfg.lineNoiseRemoval_method, false);
clear EEG_HP
% If you want to save the filteredEEG_noLN struct with the LineNoiseRemoval information for later:
%EEG_HP_noLN.etc.lineNoiseRemoval = lineNoiseOut;

% Remove out-of-interest data segments
% disp('Cropping out-of-interest data...')
% [OoI_segments_index] = select_data_of_interest(EEG_HP_noLN, 'fulldata');
% EEG_ready4BadChans = eeg_eegrej(EEG_HP_noLN, OoI_segments_index);
% clear filteredEEG_noLN

% use prep pipeline subfunction to determine noisy channels
disp('Detecting bad channels automatically')
% Get the parameters (default for now)
noisyIn = SetNoisyChansDetectionStruct(EEG_HP_noLN, cfg.globalArchitecture);
if strcmp(subject,'P1002') && contains(cfg.study_folder, 'Tumbler')
    noisyIn.correlationThreshold = 0.25;
    noisyIn.ransacOff = true;
end

switch cfg.globalArchitecture
    case 'bemobil'
        %% Check the figures folder exists
        if ~exist(fullfile(cfg.figures_folder, 'PREP_distributions'), 'dir')
            mkdir(fullfile(cfg.figures_folder, 'PREP_distributions'));
        end
        
        %% Check for bridges
        if ~(strcmp(subject, 'P1001-2') || strcmp(subject, 'P1002'))
            bridgesIn = rmfield(noisyIn, {'robustDeviationThreshold','highFrequencyNoiseThreshold',...
                'correlationWindowSeconds', 'correlationThreshold', 'badTimeThreshold',...
                'ransacSampleSize', 'ransacChannelFraction', 'ransacCorrelationThreshold', 'ransacUnbrokenTime', 'ransacWindowSeconds'});
            bridgesIn.epochWindowSeconds = 1;
            %bridgesIn.correlationThreshold = 0.99;
            %bridgesIn.badTimeThreshold = 0.95;
            if isfile(fullfile('utils', sprintf('%s_grid_display.mat', cfg.capName)))
                bridgesIn.gridPattern = load(fullfile('utils', sprintf('%s_grid_display.mat', cfg.capName)));
            else
                error('You need to create a grid pattern file for this cap to check for bridged channels.')
            end
            % Check the grid in case some channels were evicted
            chansInGrid = false(size(bridgesIn.gridPattern.grid_labels));
            for ch = 1:numel(EEG.chanlocs)
                chansInGrid = chansInGrid | strcmp(bridgesIn.gridPattern.grid_labels, EEG.chanlocs(ch).labels);
            end
            
            for x = 1:size(bridgesIn.gridPattern.grid_labels,1)
                for y = 1:size(bridgesIn.gridPattern.grid_labels,2)
                    if chansInGrid(x,y)
                        bridgesIn.gridPattern.grid_inds(x,y) =...
                            find(strcmp({EEG.chanlocs.labels}, bridgesIn.gridPattern.grid_labels{x,y}));
                    else
                        bridgesIn.gridPattern.grid_labels{x,y} = [];
                        bridgesIn.gridPattern.grid_inds(x,y) = 0;
                    end
                end
            end
            bridgesOut = findBridgedChannels(EEG_HP_noLN, bridgesIn, false);
            
            % Plot Low ED results
            if ~isempty(bridgesOut.lowElecDistClusters)
                channelsBycluster = zeros(1,EEG_HP_noLN.nbchan);
                for cl = 1:numel(bridgesOut.lowElecDistClusters)
                    channelsBycluster(bridgesOut.lowElecDistClusters{cl}) = 10+cl;
                end
                
                figure
                topoplot(channelsBycluster, EEG_HP_noLN.chanlocs,...
                    'plotgrid', bridgesOut.gridPattern.grid_inds, 'maplimits', 'absmax');
                title({sprintf('%s - Clusters of low electrical distance between channels', subject),...
                    sprintf('Electrical distance below %.3f microVolts for at least half of the recording',...
                    sqrt(bridgesOut.elecDistCutoff)),...
                    sprintf('Found %d clusters',numel(bridgesOut.lowElecDistClusters)),...
                    'Green = not connected'});
                saveCurrentFig([cfg.figures_folder 'PREP_distributions' filesep],...
                    [subject '_bridges_byElecDist'], {'png'}, [600 500]);
                
                bridgedElectrodes = [];
                coords = [[EEG_HP_noLN.chanlocs.X];[EEG_HP_noLN.chanlocs.Y];[EEG_HP_noLN.chanlocs.Z]];
                for cl = 1:numel(bridgesOut.lowElecDistClusters)
                    % Find which channel is closest to the cluster centroid
                    centroid = mean(coords(:,bridgesOut.lowElecDistClusters{cl}),2);
                    [~, ch_ind] = min(vecnorm(coords - repmat(centroid,1,EEG_HP_noLN.nbchan),2,1));
                    % Flag all electrodes in the cluster except the one we
                    % decided to keep
                    bridgedElectrodes = union(bridgedElectrodes, setdiff(bridgesOut.lowElecDistClusters{cl}, ch_ind));
                end
                noisyIn.evaluationChannels = setdiff(noisyIn.evaluationChannels, bridgedElectrodes);
            end
            
            % Plot High correlation results to compare
            if ~isempty(bridgesOut.highCorrClusters)
                channelsBycluster = zeros(1,EEG_HP_noLN.nbchan);
                for cl = 1:numel(bridgesOut.highCorrClusters)
                    channelsBycluster(bridgesOut.highCorrClusters{cl}) = 10+cl;
                end
                
                figure
                topoplot(channelsBycluster, EEG_HP_noLN.chanlocs,...
                    'plotgrid', bridgesOut.gridPattern.grid_inds, 'maplimits', 'absmax');
                title({sprintf('%s - Clusters of highly correlated channels', subject),...
                    sprintf('Correlation above %.3f for at least half of the recording',...
                    bridgesOut.corrCutoff),...
                    sprintf('Found %d clusters',numel(bridgesOut.highCorrClusters)),...
                    'Green = not correlated'});
                saveCurrentFig([cfg.figures_folder 'PREP_distributions' filesep],...
                    [subject '_bridges_byCorrelation'], {'png'}, [600 500]);
            end
        else
            warning(['NOT CHECKING FOR BRIDGES FOR SUJBECT ' subject ' (see comment)'])
            % Bridge check ended in error when executing findBridgedChannels (below), on line 79.
            % Seems it is because EEG data from amp A (that was retrieved from eego rec) has very high amplitudes (~10^6-10)
            % at some time points (at least beginning, other times also), and makes electrical distance too high for the script to handle.
            % Not sure it's the reason.
        end
        
        %% PREP detection of bad channels
        noisyOut = findNoisyChannels(EEG_HP_noLN, noisyIn);
        
        %% Add bridged channels
        if exist('bridgedElectrodes', 'var')
            noisyOut.noisyChannels.bridgedChannels = unique(bridgedElectrodes); % unique sorts the channels too
            noisyOut.noisyChannels.all = union(noisyOut.noisyChannels.all, bridgedElectrodes);
            noisyOut.noisyChannels.all = unique(noisyOut.noisyChannels.all); % unique sorts the channels too
        end
        
        %% Add bad channels flagged by user
        if ~isempty(cfg.subjects(cfg.current_subject).badElectrodes)
            badElectrodes = cfg.subjects(cfg.current_subject).badElectrodes;
            for ch = 1:numel(badElectrodes)
                i = find(strcmp({noisyOut.channelLocations.labels}, badElectrodes{ch}));
                if find(noisyOut.noisyChannels.badChannelsFromNoData == i)
                    fprintf('%s already flagged as bad by PREP pipeline\n', badElectrodes{ch})
                else
                    fprintf('Adding %s to NoData channels\n', badElectrodes{ch})
                    noisyOut.noisyChannels.badChannelsFromNoData(end+1) = i;
                    noisyOut.noisyChannels.badChannelsFromNoData = sort(noisyOut.noisyChannels.badChannelsFromNoData);
                    noisyOut.noisyChannels.all(end+1) = i;
                    noisyOut.noisyChannels.all = unique(noisyOut.noisyChannels.all); % unique sorts the channels too
                end
            end
        end
        
        % Do some plots:
        plot_distribution(noisyOut.robustChannelDeviation, noisyOut.noisyChannels.badChannelsFromDeviation,...
            'Channel', 'PREP', noisyOut.robustDeviationThreshold)
        xlabel(['Robust zscore for robust standard deviation'])
        title({[subject ' - Distribution of the deviation criterion'],...
            'for the channels inspected by PREP'})
        saveCurrentFig([cfg.figures_folder 'PREP_distributions' filesep],...
            [subject '_channels_deviation'], {'png'}, [600 500]);
        
        plot_distribution(noisyOut.zscoreHFNoise, noisyOut.noisyChannels.badChannelsFromHFNoise,...
            'Channel', 'PREP', noisyOut.highFrequencyNoiseThreshold)
        xlabel(['Robust zscore for robust estimate of the power ratio HF/LF'])
        title({[subject ' - Distribution of the HF-noisiness criterion'],...
            'for the channels inspected by PREP'})
        saveCurrentFig([cfg.figures_folder 'PREP_distributions' filesep],...
            [subject '_channels_HF_noisy'], {'png'}, [600 500]);
        
        plot_distribution(noisyOut.medianMaxCorrelation, noisyOut.noisyChannels.badChannelsFromCorrelation,...
            'Channel', 'PREP')
        xlabel(['Median maximum correlation over ' num2str(noisyOut.correlationWindowSeconds) 's windows'])
        title({[subject ' - Distribution of the correlation criterion'],...
            'for the channels inspected by PREP'})
        saveCurrentFig([cfg.figures_folder 'PREP_distributions' filesep],...
            [subject '_channels_correlation'], {'png'}, [600 500]);
        
        plot_distribution(noisyOut.ransacBadWindowFraction, noisyOut.noisyChannels.badChannelsFromRansac,...
            'Channel', 'PREP', noisyOut.ransacUnbrokenTime)
        xlabel(['Fraction of windows (' num2str(noisyOut.ransacWindowSeconds) 's) failing to correlate with the RANSAC prediction'])
        title({[subject ' - Distribution of the RANSAC criterion'],...
            'for the channels inspected by PREP'})
        saveCurrentFig([cfg.figures_folder 'PREP_distributions' filesep],...
            [subject '_channels_RANSAC'], {'png'}, [600 500]);
end
end

